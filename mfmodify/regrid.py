# IMPORT
import os
import inspect
import numpy as np
import pandas as pd
import flopy
from .utils import (
    get_parameter_set, 
    param_dict_from_list, 
    get_ts_objects, 
    get_obs_objects
)

# VARIABLES
pack_type_scale_col = {
    'drn': 'cond',
    'wel': 'q',
    'riv': 'cond',
    'ghb': 'cond',
    'rch': None,
    'evt': None,
    'chd': None
}

# FUNCTIONS

def regrid_package_sp_data(pack, convert_df, scale_col):
    """
    Regrid stress period data for a given package.

    Parameters:
    pack (object): The package containing stress period data.
    convert_df (pd.DataFrame): DataFrame used to convert the original grid to 
        the new grid. It should contain a column 'area_factor' to scale the data.
    scale_cols (list of str): List of column names in the stress period data 
        that need to be scaled by 'area_factor'.

    Returns:
    dict: A dictionary where keys are stress period identifiers and values are 
        the regridded stress period data as recarrays.
    """
    # new dictionary to hold converted sp data
    new_sp_data = {}
    # loop over all stress periods
    for sp, sp_rec in pack.stress_period_data.data.items():
        # get orig columns
        orig_sp_cols = list(sp_rec.dtype.names)
        # get orig sp data as dataframe
        orig_sp_df = (
            pd
            .DataFrame(sp_rec)
            .assign(original_order = lambda x: x.index)
            .set_index('cellid')
        )
        # associate with new grid by joining with convert_df
        new_sp_df = (
            convert_df
            .join(orig_sp_df)
            .dropna()
            .sort_values('original_order')
        )
        # apply area factor
        if (scale_col == 'multiplier') and ('multiplier' not in new_sp_df.columns):
                new_sp_df['multiplier'] = 1
                orig_sp_cols.append('multiplier')
        new_sp_df[scale_col] = new_sp_df[scale_col] * new_sp_df['area_factor']
        # reorder columns
        new_sp_df = new_sp_df.loc[:, orig_sp_cols]
        # convert to recarray
        new_sp_rec = new_sp_df.to_records(index=False)
        # assign new recarray for current sp
        new_sp_data[sp] = new_sp_rec
    return new_sp_data

def convert_array_dis_to_disv(dis_array, convert_df):
    """
    Convert a DIS array to a DISV array using a conversion DataFrame.

    Parameters:
    dis_array (np.ndarray): 3D array representing the DIS grid.
    convert_df (pd.DataFrame): DataFrame used to convert the DIS grid to the DISV 
        grid. It should contain columns 'cellid' and 'area_factor'.

    Returns:
    np.ndarray: 3D array representing the DISV grid.
    """
    # loop over values to make a dataframe
    cellid_dis = []
    value_list = []
    nlay, nrow, ncol = dis_array.shape
    for ilay in range(nlay):
        for irow in range(nrow):
            for icol in range(ncol): 
                value_list.append(dis_array[ilay, irow, icol])
                cellid_dis.append((ilay, irow, icol))
    # make the dataframe
    val_df_orig = (
        pd
        .DataFrame({'cellid_dis': cellid_dis, 'value': value_list})
        .set_index('cellid_dis')
    )
    val_df_disv = (
        convert_df
        .join(val_df_orig)
        .assign(layer = lambda x: [int(cellid[0]) for cellid in x.cellid])
        .assign(node = lambda x: [int(cellid[1]) for cellid in x.cellid])
    )
    disv_array = np.stack([
        x[1].value.values for x in val_df_disv.groupby('layer')
    ])
    return disv_array

def make_regridded_package_param_dict(pack, convert_df):
    # get list of attributes to use as parameters in instantiating object
    rem_att_set = set(['loading_package', 'stress_period_data'])
    param_set = get_parameter_set(pack, rem_att_set=rem_att_set)
    # find array data
    array_class = flopy.mf6.data.mfdataarray.MFArray
    array_params = set([
        att for att in param_set if (isinstance(getattr(pack, att), array_class)
        and getattr(pack, att).has_data())
    ])
    param_list = list(param_set - set(array_params))
    # create parameter dictionary 
    # get those from attribute list
    pack_param_dict = param_dict_from_list(pack, param_list)
    # add others manually
    pack_param_dict['pname'] = pack.package_name
    pack_param_dict['filename'] = pack.filename
    # convert grid-dependent attributes that may be present
    # array data
    nlay = convert_df.sort_values('cellid').cellid.values[-1][0] + 1
    if len(array_params) > 0:
        for att in array_params:
            orig_value = getattr(pack, att)
            file_entry = orig_value.get_file_entry()
            if 'LAYERED' in file_entry:
                if 'CONSTANT' in file_entry:
                    pack_param_dict[att] = [orig_value[x].data_const_value[0] for x in range(nlay)]
                else:
                    pack_param_dict[att] = convert_array_dis_to_disv(orig_value.data, convert_df)
            else:
                value = float(file_entry.split('CONSTANT')[1].strip())
                pack_param_dict[att] = value
    # stress-period data
    if pack.has_stress_period_data:
        # get scale column from dictionary
        sp_scale_col = pack_type_scale_col[pack.package_type]
        # check if scale_col has numeric type
        sp_data_df0 = pack.stress_period_data.dataframe[0]
        if not pd.api.types.is_numeric_dtype(sp_data_df0[sp_scale_col].dtype):
            # if not numeric, change scale col to multiplier
            sp_scale_col = 'multiplier'
            # if there isn't already a multiplier, add to the parameters
            if ('auxmultname' not in pack_param_dict) or (pack_param_dict['auxmultname'] != 'multiplier'):
                pack_param_dict['auxmultname'] = 'multiplier'
                pack_param_dict['auxiliary'] = ['multiplier']
        sp_data_new = regrid_package_sp_data(pack, convert_df, sp_scale_col)
        pack_param_dict['stress_period_data'] = sp_data_new
    return pack_param_dict

def refine_package(sim_or_gwf_orig, pack_name, sim_or_gwf_new, convert_df, manual_params={}):
    # get package
    pack_orig = sim_or_gwf_orig.get_package(pack_name)
    # get a variable to point to the class of the package
    pack_class = pack_orig.__class__

    # get package paramters
    pack_param_dict = make_regridded_package_param_dict(
        pack_orig, convert_df
    )
    for att, val in manual_params.items():
        pack_param_dict[att] = val
    # get child packages
    child_packs = pack_orig.get_package().copy()
    # instantiate package
    pack_new = pack_class(sim_or_gwf_new, **pack_param_dict)
    # convert and reassociate timeseries packages
    if hasattr(pack_orig, 'ts'):
        ts_objects = get_ts_objects(pack_orig)
        if len(ts_objects) > 0:
            i = 0
            for ts_obj in ts_objects:
                # get the package param dict
                ts_dict = copy_param_dict(ts_obj)
                if i==0:
                    ts_dict['pname'] = f'{pack_name}_ts'
                    pack_new.ts.initialize(**ts_dict)
                else:
                    ts_dict['pname'] = f'{pack_name}_ts{i}'
                    pack_new.ts.append_package(**ts_dict)
                i+=1
    # TODO: screen for hdobs packages and change cellids
    # reassociate observation packages
    if hasattr(pack_orig, 'obs'):
        obs_objects = get_obs_objects(pack_orig)
        if len(obs_objects) > 0:
            i = 0
            for obs_obj in obs_objects:
                # get the package param dict
                obs_dict = copy_param_dict(obs_obj)
                # rename and add to parent package
                if i==0:
                    obs_dict['pname'] = f'{pack_name}_obs'
                    pack_new.obs.initialize(**obs_dict)
                else:
                    obs_dict['pname'] = f'{pack_name}_obs{i}'
                    pack_new.obs.append_package(**obs_dict)
                i+=1
    return pack_new

def modify_gridgen_disv_props(disv_props_orig):
    # Make a disv package
    # fix the vertices by removing duplicate vertices
    orig_vertex_df = (
        pd
        .DataFrame(disv_props_orig['vertices'], columns=['orig_vertex', 'x', 'y'])
        .assign(x = lambda xx: xx.x.round(4))
        .assign(y = lambda xx: xx.y.round(4))
    )
    new_vertex_df = (
        orig_vertex_df
        .drop_duplicates(subset=['x', 'y'], keep='first')
        .reset_index(drop=True)
        .assign(orig_vertex = lambda x: x.index)
        .rename(columns={'orig_vertex': 'vertex'})
    )
    vertex_relate_df = (
        orig_vertex_df
        .set_index(['x', 'y'])
        .join(new_vertex_df.set_index(['x', 'y']))
        .reset_index(drop=True)
        .sort_values('orig_vertex')
    )
    vertex_lut = dict(zip(vertex_relate_df.orig_vertex, vertex_relate_df.vertex))
    vertices_new = list(list(x) for x in new_vertex_df.itertuples(index=False, name=None))
    # change the cell definitions by replacing the old vertex numbers with the new ones
    cell2d_orig = disv_props_orig['cell2d']
    cell2d_new = []
    for i_cell_info in cell2d_orig:
        node_info = i_cell_info[:4]
        vertices = i_cell_info[4:]
        # loop through vertices and remove any that are redundant
        new_vertices = []
        for vertex in vertices:
            new_vertex = vertex_lut[vertex]
            if new_vertex not in new_vertices:
                new_vertices.append(new_vertex)
        node_info[-1] = len(new_vertices)
        cell2d_new.append(node_info + new_vertices)
    # get the new number of vertices
    nvert_new = len(vertices_new)
    disv_props = disv_props_orig.copy()
    disv_props['vertices'] = vertices_new
    disv_props['nvert'] = nvert_new
    disv_props['cell2d'] = cell2d_new
    return disv_props