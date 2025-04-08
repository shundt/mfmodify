# IMPORT
import shutil
import numpy as np
import pandas as pd
import shapely
import geopandas as gpd
import flopy
from flopy.utils.gridgen import Gridgen
from .utils import (
    get_parameter_set, 
    param_dict_from_list, 
    get_ts_objects, 
    get_obs_objects,
    copy_param_dict,
    copy_empty_sim,
    copy_package
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
    """
    Create a dictionary of parameters for a regridded package.

    This function extracts parameters from an existing package, processes 
    grid-dependent attributes (e.g., arrays and stress-period data), and 
    converts them to match a new grid structure defined by a conversion 
    DataFrame.

    Parameters:
    pack (object): The original package object from which parameters are extracted.
    convert_df (pd.DataFrame): A DataFrame used to convert the original grid 
        to the new grid. It should contain columns such as 'cellid' and 
        'area_factor' for mapping and scaling.

    Returns:
    dict: A dictionary containing the parameters required to instantiate 
        a regridded package. This includes converted grid-dependent attributes 
        and stress-period data, as well as other package-specific parameters.
    """
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
                    pack_param_dict[att] = [
                        orig_value[x].data_const_value[0] for x in range(nlay)]
                else:
                    pack_param_dict[att] = convert_array_dis_to_disv(
                        orig_value.data, convert_df)
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

def refine_package(sim_or_gwf_orig, pack_name, sim_or_gwf_new, convert_df, 
    manual_params={}):
    """
    Refine and regrid a package from an original simulation or model to a new one.

    This function takes a package from an original simulation or groundwater flow 
    model (GWF), regrids its parameters to match a new grid structure, and 
    associates it with a new simulation or model. It also handles timeseries 
    and observation packages associated with the original package.

    Parameters:
    sim_or_gwf_orig (flopy.mf6.MFSimulation or flopy.mf6.ModflowGwf): The original 
        simulation or GWF model containing the package to be refined.
    pack_name (str): The name of the package to be refined.
    sim_or_gwf_new (flopy.mf6.MFSimulation or flopy.mf6.ModflowGwf): The new 
        simulation or GWF model to which the refined package will be added.
    convert_df (pd.DataFrame): A DataFrame used to convert the original grid 
        to the new grid. It should contain columns such as 'cellid' and 
        'area_factor' for mapping and scaling.
    manual_params (dict, optional): A dictionary of additional parameters to 
        override or add to the refined package. Defaults to an empty dictionary.

    Returns:
    flopy.mf6.ModflowGwfPackage: The refined package associated with the new 
        simulation or model.
    """
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
    """
    Modify DISV properties generated by Gridgen to remove duplicate vertices.

    This function processes the DISV properties dictionary created by Gridgen, 
    removes duplicate vertices, and updates the cell definitions to use the 
    new vertex indices. The resulting DISV properties are suitable for use 
    in a MODFLOW 6 DISV package.

    Parameters:
    disv_props_orig (dict): A dictionary containing the original DISV properties 
        generated by Gridgen. It should include keys such as 'vertices' and 'cell2d'.

    Returns:
    dict: A modified dictionary of DISV properties with duplicate vertices removed 
        and updated cell definitions.
    """
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

def refine_gwf_dis_to_disv(sim_orig, model_name, grid_relate, disv_props, sim_ws_new):
    """
    Refine a GWF model from a structured DIS grid to an unstructured DISV grid.

    This function converts a groundwater flow (GWF) model using a structured DIS 
    grid to an unstructured DISV grid. It handles the conversion of grid-related 
    properties, packages, and boundary conditions, ensuring compatibility with 
    the new grid structure.

    Parameters:
    sim_orig (flopy.mf6.MFSimulation): The original simulation containing the GWF model.
    model_name (str): The name of the GWF model to be refined.
    grid_relate (pd.DataFrame): A DataFrame mapping the original grid to the new grid. 
        It should include columns such as 'cellid_disv' and 'area_factor'.
    disv_props (dict): A dictionary of DISV properties, including vertices and cell2d 
        definitions, for the new grid.
    sim_ws_new (str): The working directory for the new simulation.

    Returns:
    flopy.mf6.MFSimulation: The refined simulation with the GWF model converted to DISV.
    """
    # Get original model
    gwf_orig = sim_orig.get_model(model_name)
    
    # Check if all packages in the gwf model can be converted by this function
    valid_gwf_package_set = set([
        'dis', 'sto', 'npf', 'ic', 'oc', 'chd', 'wel', 'drn', 'riv', 'ghb', 
        'obs', 'oc'])
    gwf_package_df = (
        pd
        .DataFrame(gwf_orig.name_file.packages.array)
        .assign(ftype = lambda x: [ft[:-1].lower() for ft in x.ftype])
    )
    
    gwf_package_set = set(gwf_package_df.ftype)
    unhandled_paks = gwf_package_set.difference(valid_gwf_package_set)
    if len(unhandled_paks)>0:
        raise ValueError(
            f'The gwf package types {unhandled_paks} are not yet handled')
    
    # Make conversion dataframe
    convert_df = (
        grid_relate
        .rename(columns = {'cellid_disv': 'cellid', 'prop_dis_area': 'area_factor'})
        .loc[: ,['cellid', 'area_factor']]
    )
    
    # Make a new model
    # (MFSIM) Simulation --
    sim_new = copy_empty_sim(sim_orig, sim_ws_new)
    
    # (TDIS) Temporal discretization --
    tdis_new = copy_package(sim_orig, 'tdis', sim_new)
    
    # (GWF) Groundwater flow model
    gwf_new = flopy.mf6.ModflowGwf(sim_new, modelname=model_name)
    
    # (IMS) Iterative model solution
    ims_new = copy_package(sim_orig, 'ims', sim_new)
    sim_new.register_ims_package(ims_new, [model_name])
    
    # (DISV) Discretization by vertices
    disv_props_new = modify_gridgen_disv_props(disv_props)
    # get original dis object
    dis_orig = gwf_orig.get_package('dis')
    # convert idomain, to new grid
    idomain_disv = convert_array_dis_to_disv(dis_orig.idomain.data, convert_df)
    disv_props_new['idomain'] = idomain_disv
    # add other properties from dis
    # for prop in ['length_units', 'xorigin', 'yorigin', 'angrot']:
        # disv_props[prop] = getattr(dis_orig, prop).data
    disv_props_new['length_units'] = getattr(dis_orig, 'length_units').data
    # remove dis grid and add disv grid
    disv = flopy.mf6.modflow.ModflowGwfdisv(gwf_new, **disv_props_new) 
    
    # Node properties (no multiplier, just assign original values)
    # (STO) Storage
    if gwf_orig.get_package('sto') != None:
        # get steady-state info from gwf
        steady_state = gwf_orig.modeltime.steady_state
        # make dicts and add to manual param dict
        steady_state_dict = {i:bool(val) for i,val in enumerate(steady_state)}
        transient_dict = {i:bool(val==False) for i,val in enumerate(steady_state)}
        manual_params_sto = {
            'transient': transient_dict, 'steady_state': steady_state_dict}
        # make package
        sto_new = refine_package(
            gwf_orig, 
            'sto', 
            gwf_new, 
            convert_df, 
            manual_params=manual_params_sto)
    
    # (NPF) Node property flow
    npf_new = refine_package(gwf_orig, 'npf', gwf_new, convert_df)
    
    # Boundary Conditions
    # basic boundary conditions only (for now)
    handled_boundary_packs = ['drn', 'ghb', 'riv', 'wel', 'chd']
    bound_pack_names = (
        gwf_package_df
        .query(f'ftype in {handled_boundary_packs}')
        .pname
        .to_list()
    )
    # loop over all names and refine
    for pack_name in bound_pack_names:
        _ = refine_package(gwf_orig, pack_name, gwf_new, convert_df)
    
    # (IC) Initial conditions
    ic_new = refine_package(gwf_orig, 'ic', gwf_new, convert_df)
    
    # (OC) Output control
    if gwf_orig.get_package('oc') != None:
        oc_new = copy_package(gwf_orig, 'oc', gwf_new)
    
    # (HDOBS) Head observation
    if gwf_orig.get_package('hdobs') != None:
        hdobs_orig = gwf_orig.get_package('hdobs')
        new_continuous_data = {}
        for fout, data in hdobs_orig.continuous.data.items():
            ids = [convert_df.loc[[id], 'cellid'].values[0] for id in data['id']]
            new_data = data.copy()
            new_data['id'] = ids
            new_continuous_data[fout] = new_data
        hdobs_new = copy_package(
            gwf_orig, 
            'hdobs', 
            gwf_new, 
            manual_params = {'continuous': new_continuous_data}
        )
    return sim_new

def grid_to_quadtree(modelgrid_orig, refine_gdf, refine_levels,
    exe_name='gridgen_x64.exe', layers=None, tempdir='temp'):
    """
    Generate a quadtree grid from a structured grid using refinement features.

    This function creates a quadtree grid based on a structured grid and a set of 
    refinement features. It uses the Gridgen utility to perform the refinement 
    and build the quadtree grid.

    Parameters:
    modelgrid_orig (flopy.discretization.StructuredGrid): The original structured grid.
    refine_gdf (geopandas.GeoDataFrame): A GeoDataFrame containing the refinement 
        features. Each feature should have a geometry and an associated refinement level.
    refine_levels (int or list of int): The refinement level(s) for each feature in 
        `refine_gdf`. If a single integer is provided, it is applied to all features.
    exe_name (str, optional): The name of the Gridgen executable. Defaults to 'gridgen_x64.exe'.
    layers (list of list of int, optional): A list specifying which layers to refine 
        for each feature. Defaults to all layers for all features.
    tempdir (str, optional): The temporary directory where Gridgen will perform its 
        operations. Defaults to 'temp'.

    Returns:
    flopy.utils.gridgen.Gridgen: A Gridgen object representing the refined quadtree grid.
    """
    # add type column 
    refine_gdf = refine_gdf.assign(gtype = lambda x: 
            [typ.lower().replace('string', '') for typ in x.geometry.type])
    # make a gridgen object from original model
    gridgen = Gridgen(modelgrid_orig, model_ws=tempdir, exe_name='gridgen_x64.exe')
    # # add refinement features
    if layers is None:
        layers = [list(range(modelgrid_orig.nlay))]*(refine_gdf.shape[0])
    if isinstance(refine_levels, int):
        refine_levels = [refine_levels] * (refine_gdf.shape[0])
    for i, irefine_level in enumerate(refine_levels):
        igdf = refine_gdf.iloc[i, :]
        igtype = igdf.gtype
        igeom = [igdf.geometry]
        ilayers = layers[i]
        gridgen.add_refinement_features(
            igeom,
            igtype,
            irefine_level,
            ilayers
        )
    gridgen.build()
    return gridgen

def make_grid_relate_table(gridgen):
    """
    Create a grid relate table from a Gridgen quadtree grid.

    This function generates a DataFrame that relates the original structured grid 
    cells to the refined quadtree grid cells. It includes information such as 
    cell IDs for the DIS, DISU, and DISV grids, as well as proportional areas.

    Parameters:
    gridgen (flopy.utils.gridgen.Gridgen): A Gridgen object representing the quadtree grid.

    Returns:
    pd.DataFrame: A DataFrame containing the grid relate table with columns:
        - 'cellid_dis' (index): Cell IDs for the DIS grid.
        - 'cellid_disv': Cell IDs for the DISV grid.
        - 'cellid_disu': Cell IDs for the DISU grid.
        - 'prop_dis_area': Proportional area of each DISV cell relative to the original DIS cell.
    """
    # get the quadtree grid info as a geodataframe
    qtgrid = gpd.GeoDataFrame(gridgen.qtra)
    # gridgen.export()
    layer_info = (
        qtgrid
        .groupby('layer')
        .nodenumber
        .count()
        .to_frame()
        .assign(nodes_above = lambda x: x.nodenumber.cumsum() - x.nodenumber[0])
        .nodes_above
    )
    # make a grid relate table
    grid_relate = (
        qtgrid
        .assign(nodenumber = lambda x: x.nodenumber - 1)
        .assign(cellid_dis = lambda x: [(
            int(l), int(r), int(c)) for l,r,c in zip(x.layer, x.row, x.col)])
        .assign(cellid_disu = lambda x: x.nodenumber)
        .join(layer_info, on='layer')
        .assign(cell2d = lambda x: x.nodenumber - x.nodes_above)
        .assign(cellid_disv = lambda x: [(
            int(l), int(c2)) for l,c2 in zip(x.layer, x.cell2d)])
        .assign(cell_area = lambda x: x.area)
        .assign(prop_dis_area = lambda x: 
                np.round(x.cell_area / x.groupby('cellid_dis').cell_area.transform('sum'), 5))
        .set_index('cellid_dis')
        .loc[:, ['cellid_disv', 'cellid_disu', 'prop_dis_area']]
    )
    return grid_relate

def gridgen_intersect_gdf(gridgen, gdf):
    """
    Intersect a GeoDataFrame with a Gridgen quadtree grid.

    This function computes the intersection of geometries in a GeoDataFrame with 
    the quadtree grid layers generated by Gridgen. It returns a dictionary where 
    each key corresponds to a feature ID from the GeoDataFrame, and the value is 
    a DataFrame containing the intersection properties for that feature.

    Parameters:
    gridgen (flopy.utils.gridgen.Gridgen): A Gridgen object representing the quadtree grid.
    gdf (geopandas.GeoDataFrame): A GeoDataFrame containing the features to intersect 
        with the quadtree grid. Each feature should have a geometry and an associated ID.

    Returns:
    dict: A dictionary where keys are feature IDs from the GeoDataFrame, and values 
        are DataFrames containing intersection properties for each layer of the grid.
    """
    intersect_dict = {}
    for id, geom in zip(gdf.id, gdf.geometry):
        # get type
        geom_type = geom.geom_type.lower().replace('string', '')
        # get xy coords
        if geom_type == 'polygon':
            xys = list([list([list(geom.exterior.coords)])])
        else:
            xys = list([list(geom.coords)])
        intersect_prop_list = []
        for layer in range(gridgen.get_nlay()):
            # intersect
            intersect_props = gridgen.intersect(xys, geom_type, layer)
            intersect_prop_list.append(pd.DataFrame(intersect_props))
        # add to dict
        intersect_dict[id] = pd.concat(intersect_prop_list)
    return intersect_dict

def quadtree_refine_dis_gwf(sim_orig, refine_gdf, refine_levels, layers=None, 
    model_name=None, sim_ws_new=None):
    """
    Refine a GWF model using a quadtree grid.

    This function refines a groundwater flow (GWF) model by converting its 
    structured DIS grid to a quadtree-based DISV grid. It uses refinement 
    features provided in a GeoDataFrame and handles the creation of the 
    refined grid, grid relate table, and intersection properties.

    Parameters:
    sim_orig (flopy.mf6.MFSimulation): The original simulation containing the GWF model.
    refine_gdf (geopandas.GeoDataFrame): A GeoDataFrame containing refinement features. 
        Each feature should have a geometry and an associated ID.
    refine_levels (int or list of int): The refinement level(s) for each feature in 
        `refine_gdf`. If a single integer is provided, it is applied to all features.
    layers (list of list of int, optional): A list specifying which layers to refine 
        for each feature. Defaults to all layers for all features.
    model_name (str, optional): The name of the GWF model to refine. If not provided, 
        the first model in the simulation is used.
    sim_ws_new (str, optional): The working directory for the new simulation. If not 
        provided, a default directory is created.

    Returns:
    tuple: A tuple containing:
        - flopy.mf6.MFSimulation: The refined simulation with the GWF model converted to DISV.
        - pd.DataFrame: The grid relate table mapping the original grid to the refined grid.
        - dict: A dictionary of intersection properties for each refinement feature.
    """
    # get gwf model object
    if model_name is None:
        gwf_orig = sim_orig.get_model()
        model_name = gwf_orig.name
    else:
        gwf_orig = sim_orig.get_model(model_name)
    # get modelgrid object
    modelgrid_orig = gwf_orig.modelgrid
    # assign ids to refinement features if not already
    if 'id' not in refine_gdf.columns:
        refine_gdf['id'] = np.nan
    n_ids = len(refine_gdf.id.dropna().unique())
    if n_ids != refine_gdf.shape[0]:
        print('Unnamed refinement features: all features being given generic "id" names')
        refine_gdf['id'] = [f'feature{i}' for i in range(refine_gdf.shape[0])]
    gridgen = grid_to_quadtree(modelgrid_orig, refine_gdf, refine_levels, 
        exe_name='gridgen_x64.exe', layers=layers)
    # get grid relate table
    grid_relate = make_grid_relate_table(gridgen)
    # get disv props
    disv_props = gridgen.get_gridprops_disv()
    # get intersection properties (all layers)
    intersect_prop_dict = gridgen_intersect_gdf(gridgen, refine_gdf)
    # join with grid relate cellid info
    node_ids = (
        grid_relate
        .assign(node = lambda x: x.cellid_disu)
        .set_index('node')
        .drop(['prop_dis_area'], axis=1)
    )
    # loop over all features and make new dictionary
    feature_locs = {}
    for id, int_props in intersect_prop_dict.items():
        int_props_x = int_props.join(node_ids, on='nodenumber')
        feature_locs[id] = int_props_x
    # remove temp files
    shutil.rmtree(gridgen.model_ws)
    # assign sim_ws_new if not already
    if sim_ws_new is None:
        orig_path = sim_orig.sim_path
        orig_name = orig_path.name
        new_name = f'{orig_name}_refined_to_disv'
        sim_ws_new = sim_orig.sim_path.with_name(new_name)
    sim_new = refine_gwf_dis_to_disv(
        sim_orig, model_name, grid_relate, disv_props, sim_ws_new)
    return sim_new, grid_relate, feature_locs

def refine_and_add_wel(sim_ws, well_xy, well_layer, well_refine_level, pump_rate,
    refine_shapes=None, refine_shape_levels=[], sim_ws_new=None, model_name=None,
    silent=True):
    """
    Refine a GWF model using a quadtree grid and add a well package.

    This function refines a groundwater flow (GWF) model by converting its 
    structured DIS grid to a quadtree-based DISV grid. It also adds a well 
    package to the refined model at the specified location and with the 
    specified pumping rate.

    Parameters:
    sim_ws (str): The working directory of the original simulation.
    well_xy (tuple): The (x, y) coordinates of the well.
    well_layer (int): The layer in which the well is located.
    well_refine_level (int): The refinement level for the well location.
    pump_rate (float): The pumping rate for the well.
    refine_shapes (list or shapely.geometry, optional): Additional refinement 
        features as a list of geometries or a single geometry. Defaults to None.
    refine_shape_levels (list of int, optional): Refinement levels for each 
        additional refinement feature. Defaults to an empty list.
    sim_ws_new (str, optional): The working directory for the new simulation. 
        If not provided, a default directory is created.
    model_name (str, optional): The name of the GWF model to refine. If not 
        provided, the first model in the simulation is used.
    silent (bool, optional): Whether to suppress output during simulation 
        writing and running. Defaults to True.

    Returns:
    tuple: A tuple containing:
        - flopy.mf6.MFSimulation: The refined simulation with the well package added.
        - pd.DataFrame: The grid relate table mapping the original grid to the refined grid.
        - tuple: The cell ID of the well in the refined grid.
    """
    # get original simulation
    sim_orig = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
    # make a refinement feature
    # well_xy = well_xyz[:2]
    well_pt = shapely.Point(well_xy)
    if refine_shapes is not None:
        if isinstance(refine_shapes, list):
            ids = ['pumping_well'] + [f'refine_feat{i+1}' for i in range(len(refine_shapes))]
            geoms = [well_pt] + refine_shapes
        else:
            ids = ['pumping_well', 'refine_feat']
            geoms = [well_pt, refine_shapes]
        refine_gdf = gpd.GeoDataFrame({'id': ids, 'geometry': geoms})
    elif isinstance(refine_shapes, list):
        refine_gdf = gpd.GeoDataFrame({
            'id': ['pumping_well'], 'geometry': [shapely.Point(well_xy)]})
    if (refine_shapes is not None) and (len(refine_shape_levels)==0):
        refine_shape_levels = [well_refine_level] * refine_shapes.shape[0]
    refine_levels = [well_refine_level] + refine_shape_levels
    # make a quadtree refined version of the model
    sim_new, grid_relate, _ = quadtree_refine_dis_gwf(
        sim_orig, 
        refine_gdf, 
        refine_levels,
        layers=None, 
        model_name=model_name, 
        sim_ws_new=sim_ws_new
    )
    # get gwf model object
    if model_name is None:
        gwf_new = sim_new.get_model()
        model_name = gwf_new.name
    else:
        gwf_new = sim_new.get_model(model_name)
    # add a wel object
    # get cellid
    well_cellid = (
        well_layer,
        gwf_new.modelgrid.intersect(well_xy[0], well_xy[1])
    )
    wel_spd = {0: [[well_cellid, pump_rate]]}
    # make wel file
    _ = flopy.mf6.ModflowGwfwel(gwf_new, stress_period_data=wel_spd)
    # write and run simulation
    sim_new.write_simulation(silent=silent)
    sim_new.run_simulation(silent=silent)
    return sim_new, grid_relate, well_cellid