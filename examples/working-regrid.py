#%% 
# IMPORT
import os
import numpy as np
import pandas as pd
import shapely
import geopandas as gpd
import flopy
from flopy.utils.gridgen import Gridgen
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mfmodify.regrid import (
    refine_package,
    convert_array_dis_to_disv,
    modify_gridgen_disv_props,
)
from mfmodify.utils import (
    get_parameter_set,
    param_dict_from_list,
    copy_package,
    lst_df_from_gwf,
    get_gwf_package_df
)

# INPUT
sim_ws = os.path.join('historic_models', 'model_v1-metric')
new_sim_ws = os.path.join('regrid_tests', 'model-regridded')
model_name = 'mf6-tv_hist'

#%%
# FUNCTIONS

#%%
# BODY
# load objects from existing model
sim_orig = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
gwf_orig = sim_orig.get_model(model_name=model_name)
# get some overall info
tdis_orig = sim_orig.get_package('tdis')
nper = tdis_orig.nper.data

#%%
rerun_grid_make = True
if rerun_grid_make:
    regrid_dir = os.path.join('regrid_tests', 'gridgen_export')
    # THIS COULD BE DONE BETTER
    # create gridgen_export directory if it doesn't already exist
    if not os.path.exists(regrid_dir):
        os.mkdir(regrid_dir)
    # get grid
    grid0 = gwf_orig.modelgrid
    gridgen = Gridgen(grid0, model_ws=regrid_dir, exe_name='gridgen_x64.exe')
    # make some dummy points
    middle_cell = int(grid0.size/6/2) + 30
    middle_cell_vertices_xy = grid0.get_cell_vertices(middle_cell)
    middle_center_xy = np.array((
        np.mean([x[0] for x in middle_cell_vertices_xy])
        , np.mean([x[1] for x in middle_cell_vertices_xy])
    ))
    add_point1 = shapely.geometry.Point(middle_center_xy + 4000)
    # add_point2 = shapely.geometry.Point(middle_center_xy + 10000)
    # write to shapefile
    dummy_pt_gdf = gpd.GeoDataFrame({
        # 'ptid': [0, 1], 'geometry': [add_point1, add_point2]}, 
        'ptid': [0], 'geometry': [add_point1]}, 
        geometry='geometry'
    )
    dummy_pt_gdf.to_file(os.path.join(regrid_dir, 'temp.shp'))
    # make a very simple quadtree grid
    gridgen.add_refinement_features('temp', 'point', 2, [0, 1, 2, 3, 4, 5])
    gridgen.build()
    # export to files
    gridgen.export()


#%%
# THERE IS PROBABLY A BETTER WAY TO DO THIS
# make a table to relate new ids to old ids
# right now I am doing this with the info from the qtgrid.shp file, but we may 
# find other ways of doing this with the modelgrid or gridgen objects instead

# load the qtgrid that I just exported from gridgen
qtgrid = gpd.read_file(os.path.join(regrid_dir, 'qtgrid.shp'))
# get layer information
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

#%%
# Refine gwf model from structured to vertex grid ====
# def refine_gwf_dis_to_disv(sim_orig, model_name, grid_relate, disv_props)

# Get original model
gwf_orig = sim_orig.get_model(model_name)

# Check if all packages in the gwf model can be converted by this function
valid_gwf_package_set = set([
    'dis', 'sto', 'npf', 'ic', 'oc', 'chd', 'wel', 'drn', 'riv', 'ghb'])
gwf_package_df = (
    pd
    .DataFrame(gwf_orig.name_file.packages.array)
    .assign(ftype = lambda x: [ft[:-1] for ft in x.ftype])
)
gwf_package_set = set(gwf_package_df.ftype)
unhandled_paks = gwf_package_set.difference(valid_gwf_package_set)
if len(unhandled_paks)>0:
    raise ValueError(f'The gwf package types {unhandled_paks} are not yet handled')

# Make conversion dataframe
convert_df = (
    grid_relate
    .rename(columns = {'cellid_disv': 'cellid', 'prop_dis_area': 'area_factor'})
    .loc[: ,['cellid', 'area_factor']]
)

# Make a new model
# (MFSIM) Simulation --
params = list(get_parameter_set(sim_orig))
param_dict = param_dict_from_list(sim_orig, params)
sim_new = flopy.mf6.MFSimulation(sim_ws = new_sim_ws,**param_dict)

# (TDIS) Temporal discretization --
tdis_new = copy_package(sim_orig, 'tdis', sim_new)

# (GWF) Groundwater flow model
gwf_new = flopy.mf6.ModflowGwf(sim_new, modelname=model_name)

# (IMS) Iterative model solution
ims_new = copy_package(sim_orig, 'ims', sim_new)
sim_new.register_ims_package(ims_new, [model_name])
#%%

# (DISV) Discretization by vertices
# get disv properties from gridgen
disv_props_orig = gridgen.get_gridprops_disv()
disv_props = modify_gridgen_disv_props(disv_props_orig)
# get original dis object
dis_orig = gwf_orig.get_package('dis')
# convert idomain, to new grid
idomain_disv = convert_array_dis_to_disv(dis_orig.idomain.data, convert_df)
disv_props['idomain'] = idomain_disv
# add other properties from dis
for prop in ['length_units', 'xorigin', 'yorigin', 'angrot']:
    disv_props[prop] = getattr(dis_orig, prop).data
# remove dis grid and add disv grid
disv = flopy.mf6.modflow.ModflowGwfdisv(gwf_new, **disv_props) 
# %%

# Node properties (no multiplier, just assign original values)
# (STO) Storage
# make manual transient and steady state parameter dicts
transient_dict = {0:False}
steady_state_dict = {0:True}
for iper in range(1, nper):
    transient_dict[iper] = True
    steady_state_dict[iper] = False
manual_params_sto = {'transient': transient_dict, 'steady_state': steady_state_dict}
# make package
sto_new = refine_package(gwf_orig, 'sto', gwf_new, convert_df, manual_params=manual_params_sto)

# (NPF) Node property flow
npf_new = refine_package(gwf_orig, 'npf', gwf_new, convert_df)

#%%
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

#%%

# (IC) Initial conditions
ic_new = refine_package(gwf_orig, 'ic', gwf_new, convert_df)

# (OC) Output control
oc_new = copy_package(gwf_orig, 'oc', gwf_new)

# %%
# write model files
sim_new.write_simulation()
sim_new.run_simulation()

# %%
# Compare before and after
import matplotlib.pyplot as plt

# budget obs files
# get dictionary of output filenames
obs_fileout_dict = {
    x.package_name: list(x.continuous.data.keys()) for x in gwf_orig.obs
}
# loop over and load
pack_df_dict = {}
for ipack, ifiles in obs_fileout_dict.items():
    ipack_df_list = []
    for ifile in ifiles:
        for vers, idir in zip(['original', 'regrid'], [sim_ws, new_sim_ws]):
            idf = (
                pd
                .read_csv(os.path.join(idir, ifile))
                .melt(id_vars = ['time'], var_name='bound', value_name='boundobs')
                .assign(package = ipack)
                .assign(file = ifile)
                .assign(version = vers)
            )
            ipack_df_list.append(idf)
    pack_df_dict[ipack] = pd.concat(ipack_df_list)
# plot
for ipack, ipack_df in pack_df_dict.items():
    nbounds = ipack_df.bound.nunique()
    nrows = min(nbounds, 3)
    ncols = int(np.ceil(nbounds / nrows))
    fig, axes_grid = plt.subplots(nrows, ncols, figsize=(10, 7.5))
    if isinstance(axes_grid, np.ndarray):
        axes = axes_grid.flatten()
    else:
        axes = list([axes_grid])
    i=0
    for ibound, idf in ipack_df.groupby('bound'):
        ax = axes[i]
        idf_piv = idf.pivot(index='time', columns='version', values='boundobs')
        idf_piv.plot(ax=ax, title=ibound, linewidth=1)
        i += 1
    fig.tight_layout()

#%%

from mfmodify.scenario import get_sp_data

def lst_df_from_gwf(gwf):
    # get listing file
    lst = gwf.output.list()
    # get timing information
    sim = gwf.simulation
    sp_df = get_sp_data(sim)
    start_datetime = sim.get_package('tdis').start_date_time.data
    # get the listing dataframes
    lst_df_tuples = lst.get_dataframes(start_datetime = start_datetime)
    # get and format the dataframe
    i = 0
    lst_df = (
        lst_df_tuples[i]
        .reset_index(drop=True)
        .join(sp_df.loc[:, ['sp', 'year', 'month']])
        .assign(total_recharge = lambda x: x.TOTAL_IN - x['STO-SS_IN'])
        .assign(total_discharge = lambda x: x.TOTAL_OUT - x['STO-SS_OUT'])
        .assign(net_storage = lambda x: x['STO-SS_OUT'] - x['STO-SS_IN'])
    )
    return lst_df
#%%
lst_df = lst_df_from_gwf(gwf_orig)
#%%
gwf = gwf_orig
# def lst_df_from_gwf_long(gwf):
if 1:
    # get gwf package info
    gwf_package_df = get_gwf_package_df(gwf)
    lst_df_long = (
        lst_df
        .melt(id_vars = ['sp', 'year', 'month'], var_name='pak', value_name='rate')
        .assign(paknum = lambda x: [x.split('_')[0].lower() for x in x.pak])
    )

lst_df_long.head()

#%%
gwf_package_df.head()

#%%
# list budget
lst_df_orig = lst_df_from_gwf(gwf_orig)
lst_df_new = lst_df_from_gwf(gwf_new)
#%%
cols = lst_df_orig.columns.to_list()
cols.remove('month')
cols.remove('year')
cols.remove('sp')

# plot_cols = ['total_recharge', 'total_discharge', 'net_storage']
plot_cols = cols
for plot_col in plot_cols:
    idf_orig = lst_df_orig.loc[:, ['sp', plot_col]]
    idf_new = lst_df_new.loc[:, ['sp', plot_col]]
    # get fig and axes
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(7.5, 10))
    # plot side by side
    ax0.plot(idf_orig.sp, idf_orig.iloc[:, 1], label='original')
    ax0.plot(idf_new.sp, idf_new.iloc[:, 1], label='new')
    ax0.set_title('values')
    ax0.set_xlabel('stress period')
    ax0.set_ylabel(f'{plot_col} in volume per time')
    ax0.legend()
    # plot residual
    idf_resid = idf_orig.iloc[:, 1] - idf_new.iloc[:, 1]
    ax1.plot(idf_orig.sp, idf_resid)
    ax1.set_title('difference')
    ax1.set_xlabel('stress period')
    ax1.set_ylabel(f'{plot_col} residual in volume per time')
    fig.suptitle(plot_col, fontsize=16)
    fig.tight_layout()
    # plot percent residual
    idf_pct_resid = 100 * (idf_resid / idf_orig.iloc[:, 1]).fillna(0)
    ax2.plot(idf_orig.sp, idf_pct_resid, color='black')
    ax2.set_title('percent difference')
    ax2.set_xlabel('stress period')
    ax2.set_ylabel(f'{plot_col} residual in percent')
    fig.suptitle(plot_col, fontsize=16)
    fig.tight_layout()

#%%
def get_gwf_package_df(gwf):
    gwf_package_df = (
        pd
        .DataFrame(gwf.name_file.packages.array)
        .assign(ftype = lambda x: [ft[:-1] for ft in x.ftype])
        .assign(dummy = 1)
        .assign(paknum = lambda x: x.groupby('ftype').dummy.transform('cumsum') - 1)
        .assign(paknum = lambda x: [nm if nm!=0 else '' for nm in x.paknum])
        .assign(pakname = lambda x: [f'{tp}{nm}' for nm, tp in zip(x.paknum, x.ftype)])
        .drop(['dummy', 'paknum'], axis=1)
    )
    return gwf_package_df
# gwf_package_df

# %%
