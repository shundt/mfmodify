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
    lst_df_from_gwf_long,
)

# INPUT
sim_ws_orig = os.path.join('scenario_tests', 'scenario-original_historic_subset')
sim_ws_new = os.path.join('regrid_tests', 'scenario-baseline_regrid')
sim_ws_base = os.path.join('regrid_tests', 'scenario-baseline')
model_name = 'mf6-tv_hist'


#%%
# FUNCTIONS

import math
from flopy.utils.postprocessing import get_water_table
from matplotlib import gridspec

def make_xs_line_along_dis_grid(gwf_dis, xy, line_len, along='row'):
    # get grid rotation angle
    angrot_deg = gwf_dis.get_package('dis').angrot.data
    # get angle in radians and shift quadrant if along columns
    if along == 'row':
        angrot_rad = math.radians(angrot_deg)
    elif along == 'column':
        angrot_rad = math.radians(angrot_deg + 90)
    else:
        raise(ValueError('Must be along either row or column'))
    # get line cooridates
    x0 = xy[0] + line_len * math.cos(angrot_rad)
    x1 = xy[0] + line_len * math.cos(angrot_rad + math.pi)
    y0 = xy[1] + line_len * math.sin(angrot_rad)
    y1 = xy[1] + line_len * math.sin(angrot_rad + math.pi)
    xs_line = [(x0, y0), (x1, y1)]
    return xs_line

def plot_xs_across_pt(gwf, xs_line, xs_ylim=None):
    # get data
    kstpkper = gwf.output.head().get_kstpkper()[-1]
    # get heads
    hds = gwf.output.head().get_data(kstpkper=kstpkper)
    # get water table at final step
    wt = get_water_table(hds)
    # create figure and gridspec
    fig = plt.figure(figsize=(7.5, 10))
    gs = gridspec.GridSpec(7, 1)
    # create top and bottom axes
    ax_map = fig.add_subplot(gs[:4, 0])
    ax_xs = fig.add_subplot(gs[4:, 0])
    
    # plot location of line on grid and 
    # get map plotter
    pmv = flopy.plot.PlotMapView(gwf, ax=ax_map)
    # add grid
    lc = pmv.plot_grid(lw=0.5)
    # add line
    ax_map.plot([xs_line[0][0], xs_line[1][0]], [xs_line[0][1], xs_line[1][1]], lw=0.75)
    # set equal scale
    ax_map.set_aspect('equal', adjustable='box')
    
    # plot cross-section 
    # get cross-section plotter
    xs = flopy.plot.PlotCrossSection(
        model=gwf,
        line={'line': xs_line},
        ax=ax_xs,
        geographic_coords=True
    )
    # grid
    lc = xs.plot_grid(zorder=10)
    # head fill
    pc = xs.plot_array(hds, head=hds, alpha=0.5, masked_values=[1e30])
    # head contours
    # get levels
    pc_values = pc.get_array().data
    min_val = min(pc_values)
    max_val = max(pc_values)
    val_range = max_val - min_val
    cont_int = int(np.round(val_range / 10))
    cont_levels = np.arange(np.round(min_val), np.round(max_val) + cont_int, cont_int)
    # plot
    ctr = xs.contour_array(hds, head=hds, levels=cont_levels, colors="b", masked_values=[1e30])
    # water table surface
    surf = xs.plot_surface(wt, masked_values=[1e30], color="blue", lw=2)
    # head contour labels
    labels = pmv.ax.clabel(ctr, inline=0.25, fontsize=8, inline_spacing=0)
    ax_xs.set_ylim(xs_ylim)
    fig.tight_layout()
    return fig

#%%
# BODY
# load objects from existing model
sim_orig = flopy.mf6.MFSimulation.load(sim_ws=sim_ws_orig, verbosity_level=0)
gwf_orig = sim_orig.get_model(model_name=model_name)

# Make a copy of sim to add a single head obs
sim_orig.set_sim_path(sim_ws_base)
# Make an OC package to print both head and budget
oc_orig = flopy.mf6.ModflowGwfoc(
    gwf_orig,
    budget_filerecord=f'{model_name}.bud',
    head_filerecord=f'{model_name}.hds',
    saverecord=[('HEAD', 'ALL'), ('BUDGET', 'ALL')]
)
# write the simulation
sim_orig.write_simulation()
# run simulation
sim_orig.run_simulation()
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
    refine_pt_gdf = gpd.GeoDataFrame({
        # 'ptid': [0, 1], 'geometry': [add_point1, add_point2]}, 
        'ptid': [0], 'geometry': [add_point1]}, 
        geometry='geometry'
    )
    refine_pt_gdf.to_file(os.path.join(regrid_dir, 'refinement_pt.shp'))
    # make a very simple quadtree grid
    gridgen.add_refinement_features('refinement_pt', 'point', 2, [0, 1, 2, 3, 4, 5])
    gridgen.build()
    # export to files
    gridgen.export()

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
    'dis', 'sto', 'npf', 'ic', 'oc', 'chd', 'wel', 'drn', 'riv', 'ghb', 'obs', 'oc'])
gwf_package_df = (
    pd
    .DataFrame(gwf_orig.name_file.packages.array)
    .assign(ftype = lambda x: [ft[:-1].lower() for ft in x.ftype])
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
sim_new = flopy.mf6.MFSimulation(sim_ws = sim_ws_new,**param_dict)

# (TDIS) Temporal discretization --
tdis_new = copy_package(sim_orig, 'tdis', sim_new)

# (GWF) Groundwater flow model
gwf_new = flopy.mf6.ModflowGwf(sim_new, modelname=model_name)

# (IMS) Iterative model solution
ims_new = copy_package(sim_orig, 'ims', sim_new)
sim_new.register_ims_package(ims_new, [model_name])

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
# for prop in ['length_units', 'xorigin', 'yorigin', 'angrot']:
    # disv_props[prop] = getattr(dis_orig, prop).data
disv_props['length_units'] = getattr(dis_orig, 'length_units').data
# remove dis grid and add disv grid
disv = flopy.mf6.modflow.ModflowGwfdisv(gwf_new, **disv_props) 

# Node properties (no multiplier, just assign original values)
# (STO) Storage
if gwf_orig.get_package('sto') != None:
    # get steady-state info from gwf
    steady_state = gwf_orig.modeltime.steady_state
    # make dicts and add to manual param dict
    steady_state_dict = {i:bool(val) for i,val in enumerate(steady_state)}
    transient_dict = {i:bool(val==False) for i,val in enumerate(steady_state)}
    manual_params_sto = {'transient': transient_dict, 'steady_state': steady_state_dict}
    # make package
    sto_new = refine_package(gwf_orig, 'sto', gwf_new, convert_df, manual_params=manual_params_sto)

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

# %%
# write model files
sim_new.write_simulation()
#%%
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
        for vers, idir in zip(['original', 'regrid'], [sim_ws_orig, sim_ws_new]):
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

# list budget
lst_df_long_orig = lst_df_from_gwf_long(gwf_orig)
lst_df_long_new = lst_df_from_gwf_long(gwf_new)
comps = lst_df_long_orig.bcompname.unique()

plot_comps = comps
for plot_comp in plot_comps:
    idf_orig = (
        lst_df_long_orig
        .query(f'bcompname == "{plot_comp}"')
        .loc[:, ['sp', 'rate']]
    )
    idf_new = (
        lst_df_long_new
        .query(f'bcompname == "{plot_comp}"')
        .loc[:, ['sp', 'rate']]
    )
    # get fig and axes
    fig, (ax0, ax1, ax2) = plt.subplots(3, 1, figsize=(7.5, 10))
    # plot side by side
    ax0.plot(idf_orig.sp, idf_orig.iloc[:, 1], label='original')
    ax0.plot(idf_new.sp, idf_new.iloc[:, 1], label='new')
    ax0.set_title('values')
    ax0.set_xlabel('stress period')
    ax0.set_ylabel(f'{plot_comp} in volume per time')
    ax0.legend()
    # plot residual
    idf_resid = idf_orig.iloc[:, 1] - idf_new.iloc[:, 1]
    ax1.plot(idf_orig.sp, idf_resid)
    ax1.set_title('difference')
    ax1.set_xlabel('stress period')
    ax1.set_ylabel(f'{plot_comp} residual in volume per time')
    fig.suptitle(plot_comp, fontsize=16)
    fig.tight_layout()
    # plot percent residual
    idf_pct_resid = 100 * (idf_resid / idf_orig.iloc[:, 1]).fillna(0)
    ax2.plot(idf_orig.sp, idf_pct_resid, color='black')
    ax2.set_title('percent difference')
    ax2.set_xlabel('stress period')
    ax2.set_ylabel(f'{plot_comp} residual in percent')
    fig.suptitle(plot_comp, fontsize=16)
    fig.tight_layout()


# plot results along cross-section along refinement point
# make a line about the refinement point
xy = refine_pt_gdf.geometry.values[0].coords[0]
xs_line = make_xs_line_along_dis_grid(gwf_orig, xy, 12000, along='column')
# get plots
fig_orig = plot_xs_across_pt(gwf_orig, xs_line, xs_ylim=[300, 800])
fig_new = plot_xs_across_pt(gwf_new, xs_line, xs_ylim=[300, 800])

#%%
#%%
# compare heads everywhere
# %%
