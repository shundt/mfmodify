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
    refine_gwf_dis_to_disv
)
from mfmodify.utils import (
    copy_package,
    copy_sim,
    copy_empty_sim
)

from mfmodify.plotting import (
    plot_compare_lstbud_gwfs,
    plot_compare_obs_sim
)

# INPUT
sim_ws_orig = os.path.join('scenario_tests', 'scenario-original_historic_subset')
new_sim_base_dir = os.path.join('regrid_tests')
sim_ws_new = os.path.join(new_sim_base_dir, 'scenario-baseline_regrid')
sim_ws_base = os.path.join(new_sim_base_dir, 'scenario-baseline')
model_name = 'mf6-tv_hist'

#%%
# FUNCTIONS

import math
from flopy.utils.postprocessing import get_water_table
from matplotlib import gridspec
import matplotlib.pyplot as plt

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

def plot_xs_across_pt(gwf, xs_line, xs_ylim=None, nlevs=10, zoom_rel_line=1.5):
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
    
    # plot cross-section 
    # get cross-section plotter
    xs = flopy.plot.PlotCrossSection(
        model=gwf,
        line={'line': xs_line},
        ax=ax_xs,
        geographic_coords=True
    )
    # grid
    lc = xs.plot_grid(zorder=-1)
    # head fill
    pc = xs.plot_array(hds, head=hds, alpha=0.5, masked_values=[1e30])
    # head contours
    # get levels
    pc_values = pc.get_array().data
    min_val = min(pc_values)
    max_val = max(pc_values)
    val_range = max_val - min_val
    cont_int = int(np.round(val_range / nlevs))
    cont_levels = np.arange(
        np.round(min_val) - 5*cont_int, 
        np.round(max_val) + 5*cont_int, 
        cont_int
    )
    # plot
    ctr = xs.contour_array(hds, head=hds, levels=cont_levels, colors="b", masked_values=[1e30])
    # water table surface
    surf = xs.plot_surface(wt, masked_values=[1e30], color="blue", lw=2)
    # head contour labels
    labels = xs.ax.clabel(ctr, inline=0.25, fontsize=8, inline_spacing=0)

    # plot location of line on grid and contours
    # get a zoom window for the map
    center_x = (xs_line[0][0] + xs_line[1][0]) / 2
    center_y = (xs_line[0][1] + xs_line[1][1]) / 2
    line_length = np.sqrt((xs_line[0][0] - xs_line[1][0])**2 + (xs_line[0][1] - xs_line[1][1])**2)
    map_xlim = [center_x - zoom_rel_line/2*line_length, center_x + zoom_rel_line/2*line_length]
    map_ylim = [center_y - zoom_rel_line/2*line_length, center_y + zoom_rel_line/2*line_length]
    # get map plotter
    pmv = flopy.plot.PlotMapView(gwf, ax=ax_map)
    # add grid
    lc = pmv.plot_grid(lw=0.5)
    # add water table contours
    ctr = pmv.contour_array(wt, levels=cont_levels, linewidths=0.75, cmap='winter_r')
    # add line
    ax_map.plot(
        [xs_line[0][0], xs_line[1][0]], 
        [xs_line[0][1], xs_line[1][1]], 
        lw=0.75,
        color='red'
    )
    # set zoom
    ax_map.set_xlim(map_xlim)
    ax_map.set_ylim(map_ylim)
    # set equal scale
    ax_map.set_aspect('equal', adjustable='box')
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
    gridgen.add_refinement_features('refinement_pt', 'point', 8, [0, 1, 2, 3, 4, 5])
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
# get disv properties from gridgen
disv_props = gridgen.get_gridprops_disv()
# %%
# create regrid simulation
sim_new = refine_gwf_dis_to_disv(sim_orig, model_name, grid_relate, disv_props)
# write model files
sim_new.write_simulation()
# run model
sim_new.run_simulation()
#%%

# Make models with a new pumping well
(welx, wely)  = refine_pt_gdf.geometry.values[0].coords[0]
pump_q = -5000

# Original grid
sim_pump_dis_ws = os.path.join(new_sim_base_dir, 'scenario-new_pump_orig_grid')
sim_pump_dis = copy_sim(sim_orig, sim_pump_dis_ws)
gwf_pump_dis = sim_pump_dis.get_model(model_name)
mg = sim_pump_dis.get_model(model_name).modelgrid
# add pumping well
# get cell location
cell_nolay = tuple(int(x) for x in mg.intersect(welx, wely))
cellid = (0, ) + cell_nolay
wel_spd = {0: [[cellid, pump_q]]}
# make wel file
wel_new_pump = flopy.mf6.ModflowGwfwel(gwf_pump_dis, stress_period_data=wel_spd)
# write and run simulation
sim_pump_dis.write_simulation()
sim_pump_dis.run_simulation()
#%%
# New grid
sim_pump_disv_ws = os.path.join(new_sim_base_dir, 'scenario-new_pump_regrid')
sim_pump_disv = copy_sim(sim_new, sim_pump_disv_ws)
gwf_pump_disv = sim_pump_disv.get_model(model_name)
mg = sim_pump_disv.get_model(model_name).modelgrid
# add pumping well
# get cell location
cell_nolay = mg.intersect(welx, wely)
if isinstance(cell_nolay, int): 
    cellid = (0, cell_nolay)
else:
    cellid = (0, ) + cell_nolay
wel_spd = {0: [[cellid, pump_q]]}
# make wel file
wel_new_pump = flopy.mf6.ModflowGwfwel(gwf_pump_disv, stress_period_data=wel_spd)
# write and run simulation
sim_pump_disv.write_simulation()
sim_pump_disv.run_simulation()

#%%

#%%

#%%
# Compare all obs file output
plot_compare_obs_sim(sim_pump_dis.sim_path, sim_pump_disv.sim_path)
# Compare budget time series from the list files
plot_compare_lstbud_gwfs(gwf_pump_dis, gwf_pump_disv)

#%%
# Plot results on map and along cross-section along refinement point
# make a line about the refinement point
xy = refine_pt_gdf.geometry.values[0].coords[0]
xs_line = make_xs_line_along_dis_grid(gwf_orig, xy, 1750, along='column')
# get plots
# fig_orig = plot_xs_across_pt(gwf_orig, xs_line, xs_ylim=[300, 800])
# fig_new = plot_xs_across_pt(gwf_new, xs_line, xs_ylim=[300, 800])
fig_orig = plot_xs_across_pt(gwf_pump_dis, xs_line, xs_ylim=[300, 800], nlevs=20, zoom_rel_line=0.5)
fig_new = plot_xs_across_pt(gwf_pump_disv, xs_line, xs_ylim=[300, 800], nlevs=30, zoom_rel_line=0.5)

# compare heads everywhere
# %%
