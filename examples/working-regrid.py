#%% 
# IMPORT
import os
import numpy as np
import pandas as pd
import shapely
import shutil
import geopandas as gpd
import flopy
from flopy.utils.gridgen import Gridgen
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mfmodify.regrid import (
    refine_gwf_dis_to_disv
)
from mfmodify.utils import (
    copy_sim,
)

from mfmodify.plotting import (
    plot_compare_lstbud_gwfs,
    plot_compare_obs_sim,
    make_xs_line_along_dis_grid,
    plot_xs_across_pt
)

# INPUT
sim_ws_orig = os.path.join('scenario_tests', 'scenario-original_historic_subset')
new_sim_base_dir = os.path.join('regrid_tests')
sim_ws_new = os.path.join(new_sim_base_dir, 'scenario-baseline_regrid')
sim_ws_base = os.path.join(new_sim_base_dir, 'scenario-baseline')
model_name = 'mf6-tv_hist'
# switches
rewrite_base = False
#%%
# FUNCTIONS

def grid_to_quadtree(modelgrid_orig, refine_gdf, refine_level,
    exe_name='gridgen_x64.exe', layers=None, tempdir='temp'):
    # make sure all geometries are of the same type
    geom_types = refine_gdf.geometry.type.unique()
    if len(geom_types) == 1:
        geom_type = geom_types[0].lower().replace('string', '')
    else:
        raise ValueError('Refinement features must have 1 and only 1 geometry type')
    # make a gridgen object from original model
    gridgen = Gridgen(modelgrid_orig, model_ws=tempdir, exe_name='gridgen_x64.exe')
    # # add refinement features
    if layers is None:
        layers = list(range(modelgrid_orig.nlay))
    gridgen.add_refinement_features(
        refine_gdf.geometry.to_list(),
        geom_type,
        refine_level,
        layers
    )
    gridgen.build()
    return gridgen

def make_grid_relate_table(gridgen):
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

def gridgen_intersect_gdf(gridgen, gdf, layer=0):
    intersect_dict = {}
    for id, geom in zip(gdf.id, gdf.geometry):
        # get type
        geom_type = geom.geom_type.lower().replace('string', '')
        # get xy coords
        xys = list([list(geom.coords)])
        # intersect
        intersect_props = gridgen.intersect(xys, geom_type, layer)
        # add to dict
        intersect_dict[id] = pd.DataFrame(intersect_props)
    return intersect_dict
#%%
# BODY
# load objects from existing model
sim_orig = flopy.mf6.MFSimulation.load(sim_ws=sim_ws_orig, verbosity_level=0)
gwf_orig = sim_orig.get_model(model_name=model_name)
# Make a copy of sim to add a single head obs
sim_orig.set_sim_path(sim_ws_base)
if rewrite_base:
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

# refine_feat_shp = os.path.join(
    # 'regrid_tests', 'refinement_features', 'tvgwfm', 'single_linestring.shp')
# # modelgrid_orig = gwf_orig.modelgrid
# refine_gdf = gpd.read_file(refine_feat_shp)
# refine_level=2
# layers=None
# refine_gdf


#%%

def quadtree_refine_dis(sim_orig, refine_gdf, refine_level, layers=None, model_name=None, sim_ws_new=None):
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
    gridgen = grid_to_quadtree(modelgrid_orig, refine_gdf, refine_level, 
        exe_name='gridgen_x64.exe', layers=layers)
    # get the new modelgrid
    modelgrid_new = gridgen.modelgrid
    # get grid relate table
    grid_relate = make_grid_relate_table(gridgen)
    # get disv props
    disv_props = gridgen.get_gridprops_disv()
    # get intersection properties (layer 1)
    intersect_prop_dict = gridgen_intersect_gdf(gridgen, refine_gdf, layer=0)
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
    sim_new = refine_gwf_dis_to_disv(sim_orig, model_name, grid_relate, disv_props, sim_ws_new)
    return sim_new, feature_locs

model_name=None
refine_feat_shp = os.path.join(
    'regrid_tests', 'refinement_features', 'tvgwfm', 'single_point.shp')
refine_gdf = gpd.read_file(refine_feat_shp)
refine_level = 2
pump_q = -1000

def refine_sim_with_pump_well(sim_ws, well_xyz, refine_level, pump_rate, sim_ws_new=None):
    # get original simulation
    sim_orig = flopy.mf6.MFSimulation.load(sim_ws=sim_ws_orig, verbosity_level=0)
    # make a refinement feature
    refine_gdf = gpd.GeoDataFrame({
        'id': ['pumping_well'], 'geometry': [shapely.Point(well_xyz)]})
    # make a quadtree refined version of the model
    sim_new, feature_locs = quadtree_refine_dis(sim_orig, refine_gdf, refine_level, layers=None, model_name=None, sim_ws_new=None)
    # get gwf model object
    if model_name is None:
        gwf_new = sim_new.get_model()
        model_name = gwf_new.name
    else:
        gwf_new = sim_new.get_model(model_name)
    # add a wel object
    wel_cellid = list(feature_locs.values())[0].at[0, 'cellid_disv']
    wel_spd = {0: [[wel_cellid, pump_rate]]}
    # make wel file
    wel_new_pump = flopy.mf6.ModflowGwfwel(gwf_new, stress_period_data=wel_spd)
    # write and run simulation
    sim_new.write_simulation()
    sim_new.run_simulation()
    return sim_new, wel_cellid








  




    # write model files
    # sim_new.write_simulation()
    # run model
    # sim_new.run_simulation()
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
