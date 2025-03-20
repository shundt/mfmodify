#%% 
# IMPORT
import os
import numpy as np
import matplotlib.pyplot as plt
import geopandas as gpd
import shapely
import flopy
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mfmodify.regrid import (
    quadtree_refine_dis_gwf
)
from mfmodify.plotting import (
    make_xs_line_along_dis_grid,
    plot_xs_across_pt,
    plot_compare_lstbud_gwfs
)

# INPUT
sim_ws_orig = os.path.join('regrid_tests', 'tvgwfm-average_scenario')
sim_ws_new = os.path.join('regrid_tests', 'tvgwfm-average_scenario-refined_same_stresses')
refine_xys = [(2300000, 1390000), (2320000, 1410000)]
refine_level = 5

#%%
# BODY
# load original simulation
sim_orig = flopy.mf6.MFSimulation.load(sim_ws=sim_ws_orig, verbosity_level=0)
#%%
# create a refinement geodataframe (simple line, can be point, line, or polygon)
refine_linestring = shapely.geometry.LineString(refine_xys)
refine_gdf = gpd.GeoDataFrame({'id': ['feat1'], 'geometry': [refine_linestring]})

# plot grid and refinement line
fig, ax = plt.subplots(1,1)
sim_orig.get_model().modelgrid.plot(ax=ax)
refine_gdf.plot(ax=ax)

#%%
# refine simulation at point
sim_new, grid_relate, feature_locs = quadtree_refine_dis_gwf(
    sim_orig,
    refine_gdf,
    refine_level,
    sim_ws_new=sim_ws_new
)
#%%
# plot new grid and refinement line
# plot grid and refinement line
fig, ax = plt.subplots(1,1)
sim_new.get_model().modelgrid.plot(ax=ax)
refine_gdf.plot(ax=ax)

#%%
# write files
sim_new.write_simulation()
# run simulation
sim_new.run_simulation()

#%%
# get gwf objects
gwf_orig = (
    flopy.mf6.MFSimulation
    .load(sim_ws=sim_ws_orig, verbosity_level=0)
    .get_model()
)
gwf_new = sim_new.get_model()
#%%
# Compare the lst file budgets
# plot_compare_lstbud_gwfs(gwf_orig, gwf_new, names=['original', 'refined'])
#%%
# Plot results on map and along cross-section along refinement point
# make a line about the refinement point
# get the well x,y
xy = refine_linestring.centroid.coords[0]
xs_line = make_xs_line_along_dis_grid(gwf_orig, xy, 7500, along='row')
# get plots
# plot original water table and xs
fig_orig = plot_xs_across_pt(gwf_orig, xs_line, xs_ylim=[500, 800])
# plot new water table and xs
fig_new = plot_xs_across_pt(gwf_new, xs_line, xs_ylim=[500, 800])


#%%



