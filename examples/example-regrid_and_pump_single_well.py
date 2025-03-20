#%% 
# IMPORT
import os
import numpy as np
import matplotlib.pyplot as plt
import flopy
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mfmodify import refine_and_add_wel
from mfmodify.plotting import make_xs_line_along_dis_grid, plot_xs_across_pt

# INPUT
sim_ws_base = os.path.join('regrid_tests', 'tvgwfm-average_scenario')
sim_ws_new = os.path.join('regrid_tests', 'tvgwfm-average_scenario-wel_refined_6x')
model_name = 'mf6-tv_hist'
# pumping well
well_xy = (2265074, 1392505)
well_layer = 4
pump_rate = -2071
refine_level = 6 # ~25 meters

#%%
# BODY
# loop over a bunch of refinement levels
sim_new, grid_relate, well_cellid  = refine_and_add_wel(
    sim_ws_base, 
    well_xy, 
    well_layer,
    refine_level, 
    pump_rate, 
    sim_ws_new=sim_ws_new,
    model_name=model_name
)
#%%
# Plot results on map and along cross-section along refinement point
# make a line about the refinement point
# get the well x,y
xy = well_xy
# get original gwf
gwf_orig = (
    flopy.mf6.MFSimulation
    .load(sim_ws=sim_ws_base, verbosity_level=0)
    .get_model(model_name=model_name)
)
xs_line = make_xs_line_along_dis_grid(gwf_orig, xy, 2750, along='column')
# plot original results of last sp
# plot
fig_orig = plot_xs_across_pt(gwf_orig, xs_line, xs_ylim=[300, 800])
# plot new results of last sp
# get new model
gwf_new = sim_new.get_model(model_name=model_name)
# plot
fig_orig = plot_xs_across_pt(gwf_new, xs_line, xs_ylim=[300, 800])

#%%



