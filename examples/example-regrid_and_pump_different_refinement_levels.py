#%% 
# IMPORT
import os
from glob import glob
import numpy as np
import matplotlib.pyplot as plt
import flopy
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mfmodify.regrid import (
    refine_and_add_wel
)
from mfmodify.plotting import (
    make_xs_line_along_dis_grid,
    plot_xs_across_pt,
    plot_interpolated_ws
)

# INPUT
new_sim_base_dir = os.path.join('regrid_tests')
sim_ws_base = os.path.join(new_sim_base_dir, 'scenario-baseline')
model_name = 'mf6-tv_hist'
# pumping well
well_xy = (2312990, 1400989)
well_layer = 0
pump_rate = -2500

#%%
# BODY
# loop over a bunch of refinement levels
for refine_level in range(1, 11):
    sim_ws_new = os.path.join(
        new_sim_base_dir, 
        f'scenario-baseline-pumping-refined_{refine_level}')
    sim_new, grid_relate, well_cellid  = refine_and_add_wel(
        sim_ws_base, 
        # well_xyz, 
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
xs_line = make_xs_line_along_dis_grid(gwf_orig, xy, 1750, along='row')
#%%
# get plots
# find all simulation folders
sim_folders = glob(os.path.join(new_sim_base_dir, 'scenario-baseline-pumping-refined_*'))
refine_levels = [int(f.split('_')[-1]) for f in sim_folders]
sim_folders_sorted = [x for _, x in sorted(zip(refine_levels, sim_folders))]
for sim_folder in sim_folders_sorted:
    sim_name = sim_folder.split(os.sep)[-1]
    # load simulation
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_folder, verbosity_level=0)
    # get model
    gwf = sim.get_model(model_name=model_name)
    # plot
    fig_new = plot_xs_across_pt(gwf, xs_line, xs_ylim=[300, 800])

# sort folders by refine level
# get a list of colors for each folder from the matplotlib cmap cmap_plot
cmap_plot = plt.get_cmap('tab10')
cmap_vals = np.linspace(0, 1, len(sim_folders_sorted))
colors = [cmap_plot(c) for c in cmap_vals]
# create figure and ax
fig, ax = plt.subplots(figsize=(10, 7.5))
for sim_folder, color in zip(sim_folders_sorted, colors):
    sim_name = sim_folder.split(os.sep)[-1].split('-')[-1]
    # load simulation
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_folder, verbosity_level=0)
    # get model
    gwf = sim.get_model(model_name=model_name)
    # get minimum node size
    min_areas = gwf.modelgrid.geo_dataframe.geometry.area.min()
    min_length = int(np.round(np.sqrt(min_areas),0))
    # plot name
    label = f'{sim_name} - (min node size: {min_length} m)'
    # plot
    int_surf = plot_interpolated_ws(gwf, xs_line, ax, color=color, label=label)
# add legend
ax.legend()
xs,_ = int_surf[0].get_data()
ax.set_ylim([-200, 1000])
ax.set_xlim([xs[1], xs[-3]])
fig
#%%



