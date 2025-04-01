# Minimal example to test the refine_and_add_wel function from mfmodify. This example
# uses the feature to provide multiple refinement features and refinement levels. 
# IMPORT
import os
import flopy
import shapely
from mfmodify import refine_and_add_wel
import matplotlib.pyplot as plt

# INPUT
sim_ws_base = os.path.join('regrid_tests', 'tvgwfm', 'average_scenario')
sim_ws_new = os.path.join('regrid_tests', 'tvgwfm', 'average_scenario_pump_multi_refine_features')
model_name = 'mf6-tv_hist'
# pumping well
well_xy = (2265074, 1392505)
well_layer = 4
pump_rate = -2000
well_refine_level = 6
outer_refine_level = 2
well_buffer = 50
outer_buffer = 1609*4

# BODY
# make a couple of shapes
ipt = shapely.geometry.Point(well_xy)
# create the refine shape and level inputs. First is the small area around the well
# to get a bigger area of the smallest cells. Second is a larger area for backward 
# particle tracking around the well, where a refinement is wanted, but not as much
# as the well itself. 
ipoly = ipt.buffer(well_buffer)
ipoly1 = ipt.buffer(1609*4)
refine_shapes = [ipoly, ipoly1]
refine_shape_levels = [well_refine_level, outer_refine_level]  

# call function to refine, build files, and run model
sim_new, grid_relate, well_cellid  = refine_and_add_wel(
    sim_ws_base, # existing simulation directory
    well_xy, # x,y coordinates of well
    well_layer, # layer of well (only 1)
    well_refine_level, # quadtree refinement level (number of time to divide cell in 4)
    pump_rate, # constant wel q
    refine_shapes=refine_shapes, # optional additional refinement areas (other than well point)
    refine_shape_levels=refine_shape_levels, # optional refinement levels of shapes (if different than well)
    sim_ws_new=sim_ws_new, # new simulation directory
    model_name=model_name # model name (not necessary if only one model in sim)
)

# %%
# plot head around well
fig,ax = plt.subplots(1,1,figsize=(10,10))
gwf_new = sim_new.get_model()
hds_new = gwf_new.output.head().get_data(idx=-1)
xmin,ymin,xmax,ymax = ipoly.buffer(2000).bounds
extent = (xmin, xmax, ymin, ymax)
pmv = flopy.plot.PlotMapView(gwf_new, ax=ax, layer=4, extent=extent)
pmv.plot_ibound()
arr = pmv.plot_array(hds_new, vmin=705, vmax=720)
pmv.plot_grid(linewidth=0.25)

# %%
# Compare listing files
sim_orig = flopy.mf6.MFSimulation.load(sim_ws=sim_ws_base, verbosity_level=0)
gwf_orig = sim_orig.get_model()
from mfmodify.plotting import plot_compare_final_vols
fig = plot_compare_final_vols(gwf_orig, gwf_new, names=['orig', 'regine_and_pump'])
