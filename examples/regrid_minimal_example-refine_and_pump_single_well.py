# This is a minimal example to show just what is needed to refine a model
# around a well and pump at a constant rate. If you need something other than
# a constant pumping rate, use the quadtree_refine_dis_gwf function to refine the
# grid around a point. That function will return a simulation object 
# (along with the cellid of the refinement point and a dataframe relating the
# old and new grids). You can then add a wel package to the groundwater flow 
# (gwf) object of that simulation, placing a well with whatever pumping 
# time-series you want at the location of the refinement point.

# IMPORT
import os
from mfmodify import refine_and_add_wel

# INPUT
sim_ws_base = os.path.join('regrid_tests', 'tvgwfm-average_scenario')
sim_ws_new = os.path.join('regrid_tests', 'tvgwfm-average_scenario-wel_refined_6x')
model_name = 'mf6-tv_hist'
# pumping well
well_xy = (2265074, 1392505)
well_layer = 4
pump_rate = -2000
refine_level = 6 

# BODY
# call function
sim_new, grid_relate, well_cellid  = refine_and_add_wel(
    sim_ws_base, # existing simulation directory
    well_xy, # x,y coordinates of well
    well_layer, # layer of well (only 1)
    refine_level, # quadtree refinement level (number of time to divide cell in 4)
    pump_rate, # constant wel q
    sim_ws_new=sim_ws_new, # new simulation directory
    model_name=model_name # model name (not necessary if only one model in sim)
)



