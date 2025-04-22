#%%
# IMPORT
import sys
import os
# add path of mfmodify local files to supercede those installed via pip into 
# conda environment
sys.path.insert(0, '..')
from mfmodify import refine_and_add_wel

# INPUT
orig_dir_name = 'treasure_valley_hundt_bartolino_2023'
sim_ws_base = os.path.join('original_models', orig_dir_name)
sim_ws_new = os.path.join('modified_models', f'{orig_dir_name}-test0')

 # pumping well
well_xy = (2265074, 1392505)
well_layer = 4
pump_rate = -2000
refine_level = 6 

#%%
# BODY
# call function
sim_new, grid_relate, well_cellid  = refine_and_add_wel(
    sim_ws_base, # existing simulation directory
    well_xy, # x,y coordinates of well
    well_layer, # layer of well (only 1)
    refine_level, # quadtree refinement level (number of time to divide cell in 4)
    pump_rate, # constant wel q
    sim_ws_new=sim_ws_new, # new simulation directory
    # model_name=model_name # model name (not necessary if only one model in sim)
)

# %%
