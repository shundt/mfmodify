#%%
# IMPORT
import os
from math import gcd
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import flopy
# import mf_modify
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from mfmodify import get_sp_data, scenario_from_repeat_years
from mfmodify.scenario import copy_param_dict 

#%%
# INPUT
sim_ws = os.path.join('historic_models', 'model_v1-metric')
scenario_dir = os.path.join('scenario_tests')
new_base_ws = os.path.join(scenario_dir, 'historic_base')

# FUNCTIONS
def plot_ann_summary_bars(ann_df):
    fig, axes = plt.subplots(ann_df.shape[1],1,figsize=(7, 9.5))
    for colname, ax in zip(ann_df.columns, axes):
        ann_df.loc[:, colname].plot(ax=ax, kind='bar', width=0.9, color='grey')
        ax.set_xlabel('')
        ax.set_ylabel(colname)
        ax.set_title(colname.split('_')[-1])
    fig.tight_layout()

def plot_densities(ax, data_series, color='black', hist=True, avg_bin=3, title='', xrange=None, label=''):
    nbins = len(data_series) // avg_bin
    # sample from orig
    if hist:
        data_series.hist(ax=ax, color=color, bins=nbins, edgecolor='white', density=True, alpha=0.4)
    data_series.plot.density(ax=ax, color=color, label=label)
    mean_value = data_series.mean()
    ax.axvline(x=mean_value, color=color, linestyle='dashed')
    if xrange is not None:
        ax.set_xlim(xrange[0], xrange[1])
    ax.set_xlabel('')
    ax.set_ylabel('probability density')
    ax.set_title(title)
    return ax

def plot_ann_summary_hist(ann_df):
    fig, axes = plt.subplots(ann_df.shape[1],1,figsize=(7, 9.5))
    nbins = ann_df.shape[0] // 3
    for colname, ax in zip(ann_df.columns, axes):
        xrange = (ann_df[colname].min(), ann_df[colname].max())
        plot_densities(ax, ann_df[colname], title=colname.split('_')[-1], )
    fig.tight_layout()
    return fig

def manual_reweight_series(data_series, relative_weights):
    n_bins = len(relative_weights)
    n_rows = len(data_series)
    # add a check to make sure that n_bins is >= number of data series rows
    if n_bins > n_rows:
        raise ValueError("The length of relative_weights must be less than or equal to the number of rows in data_series.")
    # find the lowest common multiple of the lengths
    lcm = abs(n_bins * n_rows) // gcd(n_bins, n_rows)
    # get the number of values in each bin
    bin_len = lcm / n_bins
    # expand the data_series to the size of the lcm
    n_dup = int(lcm / n_rows)
    data_series_exp = (
        pd.concat([data_series]*n_dup)
        .sort_values()
        .to_frame()
        .assign(bin = (np.arange(lcm) // bin_len).astype('int'))
    )
    # get the reweighted series
    series_list = []
    for i_bin, i_df in data_series_exp.groupby('bin'):
        rel_weight = relative_weights[i_bin]
        series_list.extend([i_df] * rel_weight)
    reweight_series = (
        pd
        .concat(series_list)
        .iloc[:,0]
        .sort_values()
    )
    return reweight_series

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

def annual_summary_from_gwf(gwf):
    lst_df = lst_df_from_gwf(gwf)
    # annual summary
    ann_list_df = (
        lst_df
        .groupby('year')
        .sum()
        .assign(complete_year = lambda x: x.month == x.month.max())
        .query('complete_year')
        .loc[:, ['total_recharge', 'total_discharge', 'net_storage']]
    )
    return ann_list_df

def get_idomain_df(gwf):
    # get dis and modelgrid
    dis = gwf.get_package('dis')
    grid = gwf.modelgrid
    # get idomain and flatten it
    idomain = dis.idomain.data
    idomain_flat = idomain.ravel()
    indices_flat_npint = [index.ravel() for index in np.indices(idomain.shape)]
    indices_flat = [getattr(x, 'tolist', lambda: x)() for x in indices_flat_npint]
    # get cellids and nodeids
    cellids = list(zip(*indices_flat))
    nodeids = grid.get_node(cellids)
    # make a dataframe
    idomain_df = pd.DataFrame({
        'nodeid': nodeids,
        'cellid': cellids, 
        'idomain': idomain_flat,
    })
    if len(indices_flat) > 1:
        idomain_df = idomain_df.assign(layer = indices_flat[0])
    if len(indices_flat) == 2:
        idomain_df = idomain_df.assign(icell2d = indices_flat[1])
    elif len(indices_flat) == 3:
        idomain_df = (
            idomain_df
            .assign(row = indices_flat[1])
            .assign(column = indices_flat[2])
        )
    return idomain_df

def get_spaced_cells(node_cells, ncells):
    total_cells = len(node_cells)
    step = total_cells / ncells
    indices = np.round(np.arange(step/2, total_cells, step)).astype('int')
    spaced_cells = (
        node_cells
        .iloc[indices, :]
        .cellid.tolist()
    )
    return spaced_cells

def add_new_hdobs(gwf, hdobs_continuous, digits=5):
    if gwf.get_package('hdobs') is None:
        model_name = gwf.name
        hdobs = flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            gwf,
            digits=digits,
            continuous=hdobs_continuous,
            filename=f'{model_name}.head.obs',
            pname='hdobs'
        )
    else:
        hdobs_orig = gwf.get_package('hdobs')
        hdobs_params = copy_param_dict(hdobs_orig)
        new_cont_dict = {
            **hdobs_params['continuous'],
            **new_hdobs_continuous
        }
        hdobs_params['continuous'] = new_cont_dict
        hdobs_params['digits'] = digits
        hdobs = flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            gwf,
            **hdobs_params
        )
    return hdobs

#%%
# load historical model
sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
# set mf6 file
exe_name = 'mf6.exe'
sim.exe_name = exe_name
# get model name
model_name = sim.model_names[0]
# get model
gwf = sim.get_model(model_name)

#%%
# Make a copy of sim to add a single head obs
sim.set_sim_path(new_base_ws)
# get the active grid
idomain_df = get_idomain_df(gwf)
cells_df = idomain_df.query('idomain==1').query('layer==1')
# find a few evenly spaced cells for obs points
obs_cells = get_spaced_cells(cells_df, 2)
# make a head observation continuous record
new_hdobs_continuous = {
    'demo_hdobs.csv': 
    [[f'hdobs{i}', 'head', cellid] for i, cellid in enumerate(obs_cells)]
}
hdobs = add_new_hdobs(gwf, new_hdobs_continuous)
# write the simulation
sim.write_simulation()

#%%
# run simulation
sim.run_simulation()

# %%
# get annual summary of historical model from listing file
ann_summary_df = annual_summary_from_gwf(gwf)
# get a subset
subset = (2006, 2015)
ann_sub_df = ann_summary_df.loc[subset[0]:subset[1], :]
# get sorted by storage
ann_df_sorted = ann_sub_df.sort_values('net_storage')
# Plot bar charts
plot_ann_summary_bars(ann_summary_df)
plot_ann_summary_bars(ann_sub_df)
plot_ann_summary_bars(ann_df_sorted)
# # plot annual summary histograms
plot_ann_summary_hist(ann_summary_df)
plot_ann_summary_hist(ann_sub_df)

# %%
# Reweight distribution
# make a plot to demonstrate how it works
data_series = ann_sub_df.net_storage
relative_weights = [4, 2, 2, 2, 2, 5]
reweighted_series = manual_reweight_series(data_series, relative_weights)

fig, axes = plt.subplots(3, 1, figsize = (7, 9.5))
# plot sorted years
ax = axes[0]
data_series.sort_values().plot(ax=ax, kind='bar', width=0.9, color='grey')
ax.set_xlabel('')
ax.set_xticklabels('')
ax.set_ylabel(data_series.name)
ax.set_title(data_series.name.split('_')[-1])
# plot bar chart of weights
ax = axes[1]
ax.bar(
    range(len(relative_weights))
    , relative_weights
    , width=0.98
)
ax.set_title('weights')
ax.set_ylabel('relative weight')
ax.set_xlim([-0.52, len(relative_weights) - 0.48])
# plot original and new distributions
ax = axes[2]
xrange = (data_series.min(), data_series.max())
ax = plot_densities(ax, data_series, color='black', hist=False, label='original')
ax = plot_densities(ax, reweighted_series, color='red', hist=False,
    title='original and new distributions', xrange=xrange, label='reweighted')
ax.legend()
fig.tight_layout()

#%%
# TODO: run this in parallel
sample_size = 10
random_states = [0, 1, 27]
write_and_run = True

i = 0
scen_dirs = []
for random_state in random_states:
    i+=1
    # random draws from original
    sample_orig = data_series.sample(n=sample_size, replace=True, random_state=random_state)
    # random draw from reweighted series
    sample_reweight = reweighted_series.sample(n=sample_size, replace=True, random_state=random_state)

    # plot the distributions
    fig, axes = plt.subplots(3, 1, figsize=(7, 9.5))
    nbins = sample_size // 3
    min_val = min(sample_reweight.min(), sample_orig.min())
    max_val = min(sample_reweight.max(), sample_orig.max())
    # sample from orig
    ax = axes[0]
    ax = plot_densities(ax, sample_orig, color='black', hist=True, avg_bin=3
    , title='sampled from original distribution', xrange=(min_val, max_val))
    # sample from reweighted
    ax = axes[1]
    ax = plot_densities(ax, sample_reweight, color='red', hist=True, avg_bin=3
    , title='sampled from reweighted distribution', xrange=(min_val, max_val))
    # side-by-side barchart of scenario years
    ax = axes[2]
    side_by_side = pd.DataFrame({
        'original': sample_orig.values
        , 'reweighted': sample_reweight.values
    })
    side_by_side.plot(kind='bar', ax=ax)
    ax.set_xlabel('scenario year')
    ax.set_ylabel(f'{sample_orig.name} from reused\nhistoric scenario year')
    fig.suptitle(random_state)

    # create scenarios
    if write_and_run:
        ic_year = ann_summary_df.index.tolist()[-1]
        ic_mon_year = (12, ic_year)
        # original dist
        scen_dir_name = os.path.join('scenario_tests', f'scenario-original_dist-{i}')
        scen_dirs.append(scen_dir_name)
        scen_years = sample_orig.index.to_list()
        sim_scen,_ = scenario_from_repeat_years(
            new_base_ws, scen_years, ic_mon_year, scen_dir_name)
        # remove the oc package to not print head and budget binary out
        oc = sim_scen.get_model().get_package('oc')
        if oc is not None:
            oc.remove()
        sim_scen.write_simulation()
        sim_scen.exe_name = exe_name
        sim_scen.run_simulation()
        # reweighted dist
        scen_dir_name = os.path.join('scenario_tests', f'scenario-reweighted_dist-{i}')
        scen_dirs.append(scen_dir_name)
        scen_years = sample_reweight.index.to_list()
        sim_scen,_ = scenario_from_repeat_years(
            new_base_ws, scen_years, ic_mon_year, scen_dir_name)
        oc = sim_scen.get_model().get_package('oc')
        if oc is not None:
            oc.remove()
        sim_scen.write_simulation()
        sim_scen.exe_name = exe_name
        sim_scen.run_simulation()
    
#%%
# import mf_modify
# from importlib import reload
# reload(mf_modify)
# Rerun the last 10 years as-is from the historic scenario
ic_year = ann_summary_df.index.tolist()[-1]
ic_mon_year = (12, ic_year)
# original dist
scen_dir_name = os.path.join('scenario_tests', f'scenario-original_historic_subset')
scen_years = data_series.index.to_list()
sim_scen,_ = scenario_from_repeat_years(
    new_base_ws, scen_years, ic_mon_year, scen_dir_name)
# remove the oc package to not print head and budget binary out
oc = sim_scen.get_model().get_package('oc')
if oc is not None:
    oc.remove()
sim_scen.write_simulation()
sim_scen.exe_name = exe_name
sim_scen.run_simulation()


#%%
from glob import glob
scen_dirs = glob(os.path.join('scenario_tests', '*scenario*/'))
#%%
# load and plot the output
df_list = []
for scen_dir in scen_dirs:
    ihdobs_file = os.path.join(scen_dir, 'demo_hdobs.csv')
    hdobs_df = (
        pd
        .read_csv(ihdobs_file)
        .melt(id_vars='time', var_name='obsname', value_name='head')
        .assign(scenario_name = os.path.basename(os.path.normpath(scen_dir)))
    )
    df_list.append(hdobs_df)
all_hdobs = pd.concat(df_list)

# plot by obsname
for obsname, df in all_hdobs.groupby('obsname'):
    fig, ax = plt.subplots()
    for iscen, df_ts in df.groupby('scenario_name'):
        if 'historic' in iscen:
            ax.plot(df_ts.time, df_ts['head'], color='black', lw=2, label=iscen)
        else:
            ax.plot(df_ts.time, df_ts['head'], label=iscen, lw=1)
    ax.legend(bbox_to_anchor=(1,1))
    ax.set_title(obsname)


# %%

#%%
# get and plot the list and recharge, discharge, and storage values
for scen_dir in scen_dirs:
    sim = flopy.mf6.MFSimulation.load(sim_ws=scen_dir, verbosity_level=0)
    model_name = sim.model_names[0]
    gwf = sim.get_model(model_name)
    ann_list_df = annual_summary_from_gwf(gwf)
    ann_list_df.plot(kind='bar', subplots=True)

# %%
lst_df = lst_df_from_gwf(gwf)

#%%
lst_df

# %%
lst_df.set_index('sp').iloc[:, -3:].plot()
# %%
for col in lst_df.columns:
    fig,ax = plt.subplots()
    lst_df[col].iloc[90:100].plot()
    fig.suptitle(fig)


# %%
gwf.get_package_list()

# %%
scen_years
# %%

