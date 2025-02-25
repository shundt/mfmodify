# IMPORT
import os
import inspect
import numpy as np 
import pandas as pd
import flopy

# FUNCTIONS
# functions within create_scenario_from_years
def round_date(ts, frequency='M'):
    if frequency in ['M', 'Y']:
        if frequency == 'M':
            offsetter = pd.offsets.MonthBegin
        elif frequency == 'Y':
            offsetter = pd.offsets.YearBegin
        nxt = ts + offsetter()
        prv = offsetter().rollback(ts)
        midpoint = prv + (nxt - prv) / 2
        if ts < midpoint:
            rounded = prv
        else:
            rounded = nxt
    elif frequency in ['h', 'min', 's', 'ms', 'D']:
        rounded = ts.round(frequency)  
    else:
        raise ValueError("Frequency must be 'Y', 'M', 'D', 'h', 'min', 's', or 'ms")
    return rounded

def get_sp_data(sim, snap_dates='M'):
    tdis = sim.get_package('tdis')
    start_date_time = pd.to_datetime(tdis.start_date_time.data)
    sp_df = (
        pd 
        .DataFrame(tdis.perioddata.array) 
        .assign(endtime = lambda x: x.perlen.cumsum())
        .assign(starttime = lambda x: [0] + x.endtime[:-1].to_list())
        .assign(start_date_time = lambda x: x.starttime.map(lambda x: start_date_time + pd.Timedelta(days=x)))
        .assign(start_date_time = lambda x: x.start_date_time.dt.round('min')) 
        .assign(start_date = lambda x: pd.to_datetime(x.start_date_time.dt.date))
        .assign(snap_date = lambda x: x.start_date.apply(round_date, args=snap_dates))
        .assign(year = lambda x: x.snap_date.dt.year)
        .assign(month = lambda x: x.snap_date.dt.month)
        .assign(sp = lambda x: x.index)
        .loc[:, ['sp', 'start_date', 'snap_date', 'year', 'month', 'perlen', 'nstp', 'tsmult', 'endtime', 'starttime', 'start_date_time']]
    )
    return sp_df
def scenario_years_to_sps(scenario_years, sp_df):
    scenario_sps = []
    for year in scenario_years:
        i_sps = (
            sp_df
            .query(f'year == {year}')
            .sp
            .to_list()
        )
        scenario_sps.extend(i_sps)
    return scenario_sps

def get_parameter_set(pack, rem_att_set=set([])):
    # get list of attributes
    pack_att_set = set(vars(pack).keys())
    # Get a list of parameters for instantiating the class
    constructor_params = inspect.signature(pack.__init__).parameters
    param_name_set = set(constructor_params.keys())
    # get list of attributes to extract and store in dictionary to unpack on instantiation
    attribute_trans_list = (pack_att_set.intersection(param_name_set) - rem_att_set)
    return attribute_trans_list

def param_dict_from_list(pack, param_list):
    pack_param_dict = {}
    for att in param_list:
        att_val = getattr(pack, att)
        if att_val is None:
            pack_param_dict[att] = att_val
        elif isinstance(att_val, str):
            pack_param_dict[att] = att_val
        elif isinstance(att_val, bool):
            pack_param_dict[att] = att_val
        elif att_val.has_data():
            pack_param_dict[att] = att_val.get_data()
    return pack_param_dict

def retime_package_sp_data(param_values, new_sp_df):
    # make conversion array for stress periods
    # sp_data_dict = pack.stress_period_data.data
    sp_data_dict = param_values
    sp_data_sps = list(sp_data_dict.keys())
    convert_df = (
        new_sp_df
        .assign(sp_from_orig_sp_data = [closest_value_below(x, sp_data_sps) for x in new_sp_df.sp.values])
        .sort_values(['new_sp'])
        .assign(diff_from_previous = lambda x: [True] + list(x.sp_from_orig_sp_data.values[:-1] != x.sp_from_orig_sp_data.values[1:]))
        .query('diff_from_previous==True')
        .loc[:, ['sp_from_orig_sp_data', 'new_sp']]
    )
    # new dictionary to hold converted sp data
    new_sp_data_dict = {}
    # loop over stress periods
    for orig_sp, new_sp in convert_df.itertuples(index=False, name=None):
        new_sp_data_dict[new_sp] = sp_data_dict[orig_sp]
    return new_sp_data_dict

def get_objects_from_pack(pack, attribute_name):
    object_list = []
    for i in range(1000):
        try:
            obj = getattr(pack, attribute_name)[i]  # Use getattr to access the attribute
            object_list.append(obj)
        except ValueError:
            break
    return object_list

def get_ts_objects(pack):
    return get_objects_from_pack(pack, 'ts')

def get_obs_objects(pack):
    return get_objects_from_pack(pack, 'obs')

def convert_timeseries(timeseries, convert_df):
    # make the timeseries a dataframe
    timeseries_df = pd.DataFrame(timeseries)
    # get all times as an array
    ts_times = timeseries_df.ts_time.values
    # find a maximum time
    maxtime = convert_df.endtime.max() + 1
    # new_sp_df
    new_timeseries_df_list = []
    # group the conversions
    time_adjust_groups = (
        convert_df
        .assign(time_adjust = lambda x: x.starttime - x.orig_starttime)
        .groupby('time_adjust')
        .agg({'orig_starttime':'min', 'orig_endtime':'max'})
    )
    n_groups = time_adjust_groups.shape[0]
    i = 0 
    for i_adjust, i_min, i_max in time_adjust_groups.itertuples(index=True):
        i+=1
        # find the first time to use
        first_time = ts_times[ts_times<=i_min][-1]
        # find the last time to use
        last_time = ts_times[ts_times>=i_max][0]
        # set min and max filters
        i_mintime = 0.00001
        i_maxtime = maxtime - 0.00001
        if i==1:
            i_mintime = 0
        elif i==n_groups:
            i_maxtime = maxtime
        # get all times and adjust them
        i_timeseries_df = (
            timeseries_df
            .query(f'ts_time>={first_time}')
            .query(f'ts_time<={last_time}')
            .assign(ts_time = lambda x: x.ts_time + i_adjust)
            .assign(ts_time = lambda x: x.ts_time.clip(0, maxtime))
            .query(f'ts_time >= {i_mintime}')
            .query(f'ts_time <= {i_maxtime}')
        )
        new_timeseries_df_list.append(i_timeseries_df)
    # concatenate all
    new_timeseries_df =  pd.concat(new_timeseries_df_list).sort_values(by='ts_time')
    if new_timeseries_df.ts_time.values[-1] < maxtime:
        last_row = new_timeseries_df.iloc[[-1]].copy()
        last_row['ts_time'] = maxtime
        new_timeseries_df = (
            pd
            .concat([new_timeseries_df, last_row], ignore_index=True)
            .reset_index(drop=True)
        )
    new_timeseries = new_timeseries_df.to_records(index=False)
    return new_timeseries

def make_retimed_package_param_dict(pack, convert_df=None, transient_params=['stress_period_data']):
    # get list of attributes to use as parameters in instantiating object
    # rem_att_set = set(['loading_package', 'stress_period_data'])
    rem_att_set = set(['loading_package'])
    param_set = get_parameter_set(pack, rem_att_set=rem_att_set)
    # # find array data
    # array_class = flopy.mf6.data.mfdataarray.MFArray
    # array_params = set([
    #     att for att in param_set if (isinstance(getattr(pack, att), array_class)
    #     and getattr(pack, att).has_data())
    # ])
    param_list = list(param_set - set([]))
    # create parameter dictionary 
    # get those from attribute list
    # pack_param_dict = {}
    # for att in param_list:
    #     att_val = getattr(pack, att)
    #     if att_val.has_data():
    #         pack_param_dict[att] = att_val.get_data()
    pack_param_dict = param_dict_from_list(pack, param_list)
    # add others manually
    pack_param_dict['pname'] = pack.package_name
    pack_param_dict['filename'] = pack.filename
    # stress-period data
    # transient information to retime
    for transient_param in transient_params:
        if convert_df is not None:
            if transient_param in pack_param_dict.keys():
                # print(f'Converting {transient_param} for package "{pack.package_name}"')
                param_values = pack_param_dict[transient_param]
                # convert stress periods
                sp_data_new = retime_package_sp_data(param_values, convert_df)
                pack_param_dict[transient_param] =  sp_data_new
        # else:
            # print(f'No conversion dataframe (convert_df) provided for package "{pack.package_name}"')
    return pack_param_dict

def retime_package(sim_or_gwf_orig, pack_name, sim_or_gwf_new, convert_df=None, manual_params={}, transient_params=['stress_period_data']):
    pack = sim_or_gwf_orig.get_package(pack_name)
    # return pack
    pack_class = pack.__class__
    # get package paramters
    pack_param_dict = make_retimed_package_param_dict(pack, convert_df=convert_df, transient_params=transient_params)
    # get manual parameters
    for att, val in manual_params.items():
        pack_param_dict[att] = val
    # instantiate package
    pack_new = pack_class(sim_or_gwf_new, **pack_param_dict)
    # convert and reassociate timeseries packages
    if hasattr(pack, 'ts'):
        ts_objects = get_ts_objects(pack)
        if len(ts_objects) > 0:
            i = 0
            for ts_obj in ts_objects:
                # get the package param dict
                ts_dict = make_retimed_package_param_dict(ts_obj)
                # convert the time-series
                ts_dict['timeseries'] = convert_timeseries(ts_dict['timeseries'], convert_df)
                if i==0:
                    ts_dict['pname'] = f'{pack_name}_ts'
                    pack_new.ts.initialize(**ts_dict)
                else:
                    ts_dict['pname'] = f'{pack_name}_ts{i}'
                    pack_new.ts.append_package(**ts_dict)
                i+=1
    # reassociate observation packages
    if hasattr(pack, 'obs'):
        obs_objects = get_obs_objects(pack)
        if len(obs_objects) > 0:
            i = 0
            for obs_obj in obs_objects:
                # get the package param dict
                obs_dict = make_retimed_package_param_dict(obs_obj)
                # rename and add to parent package
                if i==0:
                    obs_dict['pname'] = f'{pack_name}_obs'
                    pack_new.obs.initialize(**obs_dict)
                else:
                    obs_dict['pname'] = f'{pack_name}_obs{i}'
                    pack_new.obs.append_package(**obs_dict)
                i+=1
    return pack_new

def closest_value_below(value, list):
    list_filtered = np.array(list)[np.array(list)<=value]
    closest_value = max(list_filtered, key=lambda x: x - value)
    return closest_value

def scenario_from_repeat_years(sim_ws, scenario_years, ic_mon_year, new_sim_ws='dummy', start_date_time='2050-01-01t00:00:00'):
    """
    Create a simulation scenario by repeating specified years.

    Parameters:
    sim_ws (str): Path to the simulation workspace.
    scenario_years (list of int): List of years to repeat in the scenario.
    ic_mon_year (tuple): Initial condition month and year as a tuple (month, year).
    new_sim_ws (str, optional): Path to the new simulation workspace. Default is 'dummy'.

    Returns:
    tuple: A tuple containing:
        - sim_new (MFSimulation): The new simulation object.
        - new_pack_dict (dict): A dictionary of new package objects.

    This function loads a simulation from the specified workspace, extracts stress period
    information, and creates a new scenario by repeating the specified years. The initial
    condition is set based on the provided month and year.
    """
    # get information for defining stress periods to repeat
    # load simulation
    print(f'Loading simulation from {sim_ws}') 
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
    # get timing information
    sp_df = get_sp_data(sim)
    # get stress periods
    scenario_sps = scenario_years_to_sps(scenario_years, sp_df)
    n_sps = len(scenario_sps)
    # make new sp information
    new_sp_df = (
        sp_df
        .loc[scenario_sps, :]
        .assign(orig_starttime = lambda x: x.starttime)
        .assign(orig_endtime = lambda x: x.endtime)
        .assign(endtime = lambda x: x.perlen.cumsum())
        .assign(starttime = lambda x: [0] + x.endtime[:-1].to_list())
        .assign(new_sp = lambda x: range(x.shape[0]))
    )
    # get initial condition stress period
    ic_sp = (
        sp_df
        .query(f'month == {ic_mon_year[0]} and year == {ic_mon_year[1]}')
        .sp
        .item()
    )
    # load gwf
    gwf = sim.get_model()
    if gwf is None:
        raise ValueError("Model could not be loaded. Please check the simulation workspace.")
    # get model name
    model_name = gwf.name # type: ignore
    
    # Make a new simulation -----------
    # Simulation-level modules
    print('generating new simulation with modified sp and timeseries info')
    # (MFSIM) Simulation
    params = list(get_parameter_set(sim))
    param_dict = param_dict_from_list(sim, params)
    sim_new = flopy.mf6.MFSimulation(sim_ws = new_sim_ws,**param_dict)
    # (TDIS) Temporal discretization
    tdis = sim.get_package('tdis')
    params_tdis = list(get_parameter_set(tdis))
    param_dict_tdis = param_dict_from_list(tdis, params_tdis)
    param_dict_tdis['nper'] = n_sps
    param_dict_tdis['start_date_time'] = start_date_time
    param_dict_tdis['perioddata'] = new_sp_df.loc[:, ['perlen', 'nstp', 'tsmult']].to_records(index=False)
    new_pack_dict = {}
    new_pack_dict['tdis'] = flopy.mf6.ModflowTdis(sim_new, **param_dict_tdis)
    
    # Groundwater flow model modules -----------
    # (GWF) Groundwater Flow Model
    gwf_new = flopy.mf6.ModflowGwf(sim_new, modelname=gwf.name) # type: ignore
    # (IMS) Iterative model solution (simulation level, but I'm doing it here because of model names)
    new_pack_dict['ims'] = retime_package(sim, 'ims', sim_new)
    sim_new.register_ims_package(new_pack_dict['ims'], [gwf_new.name])
    # modules with no time info: 
    for pack_name in ['dis', 'npf']:
        new_pack_dict[pack_name] = retime_package(gwf, pack_name, gwf_new)
    # modules with time info:
    # sto
    sto = gwf.get_package('sto') # type: ignore
    param_dict_sto = make_retimed_package_param_dict(sto, transient_params=[])
    param_dict_sto['transient'] = {sp:True for sp in range(n_sps)}
    param_dict_sto['steady_state'] = {sp:False for sp in range(n_sps)}
    new_pack_dict['sto'] = flopy.mf6.ModflowGwfsto(gwf_new, **param_dict_sto)
    # oc
    if gwf.get_package('oc') is not None:
        new_pack_dict['oc'] = retime_package(gwf, 'oc', gwf_new, convert_df=new_sp_df, transient_params=['printrecord', 'saverecord'])
    # hdobs
    if gwf.get_package('hdobs') is not None:
        new_pack_dict['hdobs'] = copy_package(gwf, 'hdobs', gwf_new)
    # ic
    # get hds values from appropriate sp
    headfilename = os.path.join(f'{sim_ws}', f'{model_name}.hds')
    hds = flopy.utils.HeadFile(headfilename, precision='double')
    ic_heads = hds.get_data(kstpkper=(0,ic_sp))
    # grab package and set new strt values
    ic = copy_package(gwf, 'ic', gwf_new)
    ic.strt = ic_heads
    new_pack_dict['ic'] = ic
    # boundary files 
    # drn # ghb # riv # wel
    bound_pack_names = []
    for pack_type in ['drn', 'ghb', 'riv', 'wel']:
        packs = gwf.get_package(pack_type) # type: ignore
        if isinstance(packs, list):
            bound_pack_names.extend([pack.package_name for pack in packs])
        else:
            bound_pack_names.append(packs.package_name) # type: ignore
    for pack_name in bound_pack_names:
        new_pack_dict[pack_name] = retime_package(gwf, pack_name, gwf_new, convert_df=new_sp_df)
    # return results
    return sim_new, new_pack_dict

# functions within create_scenario_from_average_years
def get_scenario_sp_lut(scenario_years_and_weights, sp_df):
    """
    Generate a lookup table for stress periods based on scenario years and their weights.

    Parameters:
    scenario_years_and_weights (list of list of tuples): A list where each element is a list of tuples,
        with each tuple containing a year (int) and its corresponding weight (float).
    sp_df (DataFrame): DataFrame containing stress period information including 'year', 'month', 'sp',
        and others.

    Returns:
    sp_lut_df: A DataFrame with columns 'sp', 'new_sp', 'year', 'month', and 'weight' mapping the original
        stress periods to the new stress periods based on the provided scenario years and weights.

    This function creates a lookup table that maps each stress period in the scenario to the corresponding
    original stress period based on the provided scenario years and their weights. The weights are normalized
    to sum to 1 for each set of scenario years.
    """
    max_sp = -1
    all_sp_info = []
    for yr_wt in scenario_years_and_weights:
        # convert to df and enforce weights sum to 1
        yr_wt_sp =  (
            pd
            .DataFrame(yr_wt, columns=['year', 'weight'])
            .assign(weight = lambda x: x.weight / x.weight.sum())
        )
        # get the sp info for the year
        yr_sp_info = (
            yr_wt_sp
            .join(sp_df.set_index('year'), on='year')
            .assign(new_sp = lambda x: x.month + max_sp)
            .loc[:, ['sp', 'new_sp', 'year', 'month', 'weight']]
        )
        max_sp = yr_sp_info.new_sp.max()
        all_sp_info.append(yr_sp_info)
    sp_lut_df = (
        pd
        .concat(all_sp_info)
        .sort_values(['new_sp', 'sp'])
        
    )
    return sp_lut_df

def copy_param_dict(pack):
    # get list of attributes to use as parameters in instantiating object
    rem_att_set = set(['loading_package'])
    param_set = get_parameter_set(pack, rem_att_set=rem_att_set)
    param_list = list(param_set - set([]))
    # create parameter dictionary 
    # get those from attribute list
    pack_param_dict = param_dict_from_list(pack, param_list)
    # add others manually
    pack_param_dict['pname'] = pack.package_name
    pack_param_dict['filename'] = pack.filename
    return pack_param_dict

def copy_package(sim_or_gwf_orig, pack_name, sim_or_gwf_new, manual_params={}):
    pack = sim_or_gwf_orig.get_package(pack_name)
    # return pack
    pack_class = pack.__class__
    # get package paramters
    pack_param_dict = copy_param_dict(pack) 
    # instantiate package
    pack_new = pack_class(sim_or_gwf_new, **pack_param_dict)
    return pack_new

def sp_data_to_df(sp_data_dict):
    df_list = []
    for sp, ra in sp_data_dict.items():
        df = pd.DataFrame(ra).assign(sp=sp)
        df_list.append(df)
    return pd.concat(df_list)

def sp_data_df_to_dict(sp_data_df):
    sp_data_dict = {}
    for sp, df in sp_data_df.groupby('sp'):
        ra = df.drop(['sp'], axis=1).to_records(index=False)
        sp_data_dict[sp] = ra
    return sp_data_dict


def weight_mean_package_sp_data(param_values, sp_data_lut):
    """
    Generate weighted mean stress period data for a package.

    Parameters:
    param_values (dict): Dictionary containing the original stress period data.
    sp_data_lut (DataFrame): DataFrame containing the lookup table for stress periods and their weights.

    Returns:
    dict: A dictionary of weighted mean stress period data.

    This function creates weighted mean stress period data for a package by applying weighted mean calculations
    to the stress period data based on the provided lookup table. The resulting dictionary can be used to
    instantiate a new package with the modified stress period data.
    """
    # Convert stress period data to DataFrame
    sp_data_dict = param_values
    sp_data_df = sp_data_to_df(sp_data_dict)

    # Get all stress periods in the stress period data
    sp_data_sps = sp_data_df.sp.unique()

    # Find the correct stress period to reference in the new stress data prior to averaging
    orig_sp_weights = (
        sp_data_lut
        .assign(sp_from_orig_sp_data = [closest_value_below(x, sp_data_sps) for x in sp_data_lut.sp.values])
        .drop(['sp', 'year', 'starttime', 'endtime', 'orig_starttime', 'orig_endtime'], axis=1)
        .groupby(['new_sp', 'sp_from_orig_sp_data'])
        .sum()
        .reset_index()
        .sort_values(['new_sp', 'sp_from_orig_sp_data'])
    )

    # Loop through stress periods to find all needed and their weights
    orig_weight_tuple_list = []
    new_sp = []
    for sp, idf in orig_sp_weights.groupby('new_sp'):
        orig_weight_tuples = list(
            idf
            .drop(['new_sp'], axis=1)
            .itertuples(index=False, name=None)
        )
        orig_weight_tuple_list.append(orig_weight_tuples)
        new_sp.append(sp)

    # Combine all stress periods and weights and remove any if the 
    # nearest reference stress period doesn't change for the new stress 
    # period (we don't need sequential sp data if the values & ts reference columns don't change)
    needed_sp_df = (
        pd
        .DataFrame({'new_sp': new_sp, 'orig_sp_weight': orig_weight_tuple_list})
        .assign(needed = lambda x: [True] + [t0!=t1 for t0, t1 in zip(x.orig_sp_weight[1:], x.orig_sp_weight[:-1])])
        .join(orig_sp_weights.set_index('new_sp'), on='new_sp')
        .query('needed')
        .drop(['orig_sp_weight', 'needed'], axis=1)
    )

    # Find the columns to weight and those not to weight
    index_cols = ['new_sp']
    for col, coltype in sp_data_df.dtypes.items():
        # check for numeric type
        if not pd.api.types.is_numeric_dtype(coltype):
            index_cols.append(col) # type: ignore

    sort_columns = index_cols.copy()
    sort_columns.remove('cellid')
    # Join the weights DataFrame with stress period data and apply the multipliers
    needed_old_sp_join = (
        needed_sp_df
        .join(sp_data_df.set_index('sp'), on='sp_from_orig_sp_data')
        .set_index(index_cols)
    )
    new_sp_data_df = (
        needed_old_sp_join
        .drop(['sp_from_orig_sp_data', 'weight'], axis=1)
        .mul(pd.Series(needed_old_sp_join.weight.values, index=needed_old_sp_join.index), axis=0)
        .reset_index()
        .groupby(index_cols)
        .sum()
        .reset_index()
        .sort_values(by=sort_columns)
        .rename(columns={'new_sp': 'sp'})
    )
    

    # Convert the DataFrame back to a dictionary of recarrays
    new_sp_data_dict = sp_data_df_to_dict(new_sp_data_df)

    return new_sp_data_dict

def weight_mean_package_param_dict(pack, sp_data_lut, transient_params):
    """
    Generate a parameter dictionary for a package with weighted mean stress period data.

    Parameters:
    pack (flopy.mf6.mfpackage.MFPackage): The original package object.
    sp_data_lut (DataFrame): DataFrame containing the lookup table for stress periods and their weights.
    transient_params (list of str): List of transient parameters to be processed.

    Returns:
    dict: A dictionary of package parameters with weighted mean stress period data.

    This function creates a parameter dictionary for a package by applying weighted mean calculations
    to the stress period data based on the provided lookup table. The resulting dictionary can be used
    to instantiate a new package with the modified parameters.
    """
    # copy param dict
    pack_param_dict = copy_param_dict(pack)
    # stress-period data
    # transient information to retime
    for transient_param in transient_params:
        if sp_data_lut is not None:
            if transient_param in pack_param_dict.keys():
                param_values = pack_param_dict[transient_param]
                # convert stress periods
                sp_data_new = weight_mean_package_sp_data(param_values, sp_data_lut)
                pack_param_dict[transient_param] =  sp_data_new
    return pack_param_dict

def mean_weight_timeseries(timeseries, sp_data_lut):
    # make the timeseries a dataframe
    timeseries_df = pd.DataFrame(timeseries)
    # get all times as an array
    ts_times = timeseries_df.ts_time.values
    # find a maximum time
    maxtime = sp_data_lut.endtime.max() + 1
    # group the conversions
    time_adjust_groups = (
        sp_data_lut
        .assign(new_year = lambda x: x.new_sp.floordiv(12))
        .assign(time_diff = lambda x: x.starttime - x.orig_starttime)
        .assign(time_adjust = lambda x: x.groupby(['new_year', 'year']).time_diff.transform('min').round(2))
        .groupby(['time_adjust', 'weight', 'new_year', 'year'])
        .agg({'orig_starttime':'min', 'orig_endtime':'max'})
    )
    # get new year information
    all_new_years = (
        time_adjust_groups
        .reset_index()
        .new_year
        .unique()
    )
    max_new_year = np.max(all_new_years)
    # make a dictionary to store new timeseries dataframes by new year
    new_timeseries_df_dict = {new_year:[] for new_year in all_new_years}
    # loop over all adjustment groups
    for (i_adjust, i_weight, i_new_year, _), i_min, i_max in time_adjust_groups.itertuples(index=True):
        # find the first time to use
        first_time = ts_times[ts_times<=i_min][-1]
        # find the last time to use
        last_time = ts_times[ts_times>=i_max][0]
        # set min and max filters
        i_maxtime = i_max + i_adjust - 5
        i_mintime = i_min + i_adjust - 5
        if i_new_year==max_new_year:
            i_maxtime = maxtime
        if i_new_year==0:
            i_mintime = 0
        # get all times and adjust them
        i_timeseries_df = (
            timeseries_df
            .query(f'ts_time>={first_time}')
            .query(f'ts_time<={last_time}')
            .assign(ts_time = lambda x: (x.ts_time + i_adjust).round(2))
            .assign(ts_time = lambda x: x.ts_time.clip(0, maxtime))
            .query(f'ts_time <= {i_maxtime}')
            .query(f'ts_time >= {i_mintime}')
            .sort_values('ts_time')
            .set_index('ts_time')
            .mul(i_weight)
        )
        new_timeseries_df_dict[i_new_year].append(i_timeseries_df)
    # combine all timeseries (already weighted above) for each new year
    all_timeseries_list = []
    for new_year, df_list in new_timeseries_df_dict.items():
        # get all times
        i_all_times = (
            pd
            .concat(df_list)
            .sort_index()
            .index
            .unique()
            .to_list()
        )
        # create a dummy dataframe
        i_timeseries_new_year = df_list[0].reindex(i_all_times, fill_value=0).mul(0)
        for j_timeseries_df in df_list:
            j_timeseries_reindex_df = j_timeseries_df.reindex(i_all_times, method='ffill').bfill()
            i_timeseries_new_year += j_timeseries_reindex_df
        # check to see if any duplicate values
        i_timeseries_new_year_clean = (
            i_timeseries_new_year
            .assign(diff_from_last = lambda x: [True] + list(np.max(x.values[:-1,:] != x.values[1:, :],axis=1)))
            .assign(diff_from_last = lambda x: x.diff_from_last.to_list()[:-1] + [True])
            .query('diff_from_last')
            .drop(['diff_from_last'], axis=1)
        )
        # append to list
        all_timeseries_list.append(i_timeseries_new_year_clean)
    # concatenate all
    new_timeseries_df =  (
        pd
        .concat(all_timeseries_list)
        .sort_index()
        .reset_index()
        .dropna()
        .drop_duplicates()
    )
    new_timeseries = new_timeseries_df.to_records(index=False)
    return new_timeseries

def weight_mean_package(sim_or_gwf_orig, pack_name, sim_or_gwf_new, sp_data_lut, manual_params={}, transient_params=['stress_period_data']):
    pack = sim_or_gwf_orig.get_package(pack_name)
    # return pack class constructor
    pack_class = pack.__class__
    # get the parameter dictionary with weighted sp info
    pack_param_dict = weight_mean_package_param_dict(pack, sp_data_lut, transient_params)
    # instantiate package
    pack_new = pack_class(sim_or_gwf_new, **pack_param_dict)
    # convert and reassociate timeseries packages
    if hasattr(pack, 'ts'):
        ts_objects = get_ts_objects(pack)
        if len(ts_objects) > 0:
            i = 0
            for ts_obj in ts_objects:
                # get the package param dict
                ts_dict = copy_param_dict(ts_obj)
                # convert the time-series
                ts_dict['timeseries'] = mean_weight_timeseries(ts_dict['timeseries'], sp_data_lut)
                if i==0:
                    ts_dict['pname'] = f'{pack_name}_ts'
                    pack_new.ts.initialize(**ts_dict)
                else:
                    ts_dict['pname'] = f'{pack_name}_ts{i}'
                    pack_new.ts.append_package(**ts_dict)
                i+=1
    # reassociate observation packages
    if hasattr(pack, 'obs'):
        obs_objects = get_obs_objects(pack)
        if len(obs_objects) > 0:
            i = 0
            for obs_obj in obs_objects:
                # get the package param dict
                obs_dict = copy_param_dict(obs_obj)
                # rename and add to parent package
                if i==0:
                    obs_dict['pname'] = f'{pack_name}_obs'
                    pack_new.obs.initialize(**obs_dict)
                else:
                    obs_dict['pname'] = f'{pack_name}_obs{i}'
                    pack_new.obs.append_package(**obs_dict)
                i+=1
    return pack_new

def scenario_from_weighted_mean_of_years(sim_ws, new_sim_ws, ic_mon_year, scenario_years_and_weights, start_date_time='2050-01-01t00:00:00'):
    print(f'Loading simulation from {sim_ws}') 
    sim = flopy.mf6.MFSimulation.load(sim_ws=sim_ws, verbosity_level=0)
    # get timing information
    sp_df = get_sp_data(sim)
    # get sp information 
    # - get a lut of new to original stress periods and corresponding weights
    scenario_sp_lut =  get_scenario_sp_lut(scenario_years_and_weights, sp_df)
    # - get a new sp info dataframe
    new_sp_df = (
        scenario_sp_lut
        .join(sp_df.set_index('sp'), on='sp', rsuffix='xx')
        .groupby('new_sp')
        .agg({'perlen': 'mean', 'nstp': 'max', 'tsmult': 'max'})
        .reset_index()
        .assign(endtime = lambda x: x.perlen.cumsum().round(2))
        .assign(starttime = lambda x: (x.endtime - x.perlen).round(2))
    )
    sp_data_lut = (
        scenario_sp_lut
        .drop(['year', 'month'], axis=1)
        .join(new_sp_df.set_index('new_sp'), on='new_sp')
        .join(sp_df.set_index('sp'), on='sp', rsuffix='_o')
        .rename(columns={'starttime_o': 'orig_starttime', 'endtime_o': 'orig_endtime'})
        .loc[:, ['sp', 'new_sp', 'weight', 'starttime', 'endtime', 'orig_starttime', 'orig_endtime', 'year']]
    )
    n_sps = len(new_sp_df)
    # get initial condition stress period
    ic_sp = (
        sp_df
        .query(f'month == {ic_mon_year[0]} and year == {ic_mon_year[1]}')
        .sp
        .item()
    )
    # load gwf
    gwf = sim.get_model()
    if gwf is None:
        raise ValueError("Model could not be loaded. Please check the simulation workspace.")
    # get model name
    model_name = gwf.name 
    
    # Make a new simulation -----------
    # Simulation-level modules
    print('generating new simulation with modified sp and timeseries info')
    # (MFSIM) Simulation
    params = list(get_parameter_set(sim))
    param_dict = param_dict_from_list(sim, params)
    sim_new = flopy.mf6.MFSimulation(sim_ws = new_sim_ws,**param_dict)
    # # (TDIS) Temporal discretization
    tdis = sim.get_package('tdis')
    params_tdis = list(get_parameter_set(tdis))
    param_dict_tdis = param_dict_from_list(tdis, params_tdis)
    param_dict_tdis['nper'] = n_sps
    param_dict_tdis['start_date_time'] = start_date_time
    param_dict_tdis['perioddata'] = new_sp_df.loc[:, ['perlen', 'nstp', 'tsmult']].to_records(index=False)
    new_pack_dict = {}
    new_pack_dict['tdis'] = flopy.mf6.ModflowTdis(sim_new, **param_dict_tdis)
    
    # Groundwater flow model modules -----------
    # (GWF) Groundwater Flow Model
    gwf_new = flopy.mf6.ModflowGwf(sim_new, modelname=gwf.name)
    # (IMS) Iterative model solution (simulation level, but I'm doing it here because of model names)
    new_pack_dict['ims'] = copy_package(sim, 'ims', sim_new)
    sim_new.register_ims_package(new_pack_dict['ims'], [gwf_new.name]) # type: ignore
    # modules with no time info: 
    for pack_name in ['dis', 'npf']:
        new_pack_dict[pack_name] = copy_package(gwf, pack_name, gwf_new)
    
    # modules with time info:
    # sto
    sto = gwf.get_package('sto')
    param_dict_sto = copy_param_dict(sto)
    param_dict_sto['transient'] = {sp:True for sp in range(n_sps)}
    param_dict_sto['steady_state'] = {sp:False for sp in range(n_sps)}
    new_pack_dict['sto'] = flopy.mf6.ModflowGwfsto(gwf_new, **param_dict_sto)
    # # oc
    if gwf.get_package('oc') is not None:
        oc = gwf.get_package('oc')
        param_dict_oc = copy_param_dict(oc)
        param_dict_oc['saverecord'] = {0: param_dict_oc['saverecord'][0]}
        param_dict_oc['printrecord'] = {0: param_dict_oc['printrecord'][0]}
        new_pack_dict['oc'] = flopy.mf6.ModflowGwfoc(gwf_new, **param_dict_oc)
    
    # ic
    # get hds values from appropriate sp
    headfilename = os.path.join(f'{sim_ws}', f'{model_name}.hds')
    hds = flopy.utils.HeadFile(headfilename, precision='double')
    ic_heads = hds.get_data(kstpkper=(0,ic_sp))
    # grab package and set new strt values
    ic = copy_package(gwf, 'ic', gwf_new)
    ic.strt = ic_heads # type: ignore
    new_pack_dict['ic'] = ic
    
    # boundary files 
    # drn # ghb # riv # wel
    bound_pack_names = []
    for pack_type in ['drn', 'ghb', 'riv', 'wel']:
        packs = gwf.get_package(pack_type)
        if isinstance(packs, list):
            bound_pack_names.extend([pack.package_name for pack in packs])
        else:
            if packs is not None:
                bound_pack_names.append(packs.package_name)
    for pack_name in bound_pack_names:
        new_pack_dict[pack_name] = weight_mean_package(gwf, pack_name, gwf_new, sp_data_lut)
    
    # return results
    return sim_new, new_pack_dict

