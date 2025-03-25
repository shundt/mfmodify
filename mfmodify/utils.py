# IMPORT
import copy
import inspect
import numpy as np
import pandas as pd
import flopy

# VARIABLES

# FUNCTIONS
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
        raise ValueError("Frequency must be 'Y', 'M', 'D', 'h', 'min', 's', 'ms")
    return rounded

def get_sp_data(sim, snap_dates='M'):
    tdis = sim.get_package('tdis')
    start_date_time = pd.to_datetime(tdis.start_date_time.data)
    sp_df = (
        pd.DataFrame(tdis.perioddata.array) 
        .assign(endtime = lambda x: x.perlen.cumsum())
        .assign(starttime = lambda x: [0] + x.endtime[:-1].to_list())
        .assign(start_date_time = lambda x: 
            x.starttime.map(lambda x: start_date_time + pd.Timedelta(days=x)))
        .assign(start_date_time = lambda x: x.start_date_time.dt.round('min')) 
        .assign(start_date = lambda x: pd.to_datetime(x.start_date_time.dt.date))
        .assign(snap_date = lambda x: x.start_date.apply(round_date, args=snap_dates))
        .assign(year = lambda x: x.snap_date.dt.year)
        .assign(month = lambda x: x.snap_date.dt.month)
        .assign(sp = lambda x: x.index)
        .loc[:, ['sp', 'start_date', 'snap_date', 'year', 'month', 'perlen', 
            'nstp', 'tsmult', 'endtime', 'starttime', 'start_date_time']]
    )
    return sp_df

def get_parameter_set(pack, rem_att_set=set([])):
    # get list of attributes
    pack_att_set = set(vars(pack).keys())
    # Get a list of parameters for instantiating the class
    constructor_params = inspect.signature(pack.__init__).parameters
    param_name_set = set(constructor_params.keys())
    # get list of attributes to extract and store in dictionary to unpack on instantiation
    attribute_trans_list = (
        pack_att_set.intersection(param_name_set) - rem_att_set)
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

def get_objects_from_pack(pack, attribute_name):
    object_list = []
    for i in range(1000):
        try:
            # Use getattr to access the attribute
            obj = getattr(pack, attribute_name)[i]  
            object_list.append(obj)
        except ValueError:
            break
    return object_list

def get_ts_objects(pack):
    return get_objects_from_pack(pack, 'ts')

def get_obs_objects(pack):
    return get_objects_from_pack(pack, 'obs')

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
    # get manual parameters
    for att, val in manual_params.items():
        pack_param_dict[att] = val
    # instantiate package
    pack_new = pack_class(sim_or_gwf_new, **pack_param_dict)
    return pack_new

def get_gwf_package_df(gwf):
    gwf_package_df = (
        pd.DataFrame(gwf.name_file.packages.array)
        .assign(ftype = lambda x: [ft[:-1].lower() for ft in x.ftype])
        .assign(dummy = 1)
        # .sort_values(['ftype', 'pname'], ascending=True)
        .assign(paknum = lambda x: x.groupby('ftype').dummy.transform('cumsum'))
        .assign(paknum = lambda x: [nm if nm!=0 else '' for nm in x.paknum])
        .assign(pakname = lambda x: [f'{tp}{nm}' for nm, tp in zip(x.paknum, x.ftype)])
        .drop(['dummy', 'paknum'], axis=1)
    )
    return gwf_package_df

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

def lst_df_from_gwf_long(gwf):
    # get gwf package info
    gwf_package_df = get_gwf_package_df(gwf)
    pakname_lut = gwf_package_df.set_index('pakname').pname
    # get wide lst df
    lst_df = lst_df_from_gwf(gwf)
    # join with paknam lut and format
    lst_df_long_all = (
        lst_df
        .melt(id_vars = ['sp', 'year', 'month'], var_name='pak', value_name='rate')
        .assign(pakname = lambda x: [x.split('_')[0].lower() for x in x.pak])
        .join(pakname_lut, on='pakname')
        .assign(pname = lambda x: x.pname.fillna(x.pakname))
        .assign(direction = lambda x: [x.split('_')[-1].lower() for x in x.pak])
        .assign(bcompname = lambda x: 
            [f'{nm}-{di}' for nm,di in zip(x.pname, x.direction)])
    )
    # find list of components that have values
    comp_keep = (
        lst_df_long_all
        .groupby('bcompname')
        .sum()
        .query('rate!=0')
        .index
        .tolist()
    )
    # query by those components
    lst_df_long = (
        lst_df_long_all
        .query(f'bcompname in {comp_keep}')
        .drop(['pakname'], axis=1)
    )
    return lst_df_long

def get_final_total_lst(gwf):
    # get the simulation
    sim = gwf.simulation
    # get the list
    lst = gwf.output.list()
    # get cumulative volumes as df
    lst_df_cumulative = lst.get_dataframes()[1]
    # do all the manipulation
    lst_final = (
        lst_df_cumulative
        .drop(['IN-OUT', 'PERCENT_DISCREPANCY'], axis=1)
        .iloc[-1, :]
        .to_frame()
        .assign(total_volume = lambda x: x.iloc[:, 0])
        .reset_index()
        .rename(columns={'index': 'component'})
        .assign(direction = lambda x: [co.split('_')[-1].lower() for co in x.component])
        .assign(btype = lambda x: [co[:3].lower() for co in x.component])
        .groupby(['btype', 'direction'])
        .total_volume
        .sum()
        .reset_index()
        .sort_values(['direction', 'btype'])
        .set_index('direction')
    )
    return lst_final

def copy_empty_sim(sim, sim_ws):
    params = list(get_parameter_set(sim))
    param_dict = param_dict_from_list(sim, params)
    return flopy.mf6.MFSimulation(sim_ws = sim_ws,**param_dict)

def copy_sim(sim, sim_ws):
    sim_new = copy.deepcopy(sim)
    sim_new.set_sim_path(sim_ws)
    return sim_new

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
            **hdobs_continuous
        }
        hdobs_params['continuous'] = new_cont_dict
        hdobs_params['digits'] = digits
        hdobs = flopy.mf6.modflow.mfutlobs.ModflowUtlobs(
            gwf,
            **hdobs_params
        )
    return hdobs

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