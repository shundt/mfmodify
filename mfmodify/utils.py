# IMPORT
import copy
import inspect
import numpy as np
import pandas as pd
import flopy

# VARIABLES

# FUNCTIONS
def round_date(ts, frequency='M'):
    """
    Round a timestamp to the nearest specified frequency.

    This function rounds a given timestamp to the nearest month, year, or other 
    specified time frequency. For monthly and yearly frequencies, the midpoint 
    between the current and next period is used to determine rounding. For other 
    frequencies, the timestamp is rounded directly.

    Parameters:
    ts (pd.Timestamp): The timestamp to be rounded.
    frequency (str): The frequency to round to. Supported values are:
        - 'M': Month
        - 'Y': Year
        - 'D': Day
        - 'h': Hour
        - 'min': Minute
        - 's': Second
        - 'ms': Millisecond

    Returns:
    pd.Timestamp: The rounded timestamp.

    Raises:
    ValueError: If the frequency is not one of the supported values.
    """
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
    """
    Generate a stress period data table for a simulation.

    This function extracts stress period data from a MODFLOW 6 simulation and 
    organizes it into a DataFrame. It includes information such as start and 
    end times, stress period lengths, and optionally rounded or snapped dates.

    Parameters:
    sim (flopy.mf6.MFSimulation): The MODFLOW 6 simulation object.
    snap_dates (str, optional): The frequency to which dates should be rounded 
        or snapped. Supported values are:
        - 'M': Month
        - 'Y': Year
        - 'D': Day
        - 'h': Hour
        - 'min': Minute
        - 's': Second
        - 'ms': Millisecond
        Defaults to 'M' (Month).

    Returns:
    pd.DataFrame: A DataFrame containing stress period data with the following columns:
        - 'sp': Stress period index.
        - 'start_date': The start date of the stress period.
        - 'snap_date': The snapped or rounded date.
        - 'year': The year of the snapped date.
        - 'month': The month of the snapped date.
        - 'perlen': The length of the stress period.
        - 'nstp': The number of time steps in the stress period.
        - 'tsmult': The time step multiplier.
        - 'endtime': The cumulative end time of the stress period.
        - 'starttime': The cumulative start time of the stress period.
        - 'start_date_time': The exact start date and time of the stress period.
    """
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
    """
    Retrieve a set of parameters for instantiating a package.

    This function identifies the attributes of a package that can be used as 
    parameters for instantiating a new instance of the same package. It excludes 
    attributes specified in the `rem_att_set`.

    Parameters:
    pack (object): The package object from which parameters are extracted.
    rem_att_set (set, optional): A set of attribute names to exclude from the 
        parameter set. Defaults to an empty set.

    Returns:
    set: A set of attribute names that can be used as parameters for instantiating 
        the package.
    """
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
    """
    Create a dictionary of parameters from a list of attributes.

    This function extracts the values of specified attributes from a package 
    object and organizes them into a dictionary. The dictionary can be used 
    to instantiate a new instance of the package.

    Parameters:
    pack (object): The package object from which attributes are extracted.
    param_list (list): A list of attribute names to extract from the package.

    Returns:
    dict: A dictionary containing the extracted attributes and their values.
    """
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
    """
    Retrieve a list of objects from a package attribute.

    This function extracts objects from a specified attribute of a package. 
    It iterates through the attribute, collecting objects until a `ValueError` 
    is encountered, indicating the end of the collection.

    This function is used by `get_ts_objects` and `get_obs_objects` to 
    retrieve time series (`ts`) and observation (`obs`) objects, respectively.

    Parameters:
    pack (object): The package object from which the attribute is accessed.
    attribute_name (str): The name of the attribute to retrieve objects from.

    Returns:
    list: A list of objects retrieved from the specified attribute.
    """
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
    """
    Create a dictionary of parameters for copying a package.

    This function extracts the parameters of a package object and organizes 
    them into a dictionary. The dictionary can be used to instantiate a new 
    instance of the package with the same parameters.

    Parameters:
    pack (object): The package object from which parameters are extracted.

    Returns:
    dict: A dictionary containing the extracted parameters and their values.
    """
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
    """
    Copy a package from one simulation or model to another.

    This function copies a package from an original simulation or groundwater 
    flow model (GWF) to a new simulation or model. It allows for manual 
    parameter overrides during the copying process.

    Parameters:
    sim_or_gwf_orig (flopy.mf6.MFSimulation or flopy.mf6.ModflowGwf): The original 
        simulation or GWF model containing the package to be copied.
    pack_name (str): The name of the package to be copied.
    sim_or_gwf_new (flopy.mf6.MFSimulation or flopy.mf6.ModflowGwf): The new 
        simulation or GWF model to which the package will be added.
    manual_params (dict, optional): A dictionary of parameters to override or 
        add to the copied package. Defaults to an empty dictionary.

    Returns:
    flopy.mf6.ModflowGwfPackage: The copied package associated with the new 
        simulation or model.
    """
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
    """
    Generate a DataFrame of package information for a GWF model.

    This function extracts information about the packages in a groundwater flow 
    (GWF) model and organizes it into a DataFrame. The DataFrame includes details 
    such as package types, names, and formatted package identifiers.

    Parameters:
    gwf (flopy.mf6.ModflowGwf): The GWF model object.

    Returns:
    pd.DataFrame: A DataFrame containing package information with the following columns:
        - 'ftype': The package type (e.g., 'dis', 'npf').
        - 'pname': The package name.
        - 'pakname': A formatted package identifier combining the type and number.
    """
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
    """
    Generate a DataFrame from the listing file of a GWF model.

    This function extracts data from the listing file of a groundwater flow 
    (GWF) model and organizes it into a DataFrame. The DataFrame includes 
    information about stress periods, recharge, discharge, and storage changes.

    Parameters:
    gwf (flopy.mf6.ModflowGwf): The GWF model object.

    Returns:
    pd.DataFrame: A DataFrame containing listing file data with the following columns:
        - 'sp': Stress period index.
        - 'year': The year of the stress period.
        - 'month': The month of the stress period.
        - 'total_recharge': Total recharge during the stress period.
        - 'total_discharge': Total discharge during the stress period.
        - 'net_storage': Net storage change during the stress period.
    """
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
    """
    Generate an annual summary of recharge, discharge, and storage changes for a GWF model.

    This function aggregates data from the listing file of a groundwater flow 
    (GWF) model to produce an annual summary. It calculates total recharge, 
    total discharge, and net storage changes for each complete year.

    Parameters:
    gwf (flopy.mf6.ModflowGwf): The GWF model object.

    Returns:
    pd.DataFrame: A DataFrame containing the annual summary with the following columns:
        - 'total_recharge': Total recharge during the year.
        - 'total_discharge': Total discharge during the year.
        - 'net_storage': Net storage change during the year.
    """
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
    """
    Generate a long-format DataFrame from the listing file of a GWF model.

    This function converts the listing file data of a groundwater flow (GWF) model 
    into a long-format DataFrame. It includes information about stress periods, 
    recharge, discharge, and storage changes, along with package and component details.

    Parameters:
    gwf (flopy.mf6.ModflowGwf): The GWF model object.

    Returns:
    pd.DataFrame: A long-format DataFrame containing listing file data with the following columns:
        - 'sp': Stress period index.
        - 'year': The year of the stress period.
        - 'month': The month of the stress period.
        - 'pak': The package name.
        - 'rate': The rate associated with the package.
        - 'pname': The formatted package name.
        - 'direction': The flow direction (e.g., 'in', 'out').
        - 'bcompname': The combined component name (e.g., 'package-direction').
    """
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
    """
    Generate a summary of final total volumes from the listing file of a GWF model.

    This function extracts cumulative volume data from the listing file of a 
    groundwater flow (GWF) model and organizes it into a summary DataFrame. 
    It includes information about the total volume for each component, grouped 
    by flow direction and component type.

    Parameters:
    gwf (flopy.mf6.ModflowGwf): The GWF model object.

    Returns:
    pd.DataFrame: A DataFrame containing the final total volumes with the following columns:
        - 'btype': The component type (e.g., 'wel', 'riv').
        - 'direction': The flow direction (e.g., 'in', 'out').
        - 'total_volume': The total volume for each component and direction.
    """
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
    """
    Add or update a head observations (hdobs) package in a GWF model.

    This function adds a new head observations (hdobs) package to a groundwater 
    flow (GWF) model or updates an existing one (since mf6 only allows one). 
    It allows for the inclusion of new continuous observation data and adjusts 
    the package parameters accordingly.

    Parameters:
    gwf (flopy.mf6.ModflowGwf): The GWF model object.
    hdobs_continuous (dict): A dictionary containing continuous observation data 
        to be added or updated in the hdobs package.
    digits (int, optional): The number of digits to use for formatting observation 
        output. Defaults to 5.

    Returns:
    flopy.mf6.modflow.mfutlobs.ModflowUtlobs: The updated or newly created hdobs package.
    """
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
    """
    Generate a DataFrame representing the idomain array of a GWF model.

    This function extracts the idomain array from a groundwater flow (GWF) model 
    and organizes it into a DataFrame. The DataFrame includes cell IDs, idomain 
    values, and additional spatial information such as layer, row, and column 
    indices, depending on the grid type.

    Parameters:
    gwf (flopy.mf6.ModflowGwf): The GWF model object.

    Returns:
    pd.DataFrame: A DataFrame containing the idomain data with the following columns:
        - 'cellid': The cell ID.
        - 'idomain': The idomain value for the cell.
        - 'nodeid': The node ID for the cell.
        - 'layer' (optional): The layer index (for 3D grids).
        - 'row' (optional): The row index (for structured grids).
        - 'column' (optional): The column index (for structured grids).
        - 'icell2d' (optional): The cell2d index (for vertex grids).
    """
    # get dis and modelgrid
    dis = gwf.get_package('dis')
    grid = gwf.modelgrid
    # get idomain and flatten it
    idomain = dis.idomain.data
    idomain_flat = idomain.ravel()
    indices_flat_npint = [index.ravel() for index in np.indices(idomain.shape)]
    indices_flat = [getattr(x, 'tolist', lambda: x)() for x in indices_flat_npint]
    # get cellids
    cellids = list(zip(*indices_flat))
    # make a dataframe and add nodeid
    idomain_df = (
        pd
        .DataFrame({'cellid': cellids, 'idomain': idomain_flat}) 
        .sort_values('cellid')
        .assign(nodeid = lambda x: range(x.shape[0]))
    )
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

def pak_sp_prop_to_array(pak, prop, grid_relate, sp=0):
    """
    Convert a package stress period property to an array.

    This function extracts a specified property from a package's stress period 
    data and converts it into an array format. The function supports structured, 
    unstructured, and vertex grids.

    Parameters:
    pak (flopy.mf6.ModflowGwfPackage): The package object containing stress period data.
    prop (str): The property name to extract from the stress period data.
    grid_relate (pd.DataFrame): A DataFrame relating grid cell IDs to the model grid.
    sp (int, optional): The stress period index to extract data for. Defaults to 0.

    Returns:
    np.ndarray: An array representing the specified property for the given stress period.
    """
    prop_ra = pak.stress_period_data.data[sp]
    prop_df = (
        pd
        .DataFrame(prop_ra)
        .loc[:, ['cellid', prop]]
        .groupby('cellid')
        .mean()
    )
    # figure out grid type
    if isinstance(prop_df.index[0], int): # unstructured
        # unstructured
        pass
    elif isinstance(prop_df.index[0], tuple): 
        if len(prop_df.index[0]) == 3: # structured
            prop_all_cells = (
                prop_df
                .reindex(grid_relate.index)
                .assign(layer = lambda x: [i[0] for i in x.index])
                .assign(row = lambda x: [i[1] for i in x.index])
                .assign(column = lambda x: [i[2] for i in x.index])
            )
            rows = prop_all_cells.row.unique()
            columns = prop_all_cells.column.unique()
            array2d_list = []
            for layer in prop_all_cells.layer.unique():
                array2d = (
                    prop_all_cells
                    .query(f'layer == {layer}')
                    .pivot_table(index='row', columns='column', values=prop)
                    .reindex(index=rows, columns=columns)
                )
                if len(array2d) == 0:
                    array2d = np.empty((prop_all_cells.row.max()+1, prop_all_cells.column.max()+1))
                array2d_list.append(array2d)
            array_final = np.stack(array2d_list, axis=0)
        elif len(prop_df.index[0]) == 2: # vertex grid
            prop_all_cells = (
                prop_df
                .reindex(grid_relate.cellid_disv)
                .assign(layer = lambda x: [i[0] for i in x.index])
                .assign(cell2d = lambda x: [i[1] for i in x.index])
            )
            cell2ds = prop_all_cells.cell2d.unique()
            array1d_list = []
            for layer in prop_all_cells.layer.unique():
                array1d = (
                    prop_all_cells
                    .query(f'layer == {layer}')
                    .set_index('cell2d')
                    .cond
                    .reindex(cell2ds)
                )
                array1d_list.append(array1d)
            array_final = np.stack(array1d_list, axis=0)
    return array_final