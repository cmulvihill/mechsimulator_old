"""
Parses Excel files containing experimental data to obtain exp_set dictionaries
"""

import pandas as pd
import numpy as np
from mechsimulator.parser import exp_checker
from mechsimulator.parser import util
from mechsimulator.parser import spc

# Allowed physical quantities (keys) and:
# (1) corresponding allowed units
# (2) conversions factors to get to the internal units
# (3) full name and desired units of the quantity (for printing purposes)
# Note: the first units given for each quantity are the internal units
ALLOWED_UNITS = {
    'temp':         (('K', 'C', 'F', 'R'),
                     (1, None, None, 5.55555e-1),
                     'Temperature'),
    'pressure':     (('atm', 'bar', 'kPa', 'Pa', 'MPa'),
                     (1, 9.86923e-1, 9.86923e-3, 9.86923e-6, 9.86923e0),
                     'Pressure'),
    'time':         (('s', 'ms', 'micros'),
                     (1, 1.0e-3, 1.0e-6),
                     'Time'),
    'conc':         (('X', 'ppm', '%'),
                     (1, 1.0e-6, 1.0e-2),
                     'Mole fraction'),
    'length':       (('m', 'cm', 'mm'),
                     (1, 1.0e-2, 1.0e-3),
                     'Length'),
    'area':         (('m2', 'cm2', 'mm2'),
                     (1, 1.0e-4, 1.0e-6),
                     'Area'),
    'vol':          (('m3', 'cm3', 'mm3'),
                     (1, 1.0e-6, 1.0e-9),
                     'Volume'),
    'abs_coeff':    (('cm-1atm-1',),
                     (1,),
                     'Absorption coefficient'),
    'abs':          (('%', 'fraction'),
                     (1, 100),
                     'Absorption'),
    'dpdt':         (('%/ms',),
                     (1,),
                     'dP/dt'),
    'mdot':         (('kg/s', 'g/s'),
                     (1, 1.0e-3,),
                     'Mass flow rate')
}

# Mappings for alternate names (keys) to their physical quantities (values)
ALTERNATE_NAMES = {
    'path_length':  'length',
    'end_time':     'time',
    'res_time':     'time',
    'idt':          'time',
    'timestep':     'time',
}

# Number of columns in the experimental sheets
NUM_EXP_CLMS = 7


def load_exp_set(exp_filename):
    """ Reads an Excel file that contains data on an experiment set and returns
        an exp_set object

        :param exp_filename: path to .xlsx file describing an experiment set
        :type exp_filename: str
        :return exp_set: description of a set of experiments
        :rtype: dct
    """

    # Use Pandas to read the Excel file
    sheet_names = pd.ExcelFile(exp_filename).sheet_names
    assert 'info' in sheet_names, (
        f'There is no sheet called "info" in the Excel file {exp_filename}.')
    del sheet_names[sheet_names.index('info')]  # remove info
    if isinstance(sheet_names, str):  # in the case of only one exp sheet
        sheet_names = [sheet_names]

    # Read the info sheet and create the exp_set
    print(f'Reading the "info" sheet in Excel file {exp_filename}...')
    info = pd.read_excel(exp_filename, sheet_name='info', header=None,
                         skiprows=1)
    exp_set = read_info_sheet(info)

    # Read the experiment sheets and create the exp_objs
    exp_objs = []
    for sheet_name in sheet_names:
        print(f'Reading sheet {sheet_name} in Excel file {exp_filename}...')
        df = pd.read_excel(exp_filename, sheet_name=sheet_name, header=None,
                           skiprows=1)
        exp_obj = read_exp_sheet(df, exp_set)
        exp_objs.append(exp_obj)

    exp_set['exp_objs'] = exp_objs

    # Run some checks
    exp_checker.chk_exp_set(exp_set)
    for exp_obj in exp_objs:  # important to check exp_objs AFTER exp_set
        exp_checker.chk_exp_obj(exp_obj, exp_set)

    # Create the exp_data and xdata arrays
    exp_ydata, exp_xdata = get_exp_data(exp_set)
    exp_set['overall']['exp_ydata'] = exp_ydata
    exp_set['overall']['exp_xdata'] = exp_xdata

    # print(exp_ydata)
    # print(exp_xdata)

    return exp_set


def read_info_sheet(df):
    """ Reads the Excel 'info' sheet describing information on an experiment set

        :param df: a Pandas dataframe describing the contents of the Excel sheet
        :type df: Pandas dataframe
        :return: exp_set: description of a set of experiments; used for spc list
        :rtype exp_set: dct
    """

    assert df.shape[1] >= 4, f'"info" sheets must have at least 4 columns'
    exp_set = {'overall': {}, 'plot': {}, 'plot_format': {}, 'mix': {},
               'spc': {}, 'exp_objs': []}
    for _, row in df.iterrows():
        assert len(row) >= 4, f'The row {list(row)} has less than 4 entries'
        group, param, raw_val, units = row[0], row[1], row[2], row[3]
        assert group in exp_checker.ALLOWED_INFO_GROUPS, (
            f"The group '{group}' is not allowed in the info sheet. Options "
            f"are: {exp_checker.ALLOWED_INFO_GROUPS}")
        # Perform different actions based on the group
        if group == 'overall':
            exp_set[group][param] = raw_val  # no conversions required
        elif group == 'mix':
            conv_val = convert_units(raw_val, 'conc', units)
            exp_set['mix'][param] = conv_val  # param is the spc name
        elif group == 'spc':
            exp_set[group][param] = read_spc_row(row)
        elif group == 'plot':
            # If on the plotting variable info
            if param in ('start', 'end', 'inc'):
                # Check that the plotting variable has already been defined
                assert 'variable' in exp_set['plot'].keys(), (
                    f'The parameter "variable" must be defined before {param}')
                # Use the plotting variable quantity as the quantity
                fake_param = exp_set['plot']['variable']  # e.g., 'temp'
                conv_val = convert_units(raw_val, fake_param, units)
            else:
                conv_val = convert_units(raw_val, param, units)
            # If there is a comma in the input field, split the str around it
            if isinstance(conv_val, str) and ',' in conv_val:
                conv_val = [entry.strip() for entry in conv_val.split(',')]
            # If there was no comma, convert the single value from str to lst
            # for a few special parameters
            elif param in ('idt_target', 'idt_method', 'wavelength'):
                conv_val = [conv_val]
            exp_set[group][param] = conv_val  # note: don't use fake_param
        elif group == 'plot_format':
            if param in ('xlim', 'ylim'):
                conv_val = [float(entry) for entry in raw_val.split(',')]
            elif param == 'rows_cols':
                conv_val = [int(entry) for entry in raw_val.split(',')]
            elif param == 'omit_targets':
                conv_val = [entry.strip() for entry in raw_val.split(',')]
            elif param == 'plot_points':
                if raw_val == 'yes':
                    conv_val = True
                elif raw_val == 'no':
                    conv_val = False
                else:
                    conv_val = raw_val
            else:
                conv_val = raw_val
            exp_set[group][param] = conv_val
        elif group == 'ignore_row':
            pass

    return exp_set


def read_exp_sheet(df, exp_set):
    """ Reads an Excel sheet describing a single experiment

        :param df: a Pandas dataframe describing the contents of the Excel sheet
        :type df: Pandas dataframe
        :param: exp_set: description of a set of experiments; used for spc list
        :type exp_set: dct
        :return exp_obj: description of a single experiment
        :rtype: dct
    """

    assert df.shape[1] >= NUM_EXP_CLMS, (
        f'Experimental sheets must have at least {NUM_EXP_CLMS} columns')
    exp_obj = {'overall': {}, 'conds': {}, 'mix': {}, 'result': {}}

    # Loop over each row
    for _, row in df.iterrows():
        group = row[0]
        if group == 'overall':
            param = row[1]
            exp_obj[group][param] = row[2]
        elif group == 'ignore_row':
            pass
        else:
            val_tuple, param, extra_params = read_row(row, exp_set)
            if extra_params == []:  # the usual case: no commas in 'param' field
                if param in ('idt', 'abs', 'emis'):  # if a special case
                    assert exp_obj[group].get(param) is None, (
                        f"The '{param}' field has already been defined! Use "
                        f"a comma-separated integer to specify new instances.")
                    exp_obj[group][param] = {1: val_tuple}  # save with a 1
                else:  # otherwise, just save the value (most common case)
                    exp_obj[group][param] = val_tuple
            else:  # if there were commas in the parameter field
                observable_num = int(extra_params[0])  # int starting with 1
                # Create dct if needed
                if exp_obj[group].get(param) is None:
                    exp_obj[group][param] = {}
                # If on a time array of absorption data or IDT data
                if param in ('abs', 'emis', 'idt'):
                    exp_obj[group][param][observable_num] = val_tuple
                # If on an absorption coefficient
                if param == 'abs_coeff':
                    # If both a wavelength_num and spc were given
                    if len(extra_params) == 2:
                        spc = extra_params[1]
                    # Otherwise, if only spc was given, use 1 for observable_num
                    elif len(extra_params) == 1:
                        observable_num = 1
                        spc = extra_params[0]
                    # Create dct if needed
                    if exp_obj[group][param].get(observable_num) is None:
                        exp_obj[group][param][observable_num] = {}
                    exp_obj[group][param][observable_num][spc] = val_tuple

    return exp_obj


def read_row(row, exp_set):
    """ Reads a row in an experimental Excel sheet

        Note: does not read the 'spc' row (see the read_spc_row function) nor
        the 'overall' row

        :param row: an object representing a single row of an Excel sheet
        :type row: Pandas row object obtained from a dataframe
        :param exp_set:
        :type exp_set:
        :return val_tuple: a tuple describing the value, the lower and upper
            bounds on the value, and the bound type ('percent' or 'value'). Note
            that val can either be a float or a Numpy array of shape (num_vals,)
        :rtype: tuple (val, low_bnd, upp_bnd, bnd_type)
    """

    # Read data
    param, units = row[1], row[3]
    low_bnd, upp_bnd, bnd_type = row[4], row[5], row[6],
    all_spcs = exp_set['spc'].keys()  # used to check spcs names

    # Check for extra entries in the parameter column, indicated by comma(s)
    extra_params = []
    if isinstance(param, str):
        split_param = param.split(',')
        if len(split_param) > 1:  # if there were any commas
            param = split_param[0]
            extra_params.extend(split_param[1:])

    # If the "value" column is blank, assume the row is a time series
    chkd_val_clmn = util.chk_entry(row[2])
    if chkd_val_clmn is None:
        num_points = len(row) - NUM_EXP_CLMS  # num of time-resolved datapoints
        val = np.zeros(num_points)
        for idx, item in enumerate(row[NUM_EXP_CLMS:]):
            val[idx] = item
    else:
        val = row[2]  # otherwise, it's a single value

    # Convert the units of the value(s)
    if param in all_spcs:  # if param is a spcs name
        val = convert_units(val, 'conc', units)
    else:
        val = convert_units(val, param, units)
        
    # Check if any of the values were NaN (i.e., blank in the Excel file) or
    # '-'; if so, return None. If not, leave unchanged
    low_bnd = util.chk_entry(low_bnd)
    upp_bnd = util.chk_entry(upp_bnd)
    bnd_type = util.chk_entry(bnd_type)

    val_tuple = (val, low_bnd, upp_bnd, bnd_type)

    return val_tuple, param, extra_params


def read_spc_row(row):
    """ Reads a 'spc' row in an experimental Excel sheet (i.e., not the 'info'
        sheet)

        :param row: an object representing a single row of an Excel sheet
        :type row: Pandas row object obtained from a dataframe
        :return spc_dct: description of a single species
        :rtype: dct
    """

    assert len(row) >= 6, (
        f"'spc' row {list(row)} is too short. Should have smiles, multiplicity,"
        f" charge, and excited flag given.")
    smiles, mult, charge, exc_flag = row[2], row[3], row[4], row[5]
    inchi = spc.to_inchi(spc.from_smiles(smiles))
    spc_dct = {
        'inchi': inchi,
        'smiles': smiles,
        'mult': mult,
        'charge': charge,
        'exc_flag': exc_flag}

    return spc_dct


def convert_units(val, quant, units):
    """ Converts a value of some physical quantity in specified units to
        desired standard units

        Note: the desired units are found as the first entry in each value of
            the ALLOWED_UNITS

        Note: some quantities have alternate names. For example, "end_time" is
            simply a time. These alternate names are defined by the variable
            ALTERNATE_NAMES.

        :param val: the value of some physical dimension (e.g., 1.04)
        :type val: float or Numpy of shape (num_vals,)
        :param quant: the physical quantity that is quantified by val
            (e.g., 'pressure')
        :type quant: str
        :param units: the units used to quantify val (e.g., 'bar')
        :type units: str
        :return conv_val: val that has been converted to desired internal units
        :rtype: float
    """

    # If the units are blank or '-', don't convert
    chkd_entry = util.chk_entry(units)
    if chkd_entry is None:
        # Check if the quantity should have units assigned
        if val != 'bal':  # exception: 'bal' is a 'conc' but requires no units
            assert quant not in ALLOWED_UNITS.keys() and quant not in \
                ALTERNATE_NAMES.keys(), f"The quantity '{quant}' requires units"
        conv_val = val  # no conversion if units are blank or '-'

    # Otherwise, convert the units
    else:
        # Check if an alternate name exists for quant
        if quant in ALTERNATE_NAMES.keys():
            quant = ALTERNATE_NAMES[quant]  # e.g., changes "end_time" to "time"
        # Check that the quantity and the units are allowed
        allowed_quants = tuple(ALLOWED_UNITS.keys())
        assert quant in allowed_quants, (
            f'Quantity "{quant}" is not allowed. Options are: {allowed_quants}')
        allowed_units = ALLOWED_UNITS.get(quant)[0]
        assert units in allowed_units, (
            f'Units "{units}" are not allowed for the quantity "{quant}".'
            f' Options are: {allowed_units}')

        # Get the conversion factor
        conv_factors = ALLOWED_UNITS.get(quant)[1]
        units_idx = allowed_units.index(units)
        conv_factor = conv_factors[units_idx]

        # Check for special cases
        if quant == 'temp' and units == 'C':  # Celsius
            conv_val = val + 273.15
        elif quant == 'temp' and units == 'F':  # Fahrenheit
            conv_val = (val - 32) / 1.8 + 273.15
        # Otherwise, use the conversion factor
        else:
            conv_val = val * conv_factor

    return conv_val


def get_exp_data(exp_set):

    def get_unique_times(exp_set):

        all_times = []
        exp_objs = exp_set['exp_objs']
        for exp_obj in exp_objs:
            all_times.extend(exp_obj['result']['time'][0])
        unique_times = sorted(list(set(all_times)))

        return unique_times

    def interp_single_array(sing_ydata, exp_xdata, unique_xdata):
        """

        :param sing_ydata:
        :param exp_xdata:
        :param unique_xdata:
        :return:
        """

        interp_ydata = np.ndarray((len(unique_xdata)), dtype=float)
        # Loop over every unique x point
        for unique_idx, unique_xdatum in enumerate(unique_xdata):
            # If the unique x point is in the exp x data, store matching y data
            if unique_xdatum in exp_xdata:
                exp_idx = np.where(exp_xdata == unique_xdatum)
                interp_ydata[unique_idx] = sing_ydata[exp_idx]
            # Otherwise, just store as NaN
            else:
                interp_ydata[unique_idx] = np.nan

        return interp_ydata

    # Load some initial parameters
    ndims = util.get_exp_dims(exp_set)
    meas_type = exp_set['overall']['meas_type']
    plot_var = exp_set['plot']['variable']  # the variable for the exp_set
    exp_objs = exp_set['exp_objs']
    num_conds = len(exp_objs)

    # Get exp_xdata and exp_ydata based on measurement type
    if meas_type == 'outlet':
        # Get the xdata
        exp_xdata = np.ndarray(len(exp_objs))
        for cond_idx, exp_obj in enumerate(exp_objs):
            exp_xdata[cond_idx] = exp_obj['conds'][plot_var][0]

        # Get the ydata
        spcs = exp_set['spc'].keys()
        num_targets = len(spcs)
        exp_ydata = np.ndarray((num_conds, num_targets))
        for cond_idx, exp_obj in enumerate(exp_objs):
            result = exp_obj['result']
            for spc_idx, spc in enumerate(spcs):
                if spc in result:
                    # If the entry is a Numpy array, it was empty in the Excel
                    # sheet; set to NaN in this case
                    if isinstance(result[spc][0], np.ndarray):
                        exp_ydata[cond_idx, spc_idx] = np.nan
                    else:
                        exp_ydata[cond_idx, spc_idx] = result[spc][0]
                else:
                    exp_ydata[cond_idx, spc_idx] = np.nan

    elif meas_type == 'idt':
        # Get the xdata
        exp_xdata = np.ndarray(len(exp_objs))
        for cond_idx, exp_obj in enumerate(exp_objs):
            exp_xdata[cond_idx] = exp_obj['conds'][plot_var][0]

        # Get the ydata
        num_targets = exp_checker.get_num_idts(exp_set)
        exp_ydata = np.ndarray((num_conds, num_targets))
        for cond_idx, exp_obj in enumerate(exp_objs):
            result = exp_obj['result']['idt']
            for idt_idx in range(num_targets):
                if idt_idx + 1 in result:  # +1 b/c idt #s in Excel are 1-index
                    exp_ydata[cond_idx, idt_idx] = result[idt_idx + 1][0]
                else:
                    exp_ydata[cond_idx, idt_idx] = np.nan

    elif meas_type == 'conc':
        # Get the xdata
        exp_xdata = get_unique_times(exp_set)
        num_times = len(exp_xdata)

        # Get the ydata
        spcs = exp_set['spc'].keys()
        num_targets = len(spcs)
        exp_ydata = np.ndarray((num_conds, num_targets, num_times))
        for cond_idx, exp_obj in enumerate(exp_objs):
            result = exp_obj['result']
            exp_times = result['time'][0]
            for spc_idx, spc in enumerate(spcs):
                if spc in result:
                    sing_spc = result[spc][0]
                    exp_ydata[cond_idx, spc_idx] = interp_single_array(
                        sing_spc, exp_times, exp_xdata)
                else:
                    exp_ydata[cond_idx, spc_idx] = np.nan
    else:
        raise NotImplementedError(f"meas_type '{meas_type}' is not working")

    exp_xdata = np.array(exp_xdata)  # convert list to numpy array

    return exp_ydata, exp_xdata
