"""
Parses Excel files containing experimental data to obtain exp_set dictionaries
"""

import pandas as pd
import numpy as np
from mechsimulator.parser import util
from mechsimulator.parser import rdkit_

# Allowed physical quantities (keys) and:
# (1) corresponding allowed units
# (2) conversions factors to get to the desired units
# (3) full name and desired units of the quantity (for printing purposes)
ALLOWED_UNITS = {
    'temp':         (('K', 'C', 'F', 'R'),
                     (1, None, None, 5.55555e-1),
                     'Temperature (K)'),
    'pressure':     (('atm', 'bar', 'kPa', 'Pa', 'MPa'),
                     (1, 9.86923e-1, 9.86923e-3, 9.86923e-6, 9.86923e0),
                     'Pressure (atm)'),
    'time':         (('s', 'ms', 'micros'),
                     (1, 1.0e-3, 1.0e-6),
                     'Time (s)'),
    'conc':         (('X', 'ppm', '%'),
                     (1, 1.0e-6, 1.0e-2),
                     'Concentration (mole fraction)'),
    'length':       (('m', 'cm', 'mm'),
                     (1, 1.0e-2, 1.0e-3),
                     'Length (m)'),
    'area':         (('m2', 'cm2', 'mm2'),
                     (1, 1.0e-4, 1.0e-6),
                     'Area (m^2)'),
    'vol':          (('m3', 'cm3', 'mm3'),
                     (1, 1.0e-6, 1.0e-9),
                     'Volume (m^3)'),
    'abs_coeff':    (('cm-1atm-1',),
                     (1,),
                     'Absorption coefficient (cm^-1.atm^-1)'),
    'abs':          (('%', 'fraction'),
                     (1, 100),
                     'Absorption (%)'),
    'dpdt':         (('%/ms',),
                     (1,),
                     'dP/dt (%/ms)'),
    'mdot':         (('kg/s', 'g/s'),
                     (1, 1.0e-3,),
                     'Mass flow rate (kg/s)')
}

# Mappings for alternate names (keys) to their physical quantities (values)
ALTERNATE_NAMES = {
    'path_length':  'length',
    'end_time':     'time',
    'res_time':     'time',
    'idt':          'time',
    'timestep':     'time',
}

# Allowed groups (keys) and their required params (values) in the 'info' sheet
ALLOWED_INFO_GROUPS = {
    'overall':      ('set_id', 'source', 'description', 'reac_type',
                     'meas_type', 'num_exps', 'version'),
    'plot':         (),  # in some cases, none are required, so leaving blank
    'mix':          (),
    'spc':          (),
    'exp_objs':     (),
    'ignore_row':   (),
}

# Allowed reactor types (keys) and:
# (1) required inputs in the 'conds' field of the exp_sheets and the 'plot'
#       field of the 'info' sheet
# (2) optional inputs in the 'conds' field of the exp_sheets and the 'plot'
#       field of the 'info' sheet
# In the required inputs, a tuple indicates that only one of those values is
# required (e.g., either 'res_time' or 'mdot' is required for a JSR).
ALLOWED_REAC_TYPES = {
    'st':           (('temp', 'pressure', 'end_time'),
                     ('dpdt',)),
    'pfr':          (('temp', 'pressure', 'length', 'mdot', 'area'),
                     ()),
    'jsr':          (('temp', 'pressure', ('res_time', 'mdot'), 'vol'),
                     ()),
    'rcm':          (('temp', 'pressure', 'end_time', 'time', 'vol'),
                     ()),
    'const_t_p':    (('temp', 'pressure', 'end_time'),
                     ())
}

# Allowed measurement types (keys) and:
# (1) required inputs in the 'plot' field of the 'info' sheet,
# (2) required inputs in the 'conds' field of the exp_sheet, and
# (3) optional inputs.
ALLOWED_MEAS_TYPES = {
    'abs':          (('variable', 'timestep', 'wavelength'),
                     ('abs_coeff', 'path_length'),
                     ()),
    'emis':         (('variable', 'timestep', 'wavelength'),
                     (),
                     ()),
    'idt':          (('variable', 'start', 'end', 'inc', 'idt_target',
                      'idt_method', 'end_time'),
                     (),
                     ()),
    'outlet':       (('variable', 'start', 'end', 'inc'),
                     (),
                     ()),
    'ion':          (('variable', 'timestep',),
                     (),
                     ()),
    'pressure':     (('variable', 'timestep',),
                     (),
                     ()),
    'conc':         (('variable', 'timestep',),
                     (),
                     ()),
}

# Allowed 'idt_method' options for determining ignition delay time
ALLOWED_IDT_METHODS = ('baseline_extrap', 'max_slope', 'max_value')

# Number of columns in the experimental sheets
NUM_EXP_CLMS = 7


def load_exp_sets(exp_filenames):
    """ Reads multiple Excel files that each contain data on an experiment set
        and returns a list of exp_set objects

        :param exp_filenames: paths to .xlsx files, each describing an exp_set
        :type exp_filenames: list
        :return exp_sets: list of exp_set objects
        :rtype: list
    """

    exp_sets = []
    for exp_filename in exp_filenames:
        exp_set = load_exp_set(exp_filename)
        exp_sets.append(exp_set)

    return exp_sets


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
    chk_exp_set(exp_set)
    for exp_obj in exp_objs:  # important to check exp_objs AFTER exp_set
        chk_exp_obj(exp_obj, exp_set)

    return exp_set


def read_info_sheet(df):
    """ Reads the Excel 'info' sheet describing information on an experiment set

        :param df: a Pandas dataframe describing the contents of the Excel sheet
        :type df: Pandas dataframe
        :return: exp_set: description of a set of experiments; used for spc list
        :rtype exp_set: dct
    """

    assert df.shape[1] >= 4, f'"info" sheets must have at least 4 columns'
    exp_set = {'overall': {}, 'plot': {}, 'mix': {}, 'spc': {},
               'exp_objs': []}
    for _, row in df.iterrows():
        assert len(row) >= 4, f'The row {list(row)} has less than 4 entries'
        group, param, raw_val, units = row[0], row[1], row[2], row[3]
        assert group in ALLOWED_INFO_GROUPS.keys(), (
            f"The group '{group}' is not allowed in the info sheet. Options "
            f"are: {tuple(ALLOWED_INFO_GROUPS.keys())}")
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
                conv_val = conv_val.split(',')
                for idx, entry in enumerate(conv_val):
                    conv_val[idx] = entry.strip()  # remove trailing/leading WS.
            # If there was no comma, convert the single value from str to lst
            # for a few special parameters
            elif param in ('idt_target', 'idt_method', 'wavelength'):
                conv_val = [conv_val]
            exp_set[group][param] = conv_val  # note: don't use fake_param
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
    inchi = rdkit_.to_inchi(rdkit_.from_smiles(smiles))
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


# Functions for checking correctness of the created objects
def chk_exp_set(exp_set):
    """ Runs checks on an exp_set to make sure certain conditions are met. This
        function operates mostly on the 'info' sheet

        :param: exp_set: description of a set of experiments
        :type exp_set: dct
    """

    set_id = exp_set['overall']['set_id']  # for printing
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']

    # Check that the groups have the required input params
    for group, params in exp_set.items():
        for reqd_inp in ALLOWED_INFO_GROUPS[group]:
            assert reqd_inp in params, (
                f"For set ID {set_id}, the param '{reqd_inp}' is missing from"
                f" the group '{group}'")

    # Check that the reactor type is allowed
    assert reac_type in ALLOWED_REAC_TYPES.keys(), (
        f"In set ID {set_id}, the reac_type '{reac_type}' is not allowed. "
        f"Options are: {tuple(ALLOWED_REAC_TYPES.keys())}")

    # Check that the measurement type is allowed
    assert meas_type in ALLOWED_MEAS_TYPES.keys(), (
        f"In set ID {set_id}, the meas_type '{meas_type}' is not allowed. "
        f"Options are: {tuple(ALLOWED_MEAS_TYPES.keys())}")

    # Check that the number of exp_objects matches the number indicated in the
    # info sheet
    num_exps = exp_set['overall']['num_exps']
    num_exp_objs = len(exp_set['exp_objs'])
    assert num_exps == num_exp_objs, (
        f"In set ID {set_id}, the 'info' sheet indicates {num_exps}"
        f" experiments, while there are {num_exp_objs} exp_objs.")

    # Check that the required inputs exist in the plot field, both according to
    # meas_type and reac_type
    plot_info = exp_set['plot']
    reqd_meas_inps = ALLOWED_MEAS_TYPES[meas_type][0]
    reqd_reac_inps = ALLOWED_REAC_TYPES[reac_type][0]
    for reqd_meas_inp in reqd_meas_inps:
        assert reqd_meas_inp in plot_info.keys(), (
            f"In the 'plot' field of set ID {set_id}, the required input"
            f" '{reqd_meas_inp}' is missing. This is required for the meas_type"
            f" '{meas_type}'.")
    # If the meas_type requires the full 'plot' inputs (as indicated by the
    # 'inc' field), check for the proper inputs according to the reactor type
    if 'inc' in reqd_meas_inps:
        plot_var = plot_info['variable']
        for reqd_reac_inp in reqd_reac_inps:
            # If reqd_reac_inp is a string, i.e., a single value
            if isinstance(reqd_reac_inp, str):
                assert reqd_reac_inp in plot_info.keys() or reqd_reac_inp \
                    in plot_var, (
                    f"In the 'plot' field of set ID {set_id}, the required"
                    f" input '{reqd_reac_inp}' is missing. This is required for"
                    f" the reac_type '{reac_type}'.")
            # If reqd_reac_inp is a tuple with multiple optional variables
            else:
                satisfied = False
                # Make sure  one and only one of the optional inputs was given
                for opt_inp in reqd_reac_inp:
                    if opt_inp in plot_info.keys():
                        assert not satisfied, (
                            f"For set ID {set_id}, only one of the optional"
                            f" inputs, {reqd_reac_inp}, should be given.")
                        satisfied = True
                assert satisfied, (
                    f"For set ID {set_id}, one of the optional inputs,"
                    f" {reqd_reac_inp}, should be given.")

    # Check that only allowed plot instructions are given
    poss_inps = get_poss_inps(reac_type, meas_type, 'plot', remove_bad=False)
    for inp in plot_info.keys():
        assert inp in poss_inps, (
                f"For set ID {set_id}, the input '{inp}' is not permitted for"
                f" the reactor type '{reac_type}' and the measurement type"
                f" '{meas_type}. Options are {poss_inps}.")

    # Check that the number of wavelengths indicated in the 'plot' field matches
    # that given in the results of the exp_sheets
    if meas_type in ('abs', 'emis'):
        num_exp_wvlens = get_num_wvlens(exp_set)
        num_set_wvlens = len(exp_set['plot']['wavelength'])
        assert num_exp_wvlens == num_set_wvlens, (
            f"For set ID {set_id}, there are {num_exp_wvlens} wavelengths "
            f"indicated in the exp sheets, but {num_set_wvlens} indicated in"
            f" the 'plot' field")

    # Check some things regarding ignition delay time measurements
    if meas_type == 'idt':
        idt_targets = exp_set['plot']['idt_target']
        idt_methods = exp_set['plot']['idt_method']
        num_targets = len(idt_targets)
        num_methods = len(idt_methods)
        # Check that the number of targets and methods are the same
        assert num_targets == num_methods, (
            f"For set ID {set_id}, the number of IDT targets, {idt_targets}, "
            f"and IDT methods, {idt_methods}, are different.")
        # Check that the methods are allowed
        for idt_method in idt_methods:
            assert idt_method in ALLOWED_IDT_METHODS, (
                f"For set ID {set_id}, the idt_method '{idt_method}' is not "
                f"allowed. Options are: {ALLOWED_IDT_METHODS}")
        # Check that the targets are allowed
        all_spcs = tuple(exp_set['spc'].keys())
        for idt_target in idt_targets:
            if idt_target != 'pressure':
                assert idt_target in all_spcs, (
                    f"For set ID {set_id}, the idt_target '{idt_target}' is not"
                    f" allowed. Options are either 'pressure' or any of the "
                    f"defined species: {all_spcs}")
        # Check that the number of IDTs indicated in the 'plot' field matches
        # that given in the results of the exp_sheets
        num_exp_idts = get_num_idts(exp_set)
        num_set_idts = len(exp_set['plot']['idt_target'])
        assert num_exp_idts == num_set_idts, (
            f"For set ID {set_id}, there are {num_exp_idts} IDTs "
            f"indicated in the exp sheets, but {num_set_idts} indicated in"
            f" the 'plot' field")

    # Check that the end_time and timestep are compatible
    timestep = exp_set['plot'].get('timestep')
    end_time = exp_set['plot'].get('end_time')
    if timestep is not None and end_time is not None:
        assert np.isclose(end_time % timestep, 0, atol=1e-8), (
            f'For set ID {set_id}, the end_time, {end_time}, is not a multiple '
            f'of the timestep, {timestep}.')


def chk_exp_obj(exp_obj, exp_set):
    """ Runs checks on an exp_obj to make sure certain conditions are met

        :param exp_obj: description of a single experiment
        :type: dct
        :param: exp_set: description of a set of experiments
        :type exp_set: dct
    """

    set_id = exp_set['overall']['set_id']  # for printing

    # Check that the experiment ID is defined
    assert exp_obj['overall'].get('exp_id') is not None, (
        f"For set ID {set_id}, one exp_sheet is missing the 'exp_id' field")
    exp_id = exp_obj['overall']['exp_id']  # for printing

    # Check for all required inputs for the given reactor type
    reac_type = exp_set['overall']['reac_type']
    reqd_reac_inps = ALLOWED_REAC_TYPES[reac_type][0]
    for reqd_reac_inp in reqd_reac_inps:
        # If the required input is a single required value (e.g., 'temp')
        if isinstance(reqd_reac_inp, str):
            assert reqd_reac_inp in exp_obj['conds'].keys(), (
                f"For set ID {set_id}, exp ID {exp_id}, the required input"
                f" '{reqd_reac_inp}' is absent from the 'conds' field.")
        # If the required input is a tuple, only one of which is required (e.g.,
        # 'res_time' or 'mdot' for a JSR)
        else:  # in this case, reqd_reac_inp is a tuple with multiple variables
            satisfied = False
            # Make sure that one and only one of the optional inputs was given
            for opt_inp in reqd_reac_inp:
                if opt_inp in exp_obj['conds'].keys():
                    assert not satisfied, (
                        f"For set ID {set_id}, exp ID {exp_id}, only one of the"
                        f" optional inputs, {reqd_reac_inp}, should be given.")
                    satisfied = True
            assert satisfied, (
                f"For set ID {set_id}, exp ID {exp_id}, one of the"
                f" optional inputs, {reqd_reac_inp}, should be given.")

    # Check for all required inputs for the given measurement type
    meas_type = exp_set['overall']['meas_type']
    reqd_meas_inps = ALLOWED_MEAS_TYPES[meas_type][1]
    for reqd_meas_inp in reqd_meas_inps:
        assert reqd_meas_inp in exp_obj['conds'].keys(), (
            f"For set ID {set_id}, exp ID {exp_id}, the required input"
            f" {reqd_meas_inp} is absent from the 'conds' field.")

    # Check that only allowed conditions are given
    poss_inps = get_poss_inps(reac_type, meas_type, 'exps')
    for inp in exp_obj['conds'].keys():
        assert inp in poss_inps, (
            f"For set ID {set_id}, exp ID {exp_id}, the input '{inp}' is not"
            f" permitted for the reactor type '{reac_type}' and the measurement"
            f" type '{meas_type}'. Options are {poss_inps}.")

    # Check that all arrays in the result section are the same length
    lengths = []
    for name, result in exp_obj['result'].items():
        if isinstance(result, tuple):  # if a tuple, result is a val_tuple
            val = result[0]  # the array is the first entry
            if isinstance(val, np.ndarray):
                val = val[~np.isnan(val)]  # remove NaNs
                lengths.append(len(val))
        elif isinstance(result, dict):  # if dct, result is a dct of val_tuples
            for val_tuple in result.values():
                val = val_tuple[0]  # the array is the first entry
                if isinstance(val, np.ndarray):
                    val = val[~np.isnan(val)]  # remove NaNs
                    lengths.append(len(val))

    # The lengths should all be the same (i.e., 1) or empty if no arrays (0)
    assert len(set(lengths)) in (0, 1), (
        f"For set ID {set_id}, exp ID {exp_id}, the time-resolved arrays in the"
        f" 'result' section are not all the same length")

    # Check some things regarding absorption measurements
    if meas_type == 'abs':
        # Check that the keys of the result section are formatted correctly
        for name, result in exp_obj['result'].items():
            if name == 'abs':
                assert isinstance(result, dict), (
                    f"For set ID {set_id}, exp ID {exp_id}, the 'abs' result"
                    f" should be a dict.")
                for key in result.keys():
                    assert isinstance(key, int), (
                        f"For set ID {set_id}, exp ID {exp_id}, the 'abs'"
                        f" result should have integers as the keys (1, 2, ...)")

    # Check that the end_time and timestep are compatible
    timestep = exp_set['plot'].get('timestep')
    end_time = exp_obj['conds'].get('end_time')[0]
    if timestep is not None and end_time is not None:
        assert np.isclose(end_time % timestep, 0, atol=1e-8), (
            f'For set ID {set_id}, exp ID {exp_id}, the end_time, {end_time},'
            f' is not a multiple of the timestep, {timestep}.')


# Miscellaneous functions
def get_poss_inps(reac_type, meas_type, plot_or_exps, remove_bad=True):
    """ Gets all possible (i.e., allowed) inputs for either the 'plot' field of
        the 'info' sheet or the 'conds' field of the exp sheets

        :param reac_type: the reactor type
        :type reac_type: str
        :param meas_type: the measurement type
        :type meas_type: str
        :param plot_or_exps: either 'plot' or 'exps'
        :type plot_or_exps: str
        :param remove_bad: whether or not to remove the plot-specific inputs
        :type remove_bad: Bool
        :return poss_inps: all possible inputs
        :rtype: tuple (str1, str2, ...)
    """

    # Get the required inputs and flatten them
    if plot_or_exps == 'plot':
        reqd_meas_inps = ALLOWED_MEAS_TYPES[meas_type][0]
    elif plot_or_exps == 'exps':
        reqd_meas_inps = ALLOWED_MEAS_TYPES[meas_type][1]
    reqd_reac_inps = ALLOWED_REAC_TYPES[reac_type][0]
    flat_reqd_meas_inps = flatten_tup(reqd_meas_inps)
    flat_reqd_reac_inps = flatten_tup(reqd_reac_inps)

    # Get the optional inputs
    opt_meas_inps = ALLOWED_MEAS_TYPES[meas_type][2]
    opt_reac_inps = ALLOWED_REAC_TYPES[reac_type][1]

    # Add together
    poss_inps = flat_reqd_meas_inps + flat_reqd_reac_inps + opt_meas_inps +\
        opt_reac_inps

    # Remove the plotting variable stuff if indicated
    if remove_bad:
        bad_items = ('variable', 'start', 'end', 'inc')
        poss_inps = list(poss_inps)
        for bad_item in bad_items:
            if bad_item in poss_inps:
                poss_inps.remove(bad_item)
        poss_inps = tuple(poss_inps)

    return poss_inps


def flatten_tup(inp_tup):
    """ If any elements in a tuple are themselves tuples, flattens out the tuple
        into a single, flat tuple

        Example: if inp_tup = (1, (2, 3), 4), the output flat
        tuple will be flat_tup = (1, 2, 3, 4)
    """

    flat_tup = []
    for entry in inp_tup:
        if isinstance(entry, str):
            flat_tup.append(entry)
        else:  # if a tuple
            for sub_entry in entry:
                flat_tup.append(sub_entry)
    flat_tup = tuple(flat_tup)

    return flat_tup


def get_num_wvlens(exp_set):
    """ Gets the number of wavelengths specified throughout the experimental set
    """

    meas_type = exp_set['overall']['meas_type']
    num_wvlens = 0
    for exp_obj in exp_set['exp_objs']:
        current_num = len(exp_obj['result'][meas_type])
        if current_num > num_wvlens:
            num_wvlens = current_num

    return num_wvlens


def get_num_idts(exp_set):
    """ Gets the number of IDT values specified throughout the experimental set
    """

    num_idts = 0
    for exp_obj in exp_set['exp_objs']:
        current_num = len(exp_obj['result']['idt'])
        print('current_num:', current_num)
        if current_num > num_idts:
            num_idts = current_num

    return num_idts
