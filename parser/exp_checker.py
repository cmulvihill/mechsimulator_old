import numpy as np


# Allowed groups (keys) and their required params (values) in the 'info' sheet
ALLOWED_INFO_GROUPS = {
    'overall':      ('set_id', 'source', 'description', 'reac_type',
                     'meas_type', 'num_exps', 'version'),
    'plot':         (),  # in some cases, none are required, so leaving blank
    'plot_format':  (),
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
# (3) the allowed measurement types
# In the required inputs, a tuple indicates that only one of those values is
# required (e.g., either 'res_time' or 'mdot' is required for a JSR).
ALLOWED_REAC_TYPES = {
    'st':           (('temp', 'pressure', 'end_time'),
                     ('dpdt',),
                     ('abs', 'emis', 'idt', 'outlet', 'ion', 'pressure',
                      'conc')),
    'pfr':          (('temp', 'pressure', 'length', 'mdot', 'area'),
                     (),
                     ('outlet',)),
    'jsr':          (('temp', 'pressure', ('res_time', 'mdot'), 'vol'),
                     (),
                     ('outlet',)),
    'rcm':          (('temp', 'pressure', 'end_time', 'time', 'vol'),
                     (),
                     ('idt', 'pressure')),
    'const_t_p':    (('temp', 'pressure', 'end_time'),
                     (),
                     ('abs', 'emis', 'idt', 'outlet', 'ion', 'pressure',
                      'conc'))
}

# Allowed measurement types (keys) and:
# (1) required inputs in the 'plot' field of the 'info' sheet,
# (2) required inputs in the 'conds' field of the exp_sheet, and
# (3) optional inputs.
ALLOWED_MEAS_TYPES = {
    'abs':          (('variable', 'timestep', 'end_time', 'wavelength'),
                     ('abs_coeff', 'path_length'),
                     ('start', 'end', 'inc',)),
    'emis':         (('variable', 'timestep', 'end_time', 'wavelength'),
                     (),
                     ('start', 'end', 'inc',)),
    'idt':          (('variable', 'start', 'end', 'inc', 'idt_target',
                      'idt_method', 'end_time'),
                     (),
                     ()),
    'outlet':       (('variable', 'start', 'end', 'inc'),
                     (),
                     ()),
    'ion':          (('variable', 'timestep', 'end_time',),
                     (),
                     ('start', 'end', 'inc',)),
    'pressure':     (('variable', 'timestep', 'end_time',),
                     (),
                     ('start', 'end', 'inc',)),
    'conc':         (('variable', 'timestep', 'end_time',),
                     (),
                     ('start', 'end', 'inc',)),
}

# Miscellaneous plotting inputs that are allowed (although will be often unused)
MISC_PLOT_INPUTS = ('group_by', 'rows_cols', 'marker_size',
                    'exp_color', 'plot_points', 'xunit', 'xlim', 'yunit',
                    'ylim', 'omit_targets')

# Allowed 'idt_method' options for determining ignition delay time
ALLOWED_IDT_METHODS = ('baseline_extrap', 'max_slope', 'max_value')


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
    assert meas_type in ALLOWED_REAC_TYPES[reac_type][2], (
        f"In set ID {set_id}, meas_type '{meas_type}' is not allowed for "
        f"reac_type {reac_type}. Options are: "
        f"{ALLOWED_REAC_TYPES[reac_type][2]}")

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
    poss_inps = get_poss_inps(reac_type, meas_type, 'plot', rm_bad=False)
    for inp in plot_info.keys():
        assert inp in poss_inps, (
                f"For set ID {set_id}, the input '{inp}' is not permitted for"
                f" the reactor type '{reac_type}' and the measurement type"
                f" '{meas_type}'. Options are {poss_inps}.")

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

    # Check some things regarding the plot_options field
    plot_format = exp_set['plot_format']
    # Check that only allowed groups were given
    for param, val in plot_format.items():
        assert param in MISC_PLOT_INPUTS, (
            f"For set ID {set_id}, the plot_format input '{param}' is not "
            f"allowed. Options are {MISC_PLOT_INPUTS}.")
        if param == 'group_by':
            assert val in ('target', 'cond'), (
                f"'group_by' should be 'target' or 'cond', not {val}")
        if param == 'plot_points':
            assert isinstance(val, bool), (
                f"'plot_points' should be either 'yes' or 'no', not {val}")
        # Checks on the xunit and yunit are done in the plotting code
        if param in ('xunit', 'yunit'):
            assert val not in ('C', 'F'), (
                f"For set ID {set_id}, the x/y unit cannot be {val}. "
                f"Temperature must be in K.")


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


# Miscellaneous functions
def get_poss_inps(reac_type, meas_type, plot_or_exps, rm_bad=True):
    """ Gets all possible (i.e., allowed) inputs for either the 'plot' field of
        the 'info' sheet or the 'conds' field of the exp sheets

        :param reac_type: the reactor type
        :type reac_type: str
        :param meas_type: the measurement type
        :type meas_type: str
        :param plot_or_exps: either 'plot' or 'exps'
        :type plot_or_exps: str
        :param rm_bad: whether or not to remove the plot-specific inputs
        :type rm_bad: Bool
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
    if rm_bad:
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
        if current_num > num_idts:
            num_idts = current_num

    return num_idts
