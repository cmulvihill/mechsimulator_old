import numpy as np
import matplotlib.pyplot as plt
import cantera as ct
from mechsimulator import parser
from math import floor
from math import ceil
from mechsimulator.parser.exp import ALLOWED_UNITS


MEAS_TYPE_DIMS = {
    'ion':          4,
    'abs':          3,
    'emis':         3,
    'pressure':     3,
    'temp':         3,
    'burner_conc':  3,
    'burner_temp':  3,
    'conc':         3,
    'idt':          2,
    'outlet':       2,
    'lfs':          2,
}
TIME_RESOLVED_MEAS_TYPES = (
    'abs',
    'emis',
    'ion',
    'pressure',
    'conc',
)


def get_target_spcs(exp_set, rename_instr):
    """ Gets the list of target species from an exp_set (i.e., all species
        defined in the 'spc' block of the 'info' sheet) and translates them to
        the names used in the mechanism

        :param exp_set: object describing a set of experiments
        :type exp_set: dct
        :param rename_instr: spc_names in the experiment to be renamed according
            to the mechanism spc names
        :type rename_instr: dct {spc_to_be_renamed: new_spc_name, ...}
        :return target_spcs: all species defined in the 'spc' block, translated
            to the names used in the mechanism
        :rtype: list
    """

    spc_smiles_dct = exp_set['spc']
    raw_target_spcs = list(spc_smiles_dct.keys())
    target_spcs = translate_spcs(raw_target_spcs, rename_instr)

    return target_spcs


def get_rename_instr(mech_spc_dct1, mech_spc_dct2):
    """ Gets the rename_instr, which relates the species names in the Excel
        file with those in the mechanism. Does this by comparing the two
        mech_spc_dcts

        :param mech_spc_dct1: the reference mech_spc_dct
        :param mech_spc_dct2: the mech_spc_dct to be renamed
        :return rename_instr: spc_names in mech_spc_dct2 to be renamed
        :rtype: dct {spc_to_be_renamed: new_spc_name, ...}
    """

    rename_instr = {}
    # Loop over each species in the mech_spc_dct to be renamed
    for spc2, spc_dct2 in mech_spc_dct2.items():
        ich2 = spc_dct2['inchi']
        mlt2 = spc_dct2['mult']
        chg2 = spc_dct2['charge']
        exc2 = spc_dct2['exc_flag']
        # Loop over each species in the reference mech_spc_dct
        for spc1, spc_dct1 in mech_spc_dct1.items():
            ich1 = spc_dct1['inchi']
            mlt1 = spc_dct1['mult']
            chg1 = spc_dct1['charge']
            exc1 = spc_dct1['exc_flag']
            if ich2 == ich1 and mlt2 == mlt1 and chg2 == chg1 and exc1 == exc2:
                rename_instr[spc2] = spc1
                break  # no need to continue
        # If spc2 was not found in mech_spc_dct1, save value as None
        if spc2 not in rename_instr.keys():
            rename_instr[spc2] = None

    return rename_instr


def translate_spcs(names, rename_instr):
    """ Renames the species inside of an object according to the instructions
        inside rename_instr. The names input can either be a list or a dct

        :param names: an input with spc names that are to be redefined
        :type names: dct {spc1: val1, spc2: ...}, list [spc1, spc2, ...], or
            tuple (spc1, spc2, ...)
        :param rename_instr: spc_names in the experiment to be renamed according
            to the mechanism spc names
        :type rename_instr: dct {spc_to_be_renamed: new_spc_name, ...}
        :return translated_names: the input names but with all spcs renamed
        :rtype: dct, list, or tuple (depending on the input)
    """

    assert isinstance(names, (dict, list, tuple)), (
        f'Input "names" must be a dict, list, or tuple but is a {type(names)}')

    # Rename a dictionary with the format {spc1: value1, spc2: ...}; this will
    # usually be a dictionary describing the initial mixture
    if isinstance(names, dict):
        translated_names = {}
        for spc, value in names.items():
            renamed_spc = rename_instr[spc]
            assert renamed_spc is not None, (
                f'The species {spc} appears in the initial mixture but does not'
                ' appear in the mechanism. The simulation cannot proceed.')
            translated_names[renamed_spc] = value

    # Rename a list or tuple of spcs: [spc1, spc2, ...] or (spc1, spc2, ...);
    # this will often be a list of target_spcs
    else:  # if a list or tuple
        translated_names = []
        for spc in names:
            renamed_spc = rename_instr[spc]
            translated_names.append(renamed_spc)
        if isinstance(names, tuple):  # change to tuple if needed
            translated_names = tuple(translated_names)

    return translated_names


def check_time_resolved(exp_set):
    """ Checks if an exp_set has a time-resolved measurement

        :param exp_set:
        :type exp_set:
        :return is_time_resolved: whether or not an exp_set is time resolved
        :rtype: Bool
    """

    meas_type = exp_set['overall']['meas_type']
    if meas_type in TIME_RESOLVED_MEAS_TYPES:
        is_time_resolved = True
    else:
        is_time_resolved = False

    return is_time_resolved


def load_solution_objs(mech_filenames):
    """ Receives a list of Cantera .cti mechanism filenames and returns a list
        of Cantera Solution objects

        :param mech_filenames: list of .cti filenames
        :type
        :return mechs: list of Cantera Solution objects
        :rtype: list [solution1, solution2, ...]
    """

    mechs = []
    for mech_filename in mech_filenames:
        mech = ct.Solution(mech_filename)
        mech.name = mech_filename  # filename is associated with the mech
        mechs.append(mech)

    return mechs


def get_plot_variable(exp_set):
    """ Gets the identity of the plotting variable
    """
    plot_dct = exp_set['plot']
    plot_var = plot_dct['variable']

    return plot_var


def get_plot_conds(exp_set):
    """ Gets a Numpy array describing the variation of the plotting variable
    """
    plot_dct = exp_set['plot']
    start = plot_dct['start']
    end = plot_dct['end']
    inc = plot_dct['inc']
    num = int((end - start) / inc) + 1
    plot_conds = np.linspace(start, end, num)

    return plot_conds


def get_uniform_times(exp_set):
    """ Gets a uniform times array for an exp_set based on the plot instructions

    """

    max_time = exp_set['plot']['end_time']
    min_time = 0
    timestep = exp_set['plot']['timestep']
    num_times = int((max_time - min_time) / timestep) + 1
    times = np.linspace(min_time, max_time, num_times)

    return times


def interp(ydata, xdata, desired_xdata):
    """ Takes a target array with some inconsistent and/or unknown step size in
        the x variable and interpolates the y data to make the x data have a
        uniform step size. The input ydata may be a Numpy array of dimension 1
        or 2 (the second dimension must be the x dimension).

        Usually, x will be time, but it does not have to be

        :param ydata: Numpy array containing y data to be interpolated; can be
            one array (ndim=1) or multiple arrays (ndim=2) to be interpolated
        :type ydata: Numpy.ndarray
        :param xdata: Numpy array containing the x data corresponding to the y
            data; must have ndim=1 and be monotonically increasing
        :type xdata: Numpy.ndarray
        :param desired_xdata: the desired grid of x values to which the y data
            will be interpolated; must have ndim=1 and be monotonically increasing
        :return interp_ydata: the interpolated y data
        :rtype: Numpy.ndarray
    """

    num_dims = np.ndim(ydata)
    assert num_dims in (1, 2), (
        f'correct_times requires arrays of dimension 1 or 2, not {num_dims}')
    assert min(np.diff(xdata)) > 0, 'xdata should be monotonically increasing'

    # Interpolate the data; method depends slightly on the dimensionality
    if num_dims == 1:
        interp_ydata = np.interp(desired_xdata, xdata, ydata)
    else:  # num_dims = 2
        num_arrs = np.shape(ydata)[0]  # assumes time is the second dim
        num_times = len(desired_xdata)
        interp_ydata = np.ndarray((num_arrs, num_times))
        for arr_idx, arr in enumerate(ydata):
            interp_ydata[arr_idx, :] = np.interp(desired_xdata, xdata, arr)

    # Get idxs of the high and low cutoffs beyond which the data doesn't extend
    high_cutoff = np.argmin(abs(desired_xdata - xdata[-1]))  # first match
    low_cutoff = np.argmin(abs(desired_xdata - xdata[0]))  # first match
    # Get rid of interpolated y data that extend past the data
    if num_dims == 1:
        interp_ydata[high_cutoff:] = np.nan
        interp_ydata[:low_cutoff] = np.nan
    else:  # num_dims = 2
        interp_ydata[:, high_cutoff:] = np.nan
        interp_ydata[:, :low_cutoff] = np.nan

    return interp_ydata


def get_mech_info(exp_set, calc_type, x_source, conds_source, gases):
    """ Gets the shape of the Numpy array that will describe a single mechanism
        simulation of an exp_set. Works for three different calculation types:
        (1) outcome, (2) sensitivity, or (3) rate of production

    """

    def get_times(exp_set, x_source):

        if x_source == 'plot':
            num_times = get_uniform_times(exp_set)
        else:
            num_times = exp_set['overall']['xdata']

        return num_times

    def get_cond_titles(conds, units):

        cond_titles = []
        for cond in conds:
            cond_title = f'{cond} {units}'
            cond_titles.append(cond_title)

        return cond_titles

    def get_idt_targets_titles(exp_set):
        fake_targets = exp_set['plot']['idt_target']
        methods = exp_set['plot']['idt_method']
        targets = []
        titles = []
        for idx, fake_target in enumerate(fake_targets):
            targets.append((fake_target, methods[idx]))
            titles.append(fake_target + ',' + methods[idx])

        return targets, titles

    # Get the set variable and check that options were used correctly
    plot_var = get_plot_variable(exp_set)
    units = ALLOWED_UNITS[plot_var][0][0]
    if x_source == 'exps' and conds_source == 'plot':
        raise NotImplementedError(
            "x_source='exps' and conds_source='plot' is not allowed")

    # Get the conditions depending on the source
    if conds_source == 'plot':
        conds = get_plot_conds(exp_set)
        cond_titles = get_cond_titles(conds, units)
    else:  # 'exps'
        conds = np.arange(1, exp_set['overall']['num_exps'] + 1, 1)  # 1 to N
        cond_titles = [str(cond) for cond in conds]  # maybe totally useless

    # Get the targets (and possibly another variable)
    # Note: targets doesn't depend on the x_source
    meas_type = exp_set['overall']['meas_type']
    if meas_type in ('abs', 'emis'):
        targets = exp_set['plot']['wavelength']
        times = get_times(exp_set, x_source)
        shape = (len(conds), len(targets), len(times))
        x_arrays = (conds, targets, times)
        x_titles = (cond_titles, targets, None)
        xdata = times
    elif meas_type == 'conc':
        targets = list(exp_set['spc'].keys())
        times = get_times(exp_set, x_source)
        shape = (len(conds), len(targets), len(times))
        x_arrays = (conds, targets, times)
        x_titles = (cond_titles, targets, None)
        xdata = times
    elif meas_type == 'pressure':
        targets = ['Pressure']
        times = get_times(exp_set, x_source)
        shape = (len(conds), len(targets), len(times))
        x_arrays = (conds, targets, times)
        x_titles = (cond_titles, targets, None)
        xdata = times
    elif meas_type == 'idt':
        targets, target_titles = get_idt_targets_titles(exp_set)
        shape = (len(conds), len(targets))
        x_arrays = (conds, targets)
        x_titles = (cond_titles, targets)
        xdata = conds
    elif meas_type == 'outlet':
        targets = list(exp_set['spc'].keys())
        shape = (len(conds), len(targets))
        x_arrays = (conds, targets)
        x_titles = (cond_titles, targets)
        xdata = conds
    else:  # meas_type == 'ion':
        raise NotImplementedError("'ion' is not implemented yet")

    # If the calculation type is sensitivity or rate of production, add some 
    # information about the reactions
    if calc_type in ('sens', 'rop'):
        # Get the maximum number of reactions among all the mechanisms
        max_num_rxns = 0
        for gas in gases:
            if len(gas.n_reactions) > max_num_rxns:
                max_num_rxns = gas.n_reactions         
        shape = (max_num_rxns,) + shape  # prepend
        # Need to update to add reaction names to x_arrays and titles
        
    return shape, x_arrays, x_titles, xdata


def get_exp_xlims(exp_set):
    """ Gets the limits of some x variable (usually time, but could also be
        other variables) from the exp_objs in an exp_set. In the case where no
        time data is given, uses the 'end_time' entry given in the exp_sheet.

    """

    # Get the name of the x variable
    meas_type = exp_set['overall']['meas_type']
    if meas_type == 'burner':
        x_var = 'position'  # maybe should be 'xposition'?
    # Should add an option here for ion stuff
    else:
        x_var = 'time'

    # Get the limits of the x variable
    exp_objs = exp_set['exp_objs']
    timestep = exp_set['plot']['timestep']  # used to fix rounding errors
    min_x = 0
    max_x = 0
    for exp_obj in exp_objs:
        if x_var in exp_obj['result']:
            raw_min = min(exp_obj['result'][x_var][0])  # [0] omits uncertainty
            raw_max = max(exp_obj['result'][x_var][0])
            # Round raw values, in the proper direction, to nearest timestep
            current_min = floor(raw_min / timestep) * timestep  # round down
            current_max = ceil(raw_max / timestep) * timestep  # round up
            if current_min < min_x:
                min_x = current_min
            if current_max > max_x:
                max_x = current_max
        else:  # if no 'time' array, use the 'end_time' entry
            given_end_time = exp_obj['conds']['end_time'][0]
            if given_end_time > max_x:
                max_x = given_end_time

    return min_x, max_x


def plot_derivs(target, times):
    """ Toy code for visualizing derivatives

        :param target: single array of y-data
        :type target: Numpy.ndarray

    """

    # cutoff_idx = 300
    cutoff_idx = len(target)
    target = target[:cutoff_idx]
    times = times[:cutoff_idx]

    first_deriv = np.gradient(target, times)
    second_deriv = np.gradient(first_deriv, times)
    extrema_idxs = np.nonzero((np.diff(np.sign(first_deriv)) != 0) * 1)[0]
    inflection_idxs = np.nonzero((np.diff(np.sign(second_deriv)) != 0) * 1)[0]
    print('target:\n', target)
    print('diff of target:\n', np.diff(target))
    print('first_deriv:\n', first_deriv)
    print('second_deriv:\n', second_deriv)
    print('extrema_idxs:\n', extrema_idxs)
    print('length:\n', len(extrema_idxs))
    print('inflection_idxs:\n', inflection_idxs)
    print('length:\n', len(inflection_idxs))

    fig, axs = plt.subplots(3, 1)
    axs[0].plot(times, target)
    axs[1].plot(times, first_deriv)
    axs[2].plot(times, second_deriv)

    for idx in extrema_idxs:
        axs[0].plot(times[idx], target[idx], 'ro')
    # for idx in inflection_idxs:
    #     axs[0].plot(times[idx], target[idx], 'bo')

    plt.show()
