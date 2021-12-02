import matplotlib.pyplot as plt
import numpy as np
from mechsimulator.simulator import util
from mechsimulator.simulator import reactors


def get_result(conds_dct, gas, reac_type, meas_type, xdata, ydata_shape,
               prev_soln=None):

    if reac_type == 'jsr':
        mech_ydata, _ = jsr(conds_dct, gas, meas_type, ydata_shape,
                            prev_soln=prev_soln)
    elif reac_type == 'st':
        mech_ydata = st(conds_dct, gas, meas_type, xdata, ydata_shape)
    else:
        raise NotImplementedError

    return mech_ydata


def jsr(conds_dct, gas, meas_type, ydata_shape, prev_soln=None,
        output_all=False):

    # Get arrays of reactor inputs
    temps = conds_dct['temp']
    pressures = conds_dct['pressure']
    mixes = conds_dct['mix']
    res_times = conds_dct['res_time']
    vols = conds_dct['vol']
    mdots = conds_dct['mdot']  # either res_time or mdot will be used
    target_spcs = conds_dct['target_spcs']

    # Make sure that the prev_soln array is the right size (if given)
    if prev_soln is not None:
        prev_soln_shape = np.shape(prev_soln)
        assert prev_soln_shape[0] == ydata_shape[0]  # [0] is # of conditions
        assert prev_soln_shape[1] == gas.n_species

    # If output_all was indicated, initialize
    if output_all:  # all_concs has shape = (num_conds, num_species)
        all_concs = np.ndarray((ydata_shape[0], gas.n_species))
    else:
        all_concs = None

    # Loop over the conditions, performing a simulation for each
    mech_ydata = np.ndarray(ydata_shape)
    for cond_idx in range(ydata_shape[0]):
        temp = temps[cond_idx]
        pressure = pressures[cond_idx]
        mix = mixes[cond_idx]
        res_time = res_times[cond_idx]
        vol = vols[cond_idx]
        mdot = mdots[cond_idx]
        # Get the previous concentrations
        if prev_soln is not None:
            prev_concs = prev_soln[cond_idx]
        elif cond_idx == 0:        # if prev_soln not given and on first time,
            prev_concs = mixes[0]  # use mix as prev_concs
        # Run the simulation
        raw_concs, prev_concs, _, _ = reactors.jsr(
            temp, pressure, mix, gas, target_spcs, res_time, vol, mdot=mdot,
            prev_concs=prev_concs)
        # If all concentrations were requested, store them
        if output_all:
            all_concs[cond_idx, :] = prev_concs
        # "Process" the raw result
        mech_ydata[cond_idx, :] = raw_concs

    return mech_ydata, all_concs
    

def st(conds_dct, gas, meas_type, xdata, ydata_shape):
    """ Gets the raw results from a shock tube simulation

        :param conds_dct:
        :type conds_dct:
        :param gas:
        :type gas:
        :return:
    """

    # Get arrays of reactor inputs
    temps = conds_dct['temp']
    pressures = conds_dct['pressure']
    mixes = conds_dct['mix']
    end_times = conds_dct['end_time']
    target_spcs = conds_dct['target_spcs']
    p_of_ts = conds_dct['p_of_t']

    # Loop over all conditions
    mech_ydata = np.ndarray(ydata_shape)
    for cond_idx in range(ydata_shape[0]):  # [0] gives # of conditions
        temp = temps[cond_idx]
        pressure = pressures[cond_idx]
        mix = mixes[cond_idx]
        end_time = end_times[cond_idx]
        p_of_t = p_of_ts[cond_idx]
        # Perform a simulation
        raw_concs, raw_pressures, raw_temps, raw_times, _, _ = reactors.st(
            temp, pressure, mix, gas, target_spcs, end_time, p_of_t=p_of_t)
        # Process the raw results
        mech_ydata[cond_idx] = process_st(
            raw_concs, raw_pressures, raw_temps, raw_times, conds_dct, cond_idx,
            meas_type, xdata)

    return mech_ydata


def process_st(raw_concs, raw_pressures, raw_temps, raw_times, conds_dct,
               cond_idx, meas_type, uniform_times):
    """ Processes the outcome for a single condition of a shock tube simulation

        :param raw_concs:
        :param raw_pressures:
        :param raw_temps:
        :param raw_times:
        :param conds_dct:
        :param cond_idx: index of the current condition; used for things that
            are cond. dependent (e.g., absorption coefficient)
        :type cond_idx:
        :return:
    """

    if meas_type == 'conc':
        # Simply interpolate the raw concentrations to fit the uniform times
        cond_ydata = util.interp(raw_concs, raw_times, uniform_times)

    elif meas_type == 'idt':
        # Load information from conds_dct
        idt_targets = conds_dct['idt_target']
        idt_methods = conds_dct['idt_method']
        target_spcs = conds_dct['target_spcs']
        # Loop over each target/method pair and get the IDT data
        cond_ydata = np.ndarray(len(idt_targets))
        for target_idx, idt_target in enumerate(idt_targets):
            idt_method = idt_methods[target_idx]
            if idt_target == 'pressure':
                target_profile = raw_pressures
            else:  # if a species
                spc_idx = target_spcs.index(idt_target)  # idx of spc in targets
                target_profile = raw_concs[spc_idx]  # conc. of target spc
            idt, warnings = st_idt(target_profile, raw_times, method=idt_method)
            if warnings:
                print(f'Warning: {warnings} was returned for ST IDT simulation')
            cond_ydata[target_idx] = idt
    else:
        raise NotImplementedError(f"meas_type '{meas_type}' not working for st")

    return cond_ydata


# Functions for getting IDT
def _get_rcm_idt(pressures, times):
    """ Extracts ignition delay time from a rapid compression machine pressure
        time history

        :param pressures: the pressure time history used to determine IDT
        :type pressures: Numpy array of shape (num_timesteps,)
        :param times: the time values corresponding to target
        :type times: Numpy array of shape (num_timesteps,)
        :return idt: ignition delay time (s)
        :rtype: float
        :return warnings: possible warnings regarding the IDT determination
        :rtype: list of strs
    """

    # Get the extrema idxs, which are where the first_deriv changes signs
    first_deriv = np.gradient(pressures, times)
    second_deriv = np.gradient(first_deriv, times)
    extrema_idxs = np.nonzero((np.diff(np.sign(first_deriv)) != 0) * 1)[0]

    # Get the maxima; i.e., where second_deriv < 0
    maxima_idxs = []
    for extrema_idx in extrema_idxs:
        if second_deriv[extrema_idx] < 0:
            maxima_idxs.append(extrema_idx)

    # If there are two peaks, calculate IDT; otherwise, return a warning
    idt = None
    warning = None
    if len(maxima_idxs) == 2:
        time_zero = times[maxima_idxs[0]]
        ignition = times[maxima_idxs[1]]
        idt = ignition - time_zero
    else:
        warning = f'Instead of 2 peaks, {len(maxima_idxs)} peaks were found'

    return idt, warning


def st_idt(target, times, method='baseline_extrap', plot=False):
    """ Gets the ignition delay time for a shock-tube simulation. Can determine
        IDT using one of three methods.

        :param target: the target concentration used to determine IDT
        :type target: Numpy array of shape (num_timesteps,)
        :param times: the time values corresponding to target
        :type times: Numpy array of shape (num_timesteps,)
        :param method: the method by which to determine the IDT; options are as
            follows:
            1: intersection of steepest slope with baseline
            2: point of steepest slope
            3: peak value of the target profile
        :type method: str
        :return idt: ignition delay time (s)
        :rtype: float
        :return warnings: possible warnings regarding the IDT determination
        :rtype: list of strs
    """

    warnings = []

    # Get first derivative (note: np.gradient uses central differences)
    first_deriv = np.gradient(target, times)
    steepest_idx = np.argmax(first_deriv)
    steepest_slope = first_deriv[steepest_idx]
    steepest_time = times[steepest_idx]
    steepest_val = target[steepest_idx]

    # If using the baseline extrapolation or steepest slope methods, check that
    # the max slope isn't the last point
    if steepest_idx + 1 == len(times) and method in ('baseline_extrap',
                                                     'max_slope'):
        warnings.append('Max slope at last point')

    if method == 'baseline_extrap':
        # Get the slope and intercept of the baseline
        initial_slope = 0  # for now, assuming these are 0; could be better
        initial_int = 0

        # Get the y-intercept of the steepest tangent line
        steepest_int = steepest_val - steepest_slope * steepest_time

        # Find the intersection of the two lines
        idt = (initial_int - steepest_int) / (steepest_slope - initial_slope)
        int_y = steepest_slope * idt + steepest_int  # for plotting
    elif method == 'max_slope':
        idt = steepest_time
    elif method == 'max_value':
        # Check that the max value doesn't occur at the last point
        if np.argmax(target) + 1 == len(times):
            warnings.append('Peak value at last point')
        idt = times[np.argmax(target)]
    else:
        raise NotImplementedError

    # Perform plotting if indicated
    if plot:
        plt.plot(times, target)
        if method == 'baseline_extrap':
            xlim = plt.xlim()  # store current xlim and ylim
            ylim = plt.ylim()
            tangent = steepest_slope * times + steepest_int
            plt.plot(times, tangent, 'r--')
            plt.plot(idt, int_y, 'ro')
            plt.xlim(xlim)  # reset xlim and ylim to zoom in on the target data
            plt.ylim(ylim)
        if method == 'max_slope':
            steepest_value = target[np.argmax(first_deriv)]
            plt.plot(idt, steepest_value, 'ro')
        if method == 'max_value':
            peak_value = max(target)
            plt.plot(idt, peak_value, 'ro')
        plt.show()

    if warnings:
        idt = np.nan

    return idt, warnings
