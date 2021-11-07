import copy
import matplotlib.pyplot as plt
import numpy as np
from mechsimulator.sim import util
from mechsimulator.sim import reactors
from mechsimulator.parser.exp import ALLOWED_IDT_METHODS


def get_result(conds_dct, gas, reac_type, meas_type):

    raw_spcs, raw_pressures, raw_temps, uniform_times = get_raw_result(
        conds_dct, gas, reac_type)
    mech_result, uniform_times = process_raw_result(
        raw_spcs, raw_pressures, raw_temps, uniform_times, meas_type, conds_dct)

    return mech_result, uniform_times


def get_raw_result(conds_dct, gas, reac_type):

    if reac_type == 'jsr':
        raw_mech_result = jsr(conds_dct, gas)
        uniform_times = None
    elif reac_type == 'st':
        raw_spcs, raw_pressures, raw_temps, uniform_times = st(conds_dct, gas)
    else:
        raise NotImplementedError
    # elif reac_type == 'pfr':
    #     pass
        # Need to change all of this!!
        # temps, pres, mix, mdot, area, length = conds_dct
        # target_spcs = sim_util.get_target_spcs(exp_set, rename_instr)
        # raw_mech_result = np.ndarray((len(temps), len(target_spcs)))
        # for temp_idx, temp in enumerate(temps):
        #     sim_result, _, _ = reactors.pfr(  # don't need times or positions
        #         temp, pres, mix, gas, target_spcs, mdot, area, length)
        #     # Get last entry for all species, which is the outlet concentration
        #     raw_mech_result[temp_idx, :] = sim_result[-1, :]

    return raw_spcs, raw_pressures, raw_temps, uniform_times


def process_raw_result(raw_spcs, raw_pressures, raw_temps, uniform_times,
                       meas_type, conds_dct):

    if meas_type == 'conc':
        mech_result = raw_spcs  # no processing needs to occur
        uniform_times = uniform_times
    elif meas_type == 'idt':
        pass
    else:
        raise NotImplementedError(f'The meas_type {meas_type} is not working.')

    return mech_result, uniform_times


def jsr(conds_dct, gas):
    
    # Get arrays of reactor inputs
    temps = conds_dct['temp']
    pressures = conds_dct['pressure']
    mixes = conds_dct['mix']
    res_times = conds_dct['res_time']
    vols = conds_dct['vol']
    mdots = conds_dct['mdot']  # either res_time or mdot will be used
    target_spcs = conds_dct['target_spcs']

    # Initialize variables
    num_conds = conds_dct['num_conds']
    num_targets = len(target_spcs)
    raw_mech_result = np.ndarray((num_conds, num_targets))
    prev_concs = mixes[0]  # for first iter, use mix as prev_concs

    # Loop over the conditions, performing a simulation for each
    for cond_idx in range(num_conds):
        # mix = extract_mix_info(mixes[cond_idx])  # get rid of tuple
        temp = temps[cond_idx]
        pressure = pressures[cond_idx]
        mix = mixes[cond_idx]
        res_time = res_times[cond_idx]
        vol = vols[cond_idx]
        mdot = mdots[cond_idx]
        cond_result, prev_concs, _ = reactors.jsr(
            temp, pressure, mix, gas, target_spcs, res_time, vol, mdot=mdot,
            prev_concs=prev_concs)
        raw_mech_result[cond_idx, :] = cond_result

    return raw_mech_result
    

def st(conds_dct, gas):
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
    uniform_times = conds_dct['times']  # note: not used as a reactor input

    # Initialize variables
    num_conds = conds_dct['num_conds']
    num_spcs = len(target_spcs)  # note: this is num_spcs, not num_targets
    num_times = len(uniform_times)
    raw_spcs = np.ndarray((num_conds, num_spcs, num_times))
    # Note: below, the num_spcs dimension is missing since only one P and one T
    raw_pressures = np.ndarray((num_conds, num_times))
    raw_temps = np.ndarray((num_conds, num_times))

    # Perform a simulation for each condition
    for cond_idx in range(num_conds):
        temp = temps[cond_idx]
        pressure = pressures[cond_idx]
        mix = mixes[cond_idx]
        end_time = end_times[cond_idx]
        p_of_t = p_of_ts[cond_idx]
        cond_spcs, cond_pressures, cond_temps, cond_times, _, _ = reactors.st(
            temp, pressure, mix, gas, target_spcs, end_time, p_of_t=p_of_t)

        # Interpolate simulation results to fit them to the uniform time grid
        raw_spcs[cond_idx, :, :] = util.interp(
                cond_spcs, cond_times, uniform_times)
        raw_pressures[cond_idx, :] = util.interp(
            cond_pressures, cond_times, uniform_times)
        raw_temps[cond_idx, :] = util.interp(
            cond_temps, cond_times, uniform_times)

    return raw_spcs, raw_pressures, raw_temps, uniform_times


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


def _get_st_idt(target, times, method='baseline_extrap', plot=False):
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

    assert method in ALLOWED_IDT_METHODS, (
        f"IDT method should be one of {ALLOWED_IDT_METHODS}, not '{method}'")
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

    return idt, warnings
