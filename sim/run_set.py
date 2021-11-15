import copy
import numpy as np
from mechsimulator.sim import util
from mechsimulator.parser.exp_checker import get_poss_inps
from mechsimulator.sim import outcome


def mult_sets(exp_sets, gases, mech_spc_dcts,  calc_type, x_source,
              conds_source, sim_opts_lst=None):
    """ Runs multiple experimental sets against multiple mechanisms

        :param exp_sets:
        :param gases:
        :param mech_spc_dcts:
        :param calc_type:
        :param x_source: the source of the x data (e.g., temperature, time) used
            in the simulations; either 'plot' or 'exps'
        :param conds_source: the source of the conditions (e.g., temperature)
            used in the simulations; either 'plot' or 'exps'
        :param sim_opts_lst: list containing sim_options for
        :type sim_opts_lst: list
        :return:
    """

    # Check inputs
    check_source_inp(x_source, conds_source)

    # Loop over each experimental set
    set_ydata_lst = []  # using lists since array dimensions/sizes may vary
    set_xdata_lst = []
    for exp_set in exp_sets:
        set_ydata, set_xdata = mult_mechs(exp_set, gases, mech_spc_dcts,
                                          calc_type, x_source, conds_source,
                                          sim_opts_lst=sim_opts_lst)
        set_ydata_lst.append(set_ydata)
        set_xdata_lst.append(set_xdata)

    return set_ydata_lst, set_xdata_lst


def mult_mechs(exp_set, gases, mech_spc_dcts, calc_type, x_source,
               conds_source, sim_opts_lst=None):
    """ Calculates the predictions of any number of mechanisms at the conditions
        of a single exp_set

        :param exp_set: object describing a set of experiments
        :type exp_set: dct
        :param gases: Cantera objects describing kinetics mechanisms/gases
        :type gases: list
        :param mech_spc_dcts: species identifying info for each mechanisms
        :type mech_spc_dcts: list
        :return mech_results: sim results for one exp_set using many mechanisms
        :rtype: Numpy array whose shape depends on the type of measurement:
            (i) time-resolved: (num_mechs, num_conds, num_targets, num_times)
            (ii) non-time-resolved: (num_mechs, num_conds, num_targets)
    """

    # Initialize variables
    mech_ydata_shape, _, _, mech_xdata = util.get_mech_info(
        exp_set, calc_type, x_source, conds_source, gases)
    set_xdata = mech_xdata  # xdata is the same for all mechs in the set
    num_mechs = len(gases)
    set_ydata_shape = (num_mechs,) + mech_ydata_shape  # prepend mechs_dim
    set_ydata = np.ndarray(set_ydata_shape)

    # Check that the sim_opts input was correctly used
    assert sim_opts_lst is None or len(sim_opts_lst) == num_mechs, (
        'Unless None, the sim_opts_lst input should be a list of length '
        'num_mechs, {num_mechs}, not {len(sim_opts)}')

    # Loop over each mech and run a simulation with each
    for mech_idx, gas in enumerate(gases):
        # Get the mech_opts for each mech
        if sim_opts_lst is None:
            sim_opts = None  # will return a copy of exp_set
        else:
            sim_opts = sim_opts_lst[mech_idx]
        # Get the updated_exp_set for each mech
        updated_exp_set = update_exp_set(exp_set, sim_opts=sim_opts)
        # Run the simulation and store the results
        mech_spc_dct = mech_spc_dcts[mech_idx]
        mech_ydata = single_mech(updated_exp_set, gas, mech_spc_dct, calc_type,
                                 conds_source, mech_xdata, mech_ydata_shape)
        set_ydata[mech_idx] = mech_ydata

    return set_ydata, set_xdata


def single_mech(exp_set, gas, mech_spc_dct, calc_type, conds_source, xdata,
                ydata_shape):

    # Must get conds_dct for each mech since mix spc names can change
    conds_dct = get_conds_dct(exp_set, mech_spc_dct, conds_source)

    # Run the desired calculation type
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']
    if calc_type == 'outcome':
        mech_result = outcome.get_result(
            conds_dct, gas, reac_type, meas_type, xdata, ydata_shape)
    elif calc_type == 'sens':
        print('sens not ready!')
    elif calc_type == 'rop':
        print('rop not ready!')
    elif calc_type == 'pathways':
        print('pathways not ready!')

    return mech_result


def get_conds_dct(exp_set, mech_spc_dct, conds_source):

    def extract_mix_info(inp_mix):
        """ Removes the uncertainty bounds from a mixture

            :param inp_mix: mixture with uncertainty bounds
            :type inp_mix: dct {spc1: (conc1, lb1, ub1), spc2: ...}
            :return new_mix: mixture with uncertainty bounds removed
            :rtype: dct {spc1: conc1, spc2: ...}
        """

        new_mix = {}
        for spc, conc_tup in inp_mix.items():
            new_mix[spc] = conc_tup[0]

        return new_mix

    # Load initial data
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']
    if conds_source == 'plot':
        plot_dct = exp_set['plot']
        var = util.get_plot_variable(exp_set)
        conds = util.get_plot_conds(exp_set)
        num_conds = len(conds)
    else:  # 'exps'
        exp_objs = exp_set['exp_objs']
        num_conds = len(exp_objs)  # each experiment is its own condition

    # Loop over possible inputs and build the conditions dict
    conds_dct = {}
    poss_inps = get_poss_inps(reac_type, meas_type, conds_source, rm_bad=True)
    for poss_inp in poss_inps:
        conds_dct[poss_inp] = [None] * num_conds  # initialize as Nones
        if conds_source == 'plot':
            if poss_inp == var:  # if on plotting variable, use array
                conds_dct[poss_inp] = conds
            elif plot_dct.get(poss_inp) is not None:  # if a single value exists
                conds_dct[poss_inp] = [plot_dct[poss_inp]] * num_conds
            # Otherwise, just leave as Nones
        else:  # 'exps'
            for cond_idx, exp_obj in enumerate(exp_objs):
                if exp_obj['conds'].get(poss_inp) is not None:
                    conds_dct[poss_inp][cond_idx] = exp_obj['conds'][poss_inp][
                        0]
                # Otherwise, just leave as None

    # Add the mix information since it is stored separately
    exp_spc_dct = exp_set['spc']
    rename_instr = util.get_rename_instr(mech_spc_dct, exp_spc_dct)
    if conds_source == 'plot':
        raw_mix = exp_set['mix']
        mix = util.translate_spcs(raw_mix, rename_instr)
        conds_dct['mix'] = [mix] * num_conds
        # Note: if I add the ability to have the mix be varied, would need to
        # have a translator function that would convert the array of mix conds
        # to dcts defining the mix
    else:  # 'exps'
        conds_dct['mix'] = [None] * num_conds
        for cond_idx, exp_obj in enumerate(exp_objs):
            raw_mix = exp_obj['mix']
            tuple_mix = util.translate_spcs(raw_mix, rename_instr)
            mix = extract_mix_info(tuple_mix)  # remove uncertainty bounds
            conds_dct['mix'][cond_idx] = mix

    # Add the target_spcs info since it is stored separately
    raw_target_spcs = tuple(exp_spc_dct.keys())
    target_spcs = util.translate_spcs(raw_target_spcs, rename_instr)
    conds_dct['target_spcs'] = target_spcs

    # Add the P(t) information since it is stored separately
    if reac_type == 'st':
        conds_dct['p_of_t'] = [None] * num_conds
        if conds_source == 'plot':
            if 'dpdt' in plot_dct:
                end_time = plot_dct['end_time']
                dpdt = plot_dct['dpdt']
                p_of_t = create_p_of_t(end_time, dpdt)
                conds_dct['p_of_t'] = [p_of_t] * num_conds
        else:  # 'exps'
            for cond_idx, exp_obj in enumerate(exp_objs):
                if 'pressure' in exp_obj['result']:
                    times = exp_obj['result'].get('time')[0]
                    pressures = exp_obj['result'].get('pressure')[0]
                    p_of_t = np.vstack((times, pressures))
                    conds_dct['p_of_t'][cond_idx] = p_of_t
                elif 'dpdt' in exp_obj['conds']:
                    end_time = exp_obj['conds']['end_time'][0]
                    dpdt = exp_obj['conds']['dpdt'][0]
                    p_of_t = create_p_of_t(end_time, dpdt)
                    conds_dct['p_of_t'][cond_idx] = p_of_t
                # Otherwise, leave as None

    # Add ignition delay time information
    if meas_type == 'idt':
        conds_dct['idt_method'] = plot_dct['idt_method']
        raw_idt_targets = plot_dct['idt_target']
        idt_targets = []
        for raw_idt_target in raw_idt_targets:
            if raw_idt_target == 'pressure':
                idt_target = 'pressure'
            else:  # if a species
                idt_target = util.translate_spcs([raw_idt_target], rename_instr)
                assert idt_target is not None, (
                    f'The species {raw_idt_target} is an idt target but is not '
                    f'defined in the mechanism.')
                idt_target = idt_target[0]
            idt_targets.append(idt_target)
        conds_dct['idt_target'] = idt_targets

    return conds_dct


def create_p_of_t(end_time, dpdt):
    """ Creates a P(t) array from an end_time and a dP/dt

        :param end_time: final time of the simulation (seconds)
        :param dpdt: linear change in pressure during the simulation (%/ms)
        :return p_of_t: array of noramlized pressure starting at 1 and going
            to the end pressure calculated by the end_time and dpdt
        :rtype: Numpy array of shape (2,
    """

    times = np.array([0, end_time])
    end_pressure = 1 + ((end_time * 1e3) * dpdt) / 100
    pressures = np.array([1, end_pressure])
    p_of_t = np.vstack((times, pressures))

    return p_of_t


def check_source_inp(x_source, conds_source):
    """ Checks the inputs for the x_source and conds_source

    :param x_source:
    :param conds_source:
    :return:
    """

    assert x_source in ('plot', 'exps'), (
        f"'x_source' should be 'plot' or 'exps', not '{x_source}'")
    assert conds_source in ('plot', 'exps'), (
        f"'conds_source' should be 'plot' or 'exps', not '{conds_source}'")
    if x_source == 'exps':
        assert conds_source != 'plot', (
            f"The 'x_source' and 'conds_source' combination 'exps' and 'plot',"
            " respectively, is not allowed.")


def update_exp_set(exp_set, sim_opts=None):
    """ Updates an exp_set

        :param exp_set:
        :param sim_opts:
        :return:
    """

    # Doesn't actually do anything right now

    updated_exp_set = copy.deepcopy(exp_set)

    return updated_exp_set