import numpy as np
from mechsimulator.sim import util
from mechsimulator.parser.exp import get_poss_inps
from mechsimulator.sim import outcome


def mult_sets(exp_sets, gases, spc_ident_dcts,  calc_type, plot_or_exps):
    """ Runs multiple experimental sets against multiple mechanisms

        :param exp_sets:
        :param gases:
        :param spc_ident_dcts:
        :param calc_type:
        :param plot_or_exps:
        :return:
    """

    set_results = []  # using lists since array dimensions/sizes may vary
    set_times = []
    for exp_set in exp_sets:
        set_result, set_time = mult_mechs(exp_set, gases, spc_ident_dcts,
                                          calc_type, plot_or_exps)
        set_results.append(set_result)
        set_times.append(set_time)

    return set_results, set_times


def mult_mechs(exp_set, gases, mech_spc_dcts, calc_type, plot_or_exps):
    """ Plots the predictions of any number of mechanisms against a single
        exp_set

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

    assert plot_or_exps in ('plot', 'exps'), (
        f"'plot_or_exps' should be 'plot' or 'exps', not '{plot_or_exps}'")

    # Initialize array
    mech_result_shape = util.get_mech_result_shape(exp_set, plot_or_exps)
    num_mechs = len(gases)
    mech_results_shape = (num_mechs,) + mech_result_shape  # prepend mechs_dim
    mech_results = np.ndarray(mech_results_shape)

    # Loop over each mech and run a simulation with each
    for mech_idx, gas in enumerate(gases):
        mech_spc_dct = mech_spc_dcts[mech_idx]
        set_result, uniform_times = single_mech(exp_set, gas, mech_spc_dct,
                                                calc_type, plot_or_exps)
        mech_results[mech_idx] = set_result

    return mech_results, uniform_times


def single_mech(exp_set, gas, mech_spc_dct, calc_type, plot_or_exps):

    # Load the conditions dict
    if plot_or_exps == 'plot':
        conds_dct = get_conds_dct_plot(exp_set, mech_spc_dct)
    else:  # exps
        conds_dct = get_conds_dct_exps(exp_set, mech_spc_dct)

    # Run the desired calculation type
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']
    if calc_type == 'outcome':
        mech_result = outcome.get_result(conds_dct, gas, reac_type, meas_type)
    elif calc_type == 'sens':
        print('sens not ready!')
    elif calc_type == 'rop':
        print('rop not ready!')
    elif calc_type == 'pathways':
        print('pathways not ready!')

    return mech_result


def get_conds_dct_plot(exp_set, mech_spc_dct):
    """ Extracts the plotting conditions from an exp_set

        :param exp_set: object describing a set of experiments
        :type exp_set: dct
        :return conds_dct: plotting conditions; the plotting variable is a
            Numpy array of shape (num_values,) while other entries are floats
        :rtype: dct ['temp': temp(s), 'pressure': pressure(s), ...]
    """

    # Load data from the exp_set
    plot_dct = exp_set['plot']
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']
    var = util.get_plot_variable(exp_set)
    conds = util.get_plot_conds(exp_set)
    num_conds = len(conds)

    # Get the list of possible inputs for the reactor/measurement combo
    poss_inps = get_poss_inps(reac_type, meas_type, 'plot', remove_bad=True)

    # Loop over possible inputs and build the conditions dict
    conds_dct = {}
    for poss_inp in poss_inps:
        # If on the plotting variable, use the array of plotting conditions
        if poss_inp == var:
            conds_dct[poss_inp] = conds
        # If single value exists, create an array of repeats of the single value
        elif plot_dct.get(poss_inp) is not None:
            conds_dct[poss_inp] = [plot_dct[poss_inp]] * num_conds
        # If value does not exist in plot_dct, create array of Nones
        else:
            conds_dct[poss_inp] = [None] * num_conds

    # Add the mix information since it is stored separately from the plot info
    raw_mix = exp_set['mix']
    exp_spc_dct = exp_set['spc']
    rename_instr = util.get_rename_instr(mech_spc_dct, exp_spc_dct)
    mix = util.translate_spcs(raw_mix, rename_instr)
    conds_dct['mix'] = [mix] * num_conds
    # Note: if I add the ability to have the mix be varied, would need to have a
    # translator function that would convert the array of mix conds to
    # dcts defining the mix

    # Add the target_spcs info since it is stored separately
    raw_target_spcs = tuple(exp_spc_dct.keys())
    target_spcs = util.translate_spcs(raw_target_spcs, rename_instr)
    conds_dct['target_spcs'] = target_spcs

    # Add miscellaneous information
    conds_dct['num_conds'] = num_conds
    if util.check_time_resolved(exp_set):  # if time resolved
        conds_dct['times'] = util.get_uniform_times(exp_set)
        conds_dct['timestep'] = exp_set['plot']['timestep']  # might be unneeded

    return conds_dct


def get_conds_dct_exps(exp_set, mech_spc_dct):
    """

        :param exp_set:
        :param mech_spc_dct:
        :return conds_dct:
        :rtype: dct
    """

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
    exp_objs = exp_set['exp_objs']
    num_conds = len(exp_objs)  # each experiment is its own condition

    # Get the list of possible inputs for the reactor/measurement combo
    poss_inps = get_poss_inps(reac_type, meas_type, 'exps', remove_bad=True)

    # Loop over possible inputs and build the conditions dict
    conds_dct = {}
    for poss_inp in poss_inps:
        conds_dct[poss_inp] = [None] * num_conds  # initialize as Nones
        for cond_idx, exp_obj in enumerate(exp_objs):
            if exp_obj['conds'].get(poss_inp) is not None:
                conds_dct[poss_inp][cond_idx] = exp_obj['conds'][poss_inp][0]
            # Otherwise, just leave as None

    # Add the mix information since it is stored separately from other conds
    exp_spc_dct = exp_set['spc']
    rename_instr = util.get_rename_instr(mech_spc_dct, exp_spc_dct)
    conds_dct['mix'] = [None] * num_conds
    for cond_idx, exp_obj in enumerate(exp_objs):
        raw_mix = exp_obj['mix']
        # Translate to the mechanism spcs names
        tuple_mix = util.translate_spcs(raw_mix, rename_instr)
        mix = extract_mix_info(tuple_mix)  # remove the uncertainty bounds
        conds_dct['mix'][cond_idx] = mix

    # Add the target_spcs info since it is stored separately
    raw_target_spcs = tuple(exp_spc_dct.keys())
    target_spcs = util.translate_spcs(raw_target_spcs, rename_instr)
    conds_dct['target_spcs'] = target_spcs

    # Add the P(t) information since it is stored separately
    conds_dct['p_of_t'] = [None] * num_conds
    for cond_idx, exp_obj in enumerate(exp_objs):
        times = exp_obj['result'].get('time')[0]  # [0] omits uncertainty bounds
        pressures = exp_obj['result'].get('pressure')[0]
        p_of_t = np.vstack((times, pressures))
        conds_dct['p_of_t'][cond_idx] = p_of_t

    # Add miscellaneous information
    conds_dct['num_conds'] = num_conds
    if util.check_time_resolved(exp_set):  # if time resolved
        conds_dct['times'] = util.get_uniform_times(exp_set)
        conds_dct['timestep'] = exp_set['plot']['timestep']

    return conds_dct


