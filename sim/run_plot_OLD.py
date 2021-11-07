import numpy as np
import copy
from mechsimulator.sim import reactors
from mechsimulator.plotter import util as plot_util
from mechsimulator.sim import util as sim_util
from mechsimulator.parser.exp import ALLOWED_REAC_TYPES


def mult_sets(exp_sets, gases, spc_ident_dcts):
    sim_results = []
    for exp_set in exp_sets:
        sim_result = mult_mechs(exp_set, gases, spc_ident_dcts)
        sim_results.append(sim_result)

    return sim_results


def mult_mechs(exp_set, gases, mech_spc_dcts):
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

    # Initialize array based on experiment type
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']
    if reac_type == 'jsr' or reac_type == 'pfr':
        _, vals = plot_util.get_plotting_variable(exp_set)
        num_mechs = len(gases)
        num_conds = len(vals)
        num_targets = len(list(exp_set['spc'].keys()))
        mech_results = np.ndarray((num_mechs, num_conds, num_targets))
    else:
        raise NotImplementedError

    # Loop over each mech
    for mech_idx, gas in enumerate(gases):
        mech_spc_dct = mech_spc_dcts[mech_idx]
        set_result = single_mech(exp_set, gas, mech_spc_dct)
        mech_results[mech_idx] = set_result

    return mech_results


def single_mech(exp_set, gas, mech_spc_dct):
    """ Gets results for a single exp_set using a single mechanism

        :param exp_set: object describing a set of experiments
        :type exp_set: dct
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param mech_spc_dct: identifying information on species in a mech
        :type: dct {spc1: spc_dct1, spc2: ...}
        :return mech_results: sim results for an exp_set using one mechanism
        :rtype: Numpy array whose shape depends on the type of measurement:
            (i) time-resolved: (num_conds, num_targets, num_times)
            (ii) non-time-resolved: (num_conds, num_targets)
    """

    # Get the rename_instr
    exp_spc_dct = exp_set['spc']
    rename_instr = sim_util.get_rename_instr(mech_spc_dct, exp_spc_dct)

    # Run the simulations based on experiment type
    reac_type = exp_set['overall']['reac_type']
    plot_conds = get_plot_conds(exp_set, rename_instr)
    if reac_type == 'jsr':
        # Get reactor inputs
        temp = plot_conds['temp']
        pressure = plot_conds['pressure']
        mix = plot_conds['mix']
        res_time = plot_conds['res_time']
        vol = plot_conds['vol']
        mdot = plot_conds['mdot']  # either res_time or mdot will be used
        target_spcs = sim_util.get_target_spcs(exp_set, rename_instr)
        # Initialize variables
        num_conds = len(temp)  # could have used any of the array variables
        num_targets = len(target_spcs)
        mech_result = np.ndarray((num_conds, num_targets))
        prev_concs = copy.deepcopy(mix)  # for first iter, use mix as prev_concs
        # Loop over the conditions, performing a simulation for each
        for cond_idx in range(num_conds):
            sim_result, prev_concs, _ = reactors.jsr(
                temp[cond_idx], pressure[cond_idx], mix, gas, target_spcs,
                res_time[cond_idx], vol[cond_idx], mdot=mdot[cond_idx],
                prev_concs=prev_concs)
            mech_result[cond_idx, :] = sim_result

    elif reac_type == 'pfr':
        temps, pres, mix, mdot, area, length = plot_conds
        target_spcs = sim_util.get_target_spcs(exp_set, rename_instr)
        mech_result = np.ndarray((len(temps), len(target_spcs)))
        for temp_idx, temp in enumerate(temps):
            sim_result, _, _ = reactors.pfr(  # don't need times or positions
                temp, pres, mix, gas, target_spcs, mdot, area, length)
            # Get last entry for all species, which is the outlet concentration
            mech_result[temp_idx, :] = sim_result[-1, :]
    else:
        raise NotImplementedError

    return mech_result


def get_plot_conds(exp_set, rename_instr):
    """ Extracts the plotting conditions from an exp_set

        :param exp_set: object describing a set of experiments
        :type exp_set: dct
        :param rename_instr: spc_names in the experiment to be renamed according
            to the mechanism spc names
        :type rename_instr: dct {spc_to_be_renamed: new_spc_name, ...}
        :return plot_conds: plotting conditions; the plotting variable is a
            Numpy array of shape (num_values,) while other entries are floats
        :rtype: dct ['temp': temp(s), 'pressure': pressure(s), ...]
    """

    # Load initial data
    plot_dct = exp_set['plot']
    print('plot_dct:\n', plot_dct)
    raw_mix = exp_set['mix']
    reac_type = exp_set['overall']['reac_type']

    # Loop over the required inputs for the reactor type
    plot_conds = {}
    var, conds = plot_util.get_plotting_variable(exp_set)
    num_conds = len(conds)
    for req_inp in ALLOWED_REAC_TYPES[reac_type]:
        # If req_inp is a str, it's truly required
        if isinstance(req_inp, str):
            # If on the plotting variable, use the array that describes it
            if req_inp == var:
                plot_conds[req_inp] = conds
            # Otherwise, create an array of repeats of the single value
            else:
                plot_conds[req_inp] = np.ones(num_conds) * plot_dct[req_inp]
        # If req_inp is a tuple, it's contains a few optional inputs, one of
        # which must be provided
        else:
            for opt_inp in req_inp:
                # If on the plotting variable, use the array that describes it
                if opt_inp == var:
                    plot_conds[opt_inp] = conds
                # If on entry with a value, create array of repeats of the value
                elif plot_dct.get(opt_inp) is not None:
                    plot_conds[opt_inp] = np.ones(num_conds) * plot_dct[opt_inp]
                # If on entry with no value, create a list of Nones
                elif plot_dct.get(opt_inp) is None:
                    plot_conds[opt_inp] = [None] * num_conds

    # Add the mix information since it is stored separately from the plot info
    mix = sim_util.translate_spcs(raw_mix, rename_instr)
    plot_conds['mix'] = mix
    # Note: this currently only inputs a single value for 'mix'; namely, a dct.
    # If I add the ability to have the mix be varied, would need to have a
    # translator function that would convert the array of mix conds to
    # dcts defining the mix. Then, could add a list of these dcts in the 'mix'
    # field of the plot_conds dct. Would need to update the call of mix in
    # single_mech to have indexing.

    return plot_conds
