import numpy as np
import copy
from mechsimulator.sim import reactors
from mechsimulator.plotter import util


def mult_sets_plotting(exp_sets, gases, spc_ident_dcts):

    sim_results = []
    for exp_set in exp_sets:
        sim_result = mult_mechs_plotting(exp_set, gases, spc_ident_dcts)
        sim_results.append(sim_result)

    return sim_results


def mult_mechs_plotting(exp_set, gases, mech_spc_dcts):
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
            (i) time-resolved: (num_mechs,num_conds, num_targets, num_timesteps)
            (ii) non-time-resolved: (num_mechs, num_conds, num_targets)
    """

    # Initialize array based on experiment type
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']
    if reac_type == 'jsr' or reac_type == 'pfr':
        temps = util.get_plotting_variable(exp_set)
        target_spcs = list(exp_set['spc'].keys())
        mech_results = np.ndarray((len(gases), len(temps), len(target_spcs)))
    else:
        raise NotImplementedError

    # Loop over each mech
    for gas_idx, gas in enumerate(gases):
        mech_spc_dct = mech_spc_dcts[gas_idx]
        set_result = single_mech_plotting(exp_set, gas, mech_spc_dct)
        mech_results[gas_idx] = set_result

    return mech_results


def single_mech_plotting(exp_set, gas, mech_spc_dct):
    """ Gets results for a single exp_set using a single mechanism

        :param exp_set: object describing a set of experiments
        :type exp_set: dct
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param mech_spc_dct: identifying information on species in a mech
        :type: dct {spc1: spc_dct1, spc2: ...}
        :return mech_results: sim results for an exp_set using one mechanism
        :rtype: Numpy array whose shape depends on the type of measurement:
            (i) time-resolved: (num_conditions, num_targets, num_timesteps)
            (ii) non-time-resolved: (num_conditions, num_targets)
    """

    # Get the rename_instr
    exp_spc_dct = exp_set['spc']
    rename_instr = get_rename_instr(mech_spc_dct, exp_spc_dct)

    # Run the simulations based on experiment type
    reac_type = exp_set['overall']['reac_type']
    plot_conds = get_plot_conds(exp_set, rename_instr)
    if reac_type == 'jsr':
        temps, pres, mix, res_time, vol = plot_conds
        target_spcs = get_target_spcs(exp_set, rename_instr)
        prev_concs = copy.deepcopy(mix)  # for first iter, use mix as prev_concs
        mech_result = np.ndarray((len(temps), len(target_spcs)))
        for temp_idx, temp in enumerate(temps):
            sim_result, prev_concs = reactors.jsr(
                temp, pres, mix, gas, target_spcs, res_time, vol, prev_concs)
            mech_result[temp_idx, :] = sim_result
    elif reac_type == 'pfr':
        temps, pres, mix, mdot, area, length = plot_conds
        target_spcs = get_target_spcs(exp_set, rename_instr)
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
        :return plot_conds: list of plotting conditions; number of entries
            depends on the reactor type
        :rtype: list
    """

    plot_dct = exp_set['plot']
    raw_mix = exp_set['mix']
    reac_type = exp_set['overall']['reac_type']

    if reac_type == 'jsr':
        temps = util.get_plotting_variable(exp_set)
        pressure = plot_dct['pressure']
        res_time = plot_dct['res_time']
        vol = plot_dct['vol']
        mix = translate_spcs(raw_mix, rename_instr)
        conditions = [temps, pressure, mix, res_time, vol]

    elif reac_type == 'pfr':
        temps = util.get_plotting_variable(exp_set)
        pressure = plot_dct['pressure']
        mdot = plot_dct['mdot']
        area = plot_dct['area']
        length = plot_dct['length']
        mix = translate_spcs(raw_mix, rename_instr)
        conditions = [temps, pressure, mix, mdot, area, length]

    else:
        raise NotImplementedError

    return conditions


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
        file with those in the mechanism

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
    else:
        translated_names = []
        for spc in names:
            renamed_spc = rename_instr[spc]
            translated_names.append(renamed_spc)
        if isinstance(names, tuple):  # change to tuple if needed
            translated_names = tuple(translated_names)

    return translated_names
