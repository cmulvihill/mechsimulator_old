import numpy as np
from mechsimulator.simulator import reactors
from mechsimulator.simulator import outcome


def get_result(conds_dct, gas, reac_type, meas_type, xdata, ydata_shape,
               factor=0.01):

    if reac_type == 'jsr':
        sens_coeffs = jsr(conds_dct, gas, meas_type, ydata_shape, factor=factor)
    elif reac_type == 'st':
        sens_coeffs = st(conds_dct, gas, meas_type, xdata, ydata_shape,
                         factor=factor)
    else:
        raise NotImplementedError(
            f"meas_type '{meas_type}' is not implemented for sens!")

    return sens_coeffs


def st(conds_dct, gas, meas_type, xdata, ydata_shape, factor=0.01):

    # Get the reference results (used as baseline for sens calcs)
    outcome_ydata_shape = ydata_shape[1:]  # strip first dim (# of rxns)
    ref_result = outcome.st(
        conds_dct, gas, meas_type, xdata, outcome_ydata_shape)

    # Get the sensitivity coefficients
    num_rxns = gas.n_reactions
    sens_coeffs = np.ndarray(ydata_shape)
    for rxn_idx in range(num_rxns):
        current_rxn = gas.reactions()[rxn_idx].equation
        print(f'Current rxn: {current_rxn}, {rxn_idx + 1} out of {num_rxns}')
        gas.set_multiplier(1.0)  # reset all multipliers to the original values
        gas.set_multiplier(1 + factor, rxn_idx)
        # Run simulation and calculate sens coefficients
        rxn_result = outcome.st(
            conds_dct, gas, meas_type, xdata, outcome_ydata_shape)
        sens_coeff = (rxn_result - ref_result) / (ref_result * factor)
        sens_coeffs[rxn_idx, :] = sens_coeff
        # If on last rxn in this mech, fill up missing rxns with NaNs
        if rxn_idx == num_rxns - 1:
            sens_coeffs[rxn_idx + 1:, :] = np.nan

    return sens_coeffs


def jsr(conds_dct, gas, meas_type, ydata_shape, factor=0.01):
    """ Calculates sensitivity coefficients for any number of target species
        across a range of conditions for a JSR simulation

        :return sens_coeffs: sens coeffs for each rxn, temp, and target species
        :rtype: Numpy array of shape ydata_shape
        :return ref_result: the baseline experimental outcome
        :rtype: Numpy array of shape ydata_shape, except with the first
            dimension removed
    """

    # Get the reference results (used as baseline for sens calcs) and the
    # prev_soln, which helps converge to a solution
    outcome_ydata_shape = ydata_shape[1:]  # strip first dim (# of rxns)
    ref_result, prev_soln = outcome.jsr(
        conds_dct, gas, meas_type, outcome_ydata_shape, output_all=True)

    # Get the sensitivity coefficients
    num_rxns = gas.n_reactions  # might differ from sens_coeffs for this mech
    sens_coeffs = np.ndarray(ydata_shape, dtype=float)  # helps with NaNs
    for rxn_idx in range(num_rxns):
        current_rxn = gas.reactions()[rxn_idx].equation
        print(f'Current rxn: {current_rxn}, {rxn_idx + 1} out of {num_rxns}')
        gas.set_multiplier(1.0)  # reset all multipliers to the original values
        gas.set_multiplier(1 + factor, rxn_idx)
        # Run simulation and calculate sens coefficients
        rxn_result, _ = outcome.jsr(
            conds_dct, gas, meas_type, outcome_ydata_shape, prev_soln=prev_soln)
        sens_coeff = (rxn_result - ref_result) / (ref_result * factor)
        sens_coeffs[rxn_idx, :] = sens_coeff
        # If on last rxn in this mech, fill up missing rxns with NaNs
        if rxn_idx == num_rxns - 1:
            sens_coeffs[rxn_idx + 1:, :] = np.nan

    return sens_coeffs


def pfr(temps, pressure, mix, gas, target_spcs, mdot, area, length, factor=0.01):
    """ Calculates sensitivity coefficients for any number of target species
        across a range of temperatures for a PFR simulation

        :param temps: temperatures at which to perform calculations (K)
        :type temps: Numpy array of shape (num_temps,)
        :param pressure: reactor constant pressure (atm)
        :type pressure: float
        :param mix: initial species concentrations at reactor inlet
        :type mix: dct {spc1: conc1, spc2: ...}
            (example: {'H2': 0.5, 'O2': 0.5})
        :param gas: a Cantera object describing a kinetics mechanism/gas
        :type gas: Cantera Solution object
        :param target_spcs: desired species concentrations
        :type target_spcs: list [spc1, spc2, ...]
        :param mdot: reactor mass flow rate (kg/s)
        :type mdot: float
        :param area: reactor cross-sectional area (m^2)
        :type area: float
        :param length: reactor length (m)
        :type length: float
        :param factor: factor by which to perturb A factors
        :type factor: float
        :return sens_coeffs: sens coeffs for each rxn, temp, and species
        :rtype: Numpy array of shape (num_rxns, num_temps, num_spcs)
    """

    # Get the reference results (used as baseline for sens calcs)
    ref_results = np.ndarray((len(temps), len(target_spcs)))
    for temp_idx, temp in enumerate(temps):
        ref_result, _, _ = reactors.pfr(temp, pressure, mix, gas, target_spcs,
                                        mdot, area, length)
        # ref_results has shape (num_spcs, num_times); the last entry (i.e., the
        # reactor outlet) is the value to be stored in ref_results
        ref_results[temp_idx, :] = ref_result[:, -1]

    # Get the sensitivity coefficients
    num_rxns = gas.n_reactions
    sens_coeffs = np.ndarray((num_rxns, len(temps), len(target_spcs)))
    for rxn_idx in range(num_rxns):
        current_rxn = gas.reactions()[rxn_idx].equation
        print(f'Current rxn: {current_rxn}, {rxn_idx + 1} out of {num_rxns}')
        gas.set_multiplier(1.0)  # reset the multipliers to the original values
        gas.set_multiplier(1 + factor, rxn_idx)
        for temp_idx, temp in enumerate(temps):
            results, _, _ = reactors.pfr(temp, pressure, mix, gas, target_spcs,
                                         mdot, area, length)
            result = results[:, -1]  # get the last (outlet) concentration
            ref = ref_results[temp_idx, :]
            sens_coeff = (result - ref) / (ref * factor)
            sens_coeffs[rxn_idx, temp_idx, :] = sens_coeff

    return sens_coeffs, ref_results
