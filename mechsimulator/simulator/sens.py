import numpy as np
from mechsimulator.simulator import outcome


def single_mech(conds_dct, gas, reac_type, meas_type, xdata, ydata_shape,
                factor=0.01):

    if reac_type == 'jsr':
        sens_coeffs = jsr(conds_dct, gas, meas_type, ydata_shape, factor=factor)
    elif reac_type == 'st':
        sens_coeffs = st(conds_dct, gas, meas_type, xdata, ydata_shape,
                         factor=factor)
    elif reac_type == 'pfr':
        sens_coeffs = pfr(conds_dct, gas, meas_type, xdata, ydata_shape,
                          factor=factor)
    elif reac_type == 'free_flame':
        sens_coeffs = free_flame(conds_dct, gas, meas_type, ydata_shape,
                                 factor=factor)
    elif reac_type == 'const_t_p':
        sens_coeffs = const_t_p(conds_dct, gas, meas_type, xdata, ydata_shape,
                                factor=factor)
    else:
        raise NotImplementedError(
            f"reac_type '{reac_type}' is not implemented for sens!")

    return sens_coeffs


def const_t_p(conds_dct, gas, meas_type, xdata, ydata_shape, factor=0.01):

    # Get the reference results (used as baseline for sens calcs)
    outcome_ydata_shape = ydata_shape[1:]  # strip first dim (# of rxns)
    ref_result = outcome.const_t_p(
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
        rxn_result = outcome.const_t_p(
            conds_dct, gas, meas_type, xdata, outcome_ydata_shape)
        sens_coeff = (rxn_result - ref_result) / (ref_result * factor)
        sens_coeffs[rxn_idx, :] = sens_coeff
        # If on last rxn in this mech, fill up missing rxns with NaNs
        if rxn_idx == num_rxns - 1:
            sens_coeffs[rxn_idx + 1:, :] = np.nan

    return sens_coeffs


def free_flame(conds_dct, gas, meas_type, ydata_shape, factor=0.01):

    # Get the reference results (used as baseline for sens calcs)
    outcome_ydata_shape = ydata_shape[1:]  # strip first dim (# of rxns)
    ref_result, prev_solns = outcome.free_flame(
        conds_dct, gas, meas_type, outcome_ydata_shape)

    # Get the sensitivity coefficients
    num_rxns = gas.n_reactions
    sens_coeffs = np.ndarray(ydata_shape)
    for rxn_idx in range(num_rxns):
        current_rxn = gas.reactions()[rxn_idx].equation
        print(f'Current rxn: {current_rxn}, {rxn_idx + 1} out of {num_rxns}')
        gas.set_multiplier(1.0)  # reset all multipliers to the original values
        gas.set_multiplier(1 + factor, rxn_idx)
        # Run simulation and calculate sens coefficients
        rxn_result, _ = outcome.free_flame(
            conds_dct, gas, meas_type, outcome_ydata_shape,
            prev_solns=prev_solns)
        sens_coeff = (rxn_result - ref_result) / (ref_result * factor)
        sens_coeffs[rxn_idx, :] = sens_coeff
        # If on last rxn in this mech, fill up missing rxns with NaNs
        if rxn_idx == num_rxns - 1:
            sens_coeffs[rxn_idx + 1:, :] = np.nan

    return sens_coeffs


def st(conds_dct, gas, meas_type, xdata, ydata_shape, factor=0.01):

    # Get the reference results (used as baseline for sens calcs)
    outcome_ydata_shape = ydata_shape[1:]  # strip first dim (# of rxns)
    ref_result = outcome.st(
        conds_dct, gas, meas_type, xdata, outcome_ydata_shape)

    # Get the sensitivity coefficients
    num_rxns = gas.n_reactions
    sens_coeffs = np.full(ydata_shape, np.nan)

    # Hardcoding for Jaeyoung's C2H6 simulations
    # hardcoded_rxn_idxs = [279, 284, 253, 285, 59, 60, 48, 204, 231, 214, 289, 276, 205,
    #             158, 254, 2, 3, 326, 152, 271, 262, 250, 156, 210, 568, 1]

    # Hardcoding for Curran C3H8 paper (using Luna's mech!)
    # hardcoded_rxn_idxs = [9,596,118,656,472,365,417,10,364,169,369,164,393,14]

    # Hardcoding for Curran C3H8 paper (using full NUIG1.1 mech!)
    hardcoded_rxn_idxs = [642,37,38,177,776,125,739,693,748,641,687,679,688]

    for rxn_idx in range(num_rxns):
        if rxn_idx in hardcoded_rxn_idxs:
            current_rxn = gas.reactions()[rxn_idx].equation
            print(f'Current rxn: {current_rxn}, {rxn_idx + 1} out of {num_rxns}')
            gas.set_multiplier(1.0)  # reset all multipliers to the original values
            gas.set_multiplier(1 + factor, rxn_idx)
            # Run simulation and calculate sens coefficients
            rxn_result = outcome.st(
                conds_dct, gas, meas_type, xdata, outcome_ydata_shape)
            sens_coeff = (rxn_result - ref_result) / (ref_result * factor)
            sens_coeffs[rxn_idx, :] = sens_coeff
            # # If on last rxn in this mech, fill up missing rxns with NaNs
            # if rxn_idx == num_rxns - 1:
            #     sens_coeffs[rxn_idx + 1:, :] = np.nan

    return sens_coeffs


def jsr(conds_dct, gas, meas_type, ydata_shape, factor=0.01):
    """ Calculates sensitivity coefficients for any number of target species
        across a range of conditions for a JSR simulation

    """

    # Get the reference results (used as baseline for sens calcs) and the
    # prev_soln, which helps converge to a solution
    outcome_ydata_shape = ydata_shape[1:]  # strip first dim (# of rxns)
    ref_result, prev_solns = outcome.jsr(
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
        rxn_result, _ = outcome.jsr(conds_dct, gas, meas_type,
                                    outcome_ydata_shape, prev_solns=prev_solns)
        sens_coeff = (rxn_result - ref_result) / (ref_result * factor)
        sens_coeffs[rxn_idx, :] = sens_coeff
        # If on last rxn in this mech, fill up missing rxns with NaNs
        if rxn_idx == num_rxns - 1:
            sens_coeffs[rxn_idx + 1:, :] = np.nan

    return sens_coeffs


def pfr(conds_dct, gas, meas_type, xdata, ydata_shape, factor=0.01):
    """ Calculates sensitivity coefficients for any number of target species
        across a range of temperatures for a PFR simulation

        :param factor: factor by which to perturb A factors
        :type factor: float
        :return sens_coeffs: sens coeffs for each rxn, temp, and species
        :rtype: Numpy array of shape (num_rxns, num_temps, num_spcs)
    """

    # Get the reference results (used as baseline for sens calcs) and the
    # prev_soln, which helps converge to a solution
    outcome_ydata_shape = ydata_shape[1:]  # strip first dim (# of rxns)
    ref_result = outcome.pfr(conds_dct, gas, meas_type, xdata,
                             outcome_ydata_shape)

    # Get the sensitivity coefficients
    num_rxns = gas.n_reactions  # might differ from sens_coeffs for this mech
    sens_coeffs = np.ndarray(ydata_shape, dtype=float)  # helps with NaNs
    for rxn_idx in range(num_rxns):
        current_rxn = gas.reactions()[rxn_idx].equation
        print(f'Current rxn: {current_rxn}, {rxn_idx + 1} out of {num_rxns}')
        gas.set_multiplier(1.0)  # reset all multipliers to the original values
        gas.set_multiplier(1 + factor, rxn_idx)
        # Run simulation and calculate sens coefficients
        rxn_result = outcome.pfr(
            conds_dct, gas, meas_type, xdata, outcome_ydata_shape)
        sens_coeff = (rxn_result - ref_result) / (ref_result * factor)
        sens_coeffs[rxn_idx, :] = sens_coeff
        # If on last rxn in this mech, fill up missing rxns with NaNs
        if rxn_idx == num_rxns - 1:
            sens_coeffs[rxn_idx + 1:, :] = np.nan

    return sens_coeffs
