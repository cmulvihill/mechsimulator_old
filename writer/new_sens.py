import pandas as pd
import numpy as np


def mult_mechs(sorted_sens_coeffs, sorted_rxn_names, target_spcs, temps,
               ref_results):
    print('shape of ref_results: ', np.shape(ref_results))

    nmechs = np.shape(sorted_sens_coeffs)[0]  # first dim is nmechs
    for mech_idx in range(nmechs):
        single_mech_sens_coeffs = sorted_sens_coeffs[mech_idx]
        single_mech_rxn_names = sorted_rxn_names[mech_idx]
        single_mech_ref_results = ref_results[mech_idx]
        # NOTE: will probably need to fix sorted_rxn_names to be specific to the
        # mechanism as well
        print('shape of single_mech_sens_coeffs: ', np.shape(single_mech_sens_coeffs))
        print('shape of single_mech_rxn_names: ', np.shape(single_mech_rxn_names))

        write_file(single_mech_sens_coeffs, single_mech_rxn_names, target_spcs,
                   temps, single_mech_ref_results)


def write_file(single_mech_sens_coeffs, sorted_rxn_names, target_spcs, temps,
               ref_results, filename='sens_results.xlsx'):
    """

    :param single_mech_sens_coeffs: sens coeffs for a single mechanism
    :type single_mech_sens_coeffs: np.ndarray of shape (nrxns, nconds, ntargs,
        ntimes) if time-resolved or (nrxns, nconds, ntargs)
    :param sorted_rxn_names: not sure...
    :param target_spcs:
    :param temps:
    :param ref_results:
    :param filename:
    :return:
    """

    print('Writing results to Excel...')

    # Get the array of sorted sens_coeffs
    # sorted_idxs = sim_sens.get_sorted_idxs(sens_coeffs)
    # sorted_sens_coeffs = sim_sens.get_sorted_sens_coeffs(sens_coeffs,
    #                                                      sorted_idxs)
    # sorted_rxn_names = sim_sens.get_sorted_rxn_names(sorted_idxs, rxn_names)

    with pd.ExcelWriter(filename) as writer:
        # Each species gets three sheets: sorted, unsorted, and ref_results
        # ONLY TWO! for now...
        for spc_idx, spc in enumerate(target_spcs):
            # # Unsorted sens
            # spc_sens_array = sens_coeffs[:, :, spc_idx]
            # sens_df = pd.DataFrame(spc_sens_array, index=rxn_names)
            # sheet_name = f'{spc}, unsorted'
            # sens_df.to_excel(writer, sheet_name=sheet_name, header=temps)

            # Sorted sens

            sorted_spc_sens_array = single_mech_sens_coeffs[:, :, spc_idx]
            sorted_spc_rxn_names = sorted_rxn_names[:, spc_idx]

            print('shape of sorted_spc_sens_array: ',
                  np.shape(sorted_spc_sens_array))
            print('shape of sorted_spc_rxn_names: ',
                  np.shape(sorted_spc_rxn_names))
            # Note:
            # sorted_sens_df = pd.DataFrame(
            #     sorted_spc_sens_array, index=sorted_spc_rxn_names)
            # sheet_name = f'{spc}, sorted'
            # sorted_sens_df.to_excel(writer, sheet_name=sheet_name, header=temps)
            sorted_spc_sens_array = np.transpose(sorted_spc_sens_array)
            sorted_sens_df = pd.DataFrame(
                sorted_spc_sens_array, index=temps,)

            sheet_name = f'{spc}, sorted'
            sorted_sens_df.to_excel(writer, sheet_name=sheet_name, header=sorted_spc_rxn_names )
            # Reference concentrations
            spc_ref_results = ref_results[:, spc_idx]
            print('shape of spc_ref_results: ',
                  np.shape(spc_ref_results))
            ref_results_df = pd.DataFrame(spc_ref_results, index=temps)
            sheet_name = f'{spc}, ref_results'
            ref_results_df.to_excel(writer, sheet_name=sheet_name)
