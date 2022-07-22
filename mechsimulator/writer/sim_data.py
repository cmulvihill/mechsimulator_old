import numpy as np
import pandas as pd
from mechsimulator.plotter import util
# from mechsimulator.plotter import comparisons


def write_mech_results(exp_set, mech_results, mech_nicknames,
                       pathname='results'):
    """

    :param exp_set:
    :param mech_results: contains simulation data for any number of mechanisms
    :type mech_results: Numpy array of shape (num_mechs, num_temps, num_spcs)
        (this is only true for JSR; will be different for other experiments)
    :return:
    """
    num_mechs, num_temps, num_spcs = np.shape(mech_results)
    assert num_mechs == len(mech_nicknames)
    target_spcs = list(exp_set['spc'].keys())
    var, temps = util.get_plotting_variable(exp_set)

    for mech_idx, mech_nickname in enumerate(mech_nicknames):
        mech_result = mech_results[mech_idx, :, :]
        filename = mech_nickname + '.xlsx'
        with pd.ExcelWriter(pathname + '/' + filename) as writer:
            mech_result_df = pd.DataFrame(mech_result, index=temps)
            mech_result_df.to_excel(writer, header=target_spcs)


# def write_exp_data(exp_set, pathname='results'):
#
#     exp_data, temps = comparisons.get_exp_data_jsr(exp_set)
#     source = exp_set['overall']['source']
#     description = exp_set['overall']['description']
#     target_spcs = list(exp_set['spc'].keys())
#     filename = source + ', ' + description + '.xlsx'
#
#     # Each species gets its own sheet
#     with pd.ExcelWriter(pathname + '/' + filename) as writer:
#         for spc_idx, target_spc in enumerate(target_spcs):
#             spc_data = exp_data[spc_idx, :]  # all data for one spc
#             # If all entries for the species are empty, print a notice
#             if all(np.isnan(spc_data)):
#                 print(f'No experimental data found for spc {target_spc}')
#             # Otherwise, write to the file
#             else:
#                 spc_data = exp_data[spc_idx, :]
#                 spc_data_df = pd.DataFrame(spc_data, index=temps)
#                 sheet_name = target_spc
#                 spc_data_df.to_excel(writer, sheet_name=sheet_name)
