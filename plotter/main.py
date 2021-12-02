from mechsimulator.plotter import outcome
from mechsimulator.plotter import sens
from mechsimulator import simulator
import numpy as np

NUM_TOP_RXNS = 25  # used for sensitivity and ROP plotting


def mult_sets(exp_sets, gases, mech_spc_dcts, calc_types, x_sources,
              conds_sources, mech_opts_lst=None):
    """ Plots any number of sets (against any number of mechanisms)

        :param exp_sets:
        :param gases:
        :param mech_spc_dcts:
        :param calc_types:
        :param x_sources:
        :param conds_sources:
        :param mech_opts_lst:
        :return:
    """

    mult_sets_figs_axes = []
    for idx, exp_set in enumerate(exp_sets):
        set_figs_axes = single_set(
            exp_set, gases, mech_spc_dcts, calc_types[idx], x_sources[idx],
            conds_sources[idx], mech_opts_lst=mech_opts_lst)
        # Later, will add some kind of header page for each set...for now, just
        # add the lists of figures together
        mult_sets_figs_axes.extend(set_figs_axes)

    return mult_sets_figs_axes


def single_set(exp_set, gases, mech_spc_dcts, calc_type, x_source,
               conds_source, mech_opts_lst=None):
    """ Plots a single set (i.e., with any number of mechanisms)

        :param exp_set:
        :param gases:
        :param mech_spc_dcts:
        :param calc_type:
        :param x_source:
        :param conds_source:
        :param mech_opts_lst:
        :return:
    """

    if calc_type == 'outcome':
        # Calculate the simulation data
        set_ydata, set_xdata = simulator.main.single_set(
            exp_set, gases, mech_spc_dcts, 'outcome', x_source, conds_source,
            mech_opts_lst=mech_opts_lst)
        # Plot the data
        figs_axes = outcome.single_set(set_ydata, set_xdata, exp_set,
                                       conds_source)
    elif calc_type == 'sens':
        # Calculate the sensitivity coefficients
        set_sens, set_xdata = simulator.main.single_set(
            exp_set, gases, mech_spc_dcts, 'sens', x_source, conds_source,
            mech_opts_lst=mech_opts_lst)
        # Sort the sensitivity coefficients
        sorted_set_sens, sorted_set_rxns = simulator.sort_sens.sort_single_set(
            set_sens, gases)
        # Calculate the reference results
        set_ref_results, _ = simulator.main.single_set(
            exp_set, gases, mech_spc_dcts, 'outcome', x_source, conds_source,
            mech_opts_lst=mech_opts_lst)
        # Plot the data
        figs_axes = sens.single_set(sorted_set_sens[:, :NUM_TOP_RXNS],
                                    set_xdata,
                                    sorted_set_rxns[:, :NUM_TOP_RXNS],
                                    set_ref_results, exp_set, conds_source)
    else:
        raise NotImplementedError(f"calc_type {calc_type} not implemented!")

    return figs_axes
