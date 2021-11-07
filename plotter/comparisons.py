import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as plt_pdf
from matplotlib.ticker import FormatStrFormatter
from mechsimulator.plotter import util

COLORS = ['Red', 'Blue', 'Green', 'Black']  # for plot formatting
LINES = ['-', '--', '-.', ':']


def plot_jsr_experiments(sim_results, exp_sets):

    # num_mechs, num_sets, num_temps, num_spcs = np.shape(sim_results)

    for set_idx, set_result in enumerate(sim_results):

        # Run some function here to get the temp range for the set

        for mech_idx, mech_result in enumerate(set_result):

            'hi'


def plot_single_set_jsr(set_result, mech_nicknames, exp_set):
    """

        :param sim_results: sim results for a single set of experiments
        :param exp_set:
        :return:
    """
    assert len(mech_nicknames) == len(set_result), (
        f'# of mech nicknames should match length of first dim of set_results ')

    var, temps = util.get_plotting_variable(exp_set)
    target_spcs = list(exp_set['spc'].keys())
    figs_axes = build_figs_and_axes(target_spcs, exp_set, mech_nicknames)

    # Plot the simulated data
    for mech_idx, mech_nickname in enumerate(mech_nicknames):
        for spc_idx, spc in enumerate(target_spcs):
            spc_result = set_result[mech_idx, :, spc_idx]  # get results
            fig_idx = int(spc_idx / 12)  # current fig (i.e., page) index
            fig, axs = figs_axes[fig_idx]
            ax_idx = spc_idx % 12  # each fig (i.e., page) only has 12 axes
            axs[ax_idx].plot(temps, spc_result, label=mech_nickname,
                             color=COLORS[mech_idx], linestyle=LINES[mech_idx])

    # Plot the experimental data
    exp_data, temps = get_exp_data_jsr(exp_set)
    for spc_idx, spc in enumerate(target_spcs):
        spc_data = exp_data[spc_idx, :]  # get exp_data for this species
        if spc_data.size != 0:  # if the species data is not empty
            # temps = spc_data[:, 0]
            # concs = spc_data[:, 1]
            fig_idx = int(spc_idx / 12)  # current fig (i.e., page) index
            fig, axs = figs_axes[fig_idx]
            ax_idx = spc_idx % 12  # each fig (i.e., page) only has 12 axes
            axs[ax_idx].plot(temps, spc_data, 'ko', label='Experiment',
                             markersize=3)

    return figs_axes


def get_exp_data_jsr(exp_set):
    """ Extracts the experimental data from an experimental set

        Note: at some point, should change this to be generic for any exp_set.
        However, I think the distinction will be by meas_type and not by
        reac_type? Not sure about that...if true, I'll need to make some
        changes in other parts of the code also. I think this might be false? I
        think what I care about is what is being varied from experiment to
        experiment

        Ok, another update: it depends on both reac_type and meas_type. For
        example, for ST IDT, I would want the results as a function of temperature.
        However, for ST absorption, I would want the results for each specific experiment
        as a separate entry. Plus, for this latter case, I would have to allow
        the arrays to be of different lengths due to different timespans of the
        data

        :param exp_set:
        :type exp_set:
        :return exp_data: experimental data for each spc
        :rtype: Numpy array of shape (num_spcs, num_exps)
        :return temps: temperatures for each experiment
        :rtype: Numpy array of shape (num_exps,)
    """

    # Read some data
    exp_objs = exp_set['exp_objs']
    target_spcs = list(exp_set['spc'].keys())
    num_exps = exp_set['overall']['num_exps']
    num_spcs = len(target_spcs)

    # Loop over each spc
    exp_data = np.zeros([num_spcs, num_exps])
    temps = np.zeros([num_exps])
    for spc_idx, target_spc in enumerate(target_spcs):
        for exp_idx, exp_obj in enumerate(exp_objs):
            temps[exp_idx] = exp_obj['conds']['temp'][0]
            if target_spc in exp_obj['result'].keys():
                exp_datum = exp_obj['result'][target_spc][0]
                # Set to NaN if ndarray; occurs if spc is present in the Excel
                # but is blank
                if isinstance(exp_datum, np.ndarray):
                    exp_datum = np.nan
            # If the target_spc is not in the exp_obj, return an empty array
            else:
                exp_datum = np.nan
            exp_data[spc_idx, exp_idx] = exp_datum

    return exp_data, temps


def build_figs_and_axes(target_spcs, exp_set, mech_nicknames):

    source = exp_set['overall']['source']
    description = exp_set['overall']['description']
    num_pgs = int(len(target_spcs) / 12) + 1
    num_plts = len(target_spcs)

    figs_axes = []
    for fig_idx in range(num_pgs):
        fig = plt.figure(figsize=(8.5, 11))
        axes = []
        plt_idx_start = 12 * fig_idx

        for plt_idx in range(12):
            if (plt_idx + plt_idx_start) < num_plts:
                axis = fig.add_subplot(4, 3, plt_idx + 1)
                plt.title(target_spcs[plt_idx + plt_idx_start], fontsize=10,
                          loc='center', y=0.97)
                plt.xticks(fontsize=8)
                plt.yticks(fontsize=8)
                plt.subplots_adjust(left=None, bottom=None, right=None,
                                    top=None, wspace=0.38, hspace=0.23)
                axis.yaxis.set_major_formatter(FormatStrFormatter('%0.2E'))
                axes.append(axis)

            if plt_idx == 0:  # if first species of a new page, add title
                title = source + ', ' + description + \
                        f' (pg. {fig_idx + 1} of {num_pgs})'
                fig.suptitle(title, y=0.97, fontsize=20)

        add_headers(mech_nicknames)
        add_footers('Kelvin', 'mole fraction')
        figs_axes.append([fig, axes])

    return figs_axes


def add_headers(mech_nicknames):

    header_x_positions = [0.06, 0.37, 0.68]
    header_y_positions = [0.91, 0.93]
    for mech_idx, mech_nickname in enumerate(mech_nicknames):
        header = f'{COLORS[mech_idx % 4]} lines: {mech_nickname}'
        if mech_idx < 3:
            y_idx = 0
        else:
            y_idx = 1
        plt.figtext(header_x_positions[mech_idx % 3], header_y_positions[y_idx],
                    header, fontsize=12, color=COLORS[mech_idx % 4])


def add_footers(x_units, y_units):

    footnote1 = f'Y-axis units: {x_units}\n'
    footnote2 = f'X-axis units: {y_units}\n'
    footnotes = footnote1 + footnote2
    plt.figtext(0.11, 0.06, footnotes, fontsize=8, va="top", ha="left")


def build_pdf(figs_axes, filename=None, path=None):
    """ Produce a PDF with one figure object per page

        :param figs: list of figure objects
        :type figs: list [fig1, fig2, ...]
        :param filename: filename for the output pdf; default is 'rate_comparison.pdf'
        :type filename: str
        :param path: path for the output pdf; default is the current directory
        :type path: str
    """
    print('Producing PDF...')
    if filename is None:
        filename = 'experiments_vs_sims.pdf'
    if path is not None:
        filename = os.path.join(path, filename)
    pdf = plt_pdf.PdfPages(filename)
    for (fig, _) in figs_axes:
        pdf.savefig(fig)
    pdf.close()
