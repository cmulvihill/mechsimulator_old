import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mechsimulator.plotter import util
from mechsimulator.sim import run_set
from mechsimulator.parser.exp import ALLOWED_UNITS
from mechsimulator.sim import util as sim_util

# Define some stuff for plot formatting
COLORS = ['Red', 'Blue', 'Green', 'Black']
LINES = ['-', '--', '-.', ':']
REAC = {'st': 'ST',
        'rcm': 'RCM',
        'jsr': 'JSR',
        'pfr': 'PFR',
        'const_t_p': 'Const. TP'}
MEAS = {'abs': 'Absorption',
        'emis': 'Emission',
        'conc': 'Concentration',
        'ion': 'Ion',
        'pressure': 'Pressure',
        'idt': 'IDT',
        'outlet': 'Outlet'}
DEFAULT_OPTIONS = {
    'group_by':     'target',
    'ignore_exps':  False,
    'rows_cols':    (4, 3),
}


def single_set(set_result, uniform_times, exp_set, mech_names=None,
               options=None):

    options = options or DEFAULT_OPTIONS  # sets to default if None
    num_mechs = len(set_result)
    mech_names = [f'mech {idx+1}' for idx in range(num_mechs)] or mech_names
    allowed_sort_methods = ('cond', 'target')
    group_by = options['group_by']
    assert group_by in allowed_sort_methods, (
        f"'group_by' should be in {allowed_sort_methods}, not {group_by}")

    # Build the empty figures and axes
    figs_axes = build_figs_and_axes(exp_set, mech_names, options)

    # Loop over each mechanism and plot
    for mech_idx in range(num_mechs):
        mech_result = set_result[mech_idx]
        figs_axes = single_mech(mech_result, uniform_times, exp_set, figs_axes,
                                options)

    # Add experimental data if present
    exp_data = get_exp_data(exp_set)


    return figs_axes


def single_mech(mech_result, uniform_times, exp_set, figs_axes, options):

    meas_type = exp_set['overall']['meas_type']
    if meas_type in ('abs', 'emis'):
        pass
    elif meas_type == 'idt':
        pass
    elif meas_type == 'conc':
        figs_axes = conc(mech_result, uniform_times, figs_axes, options)

    return figs_axes


def conc(mech_result, uniform_times, figs_axes, options):

    # Get some initial data
    nrows, ncols = options['rows_cols']
    plts_per_pg = nrows * ncols
    assert np.ndim(mech_result) == 3, (
        f"'mech_result' should have 3 dimensions, not {np.ndim(mech_result)}")
    group_by = options['group_by']
    num_conds, num_targets, _ = np.shape(mech_result)
    if group_by == 'cond':
        num_grps = num_conds  # number of plotting groups
        num_plts = num_targets  # number of plots per group
    else:  # 'target'
        num_grps = num_targets
        num_plts = num_conds
    pgs_per_grp = int(num_plts / plts_per_pg) + 1

    # Loop over each group of pages (i.e., either each cond or each target)
    for grp_idx in range(num_grps):
        # Loop over pages in each group
        for pg_idx in range(pgs_per_grp):
            fig_idx = grp_idx * pgs_per_grp + pg_idx
            fig, axs = figs_axes[fig_idx]
            # Loop over plots on each page
            for plt_idx_pg in range(plts_per_pg):
                plt_idx_grp = pg_idx * plts_per_pg + plt_idx_pg  # overall idx
                if plt_idx_grp < num_plts:
                    if group_by == 'cond':
                        y_data = mech_result[grp_idx, plt_idx_grp, :]
                    else:  # 'target'
                        y_data = mech_result[plt_idx_grp, grp_idx, :]
                    current_ax = axs[plt_idx_pg]
                    current_ax.plot(uniform_times, y_data)

    return figs_axes


def build_figs_and_axes(exp_set, mech_names, options):
    """

        :param exp_set:
        :param mech_names:
        :param group_by: either 'cond' or 'target'
        :param ignore_exps:
        :return:
    """

    def create_figs(grp_titles, plt_titles, xlabel, ylabel, mech_names):
        """

            :param pg_titles: page titles
            :type pg_titles: list
            :param plt_titles: plot titles
            :type plt_titles: list
            :param xlabel: units for the x axes
            :type xlabel: str
            :param ylabel: units for the y axes
            :type ylabel: str
            :param mech_names: nicknames for the mechanisms
            :type mech_names: list
            :return:
        """

        # Load/calculate various data
        nrows, ncols = options['rows_cols']
        plts_per_pg = nrows * ncols
        num_grps = len(grp_titles)
        plts_per_grp = len(plt_titles)
        pgs_per_grp = int(plts_per_grp / plts_per_pg) + 1  # int rounds down

        # Loop over each group and build the figures and axes
        figs_axes = []
        for grp_idx in range(num_grps):
            # Loop over each page within the group
            for pg_idx in range(pgs_per_grp):
                fig = plt.figure(figsize=(8.5, 11))
                axes = []
                # Loop over each plot on the page
                for plt_idx_pg in range(plts_per_pg):
                    plt_idx_grp = pg_idx * plts_per_pg + plt_idx_pg
                    if plt_idx_grp < plts_per_grp:
                        axis = fig.add_subplot(4, 3, plt_idx_pg + 1)
                        plt_title = plt_titles[plt_idx_grp]
                        plt.title(plt_title, fontsize=10, loc='center', y=0.97)
                        plt.xticks(fontsize=8)
                        plt.yticks(fontsize=8)
                        plt.subplots_adjust(left=None, bottom=None, right=None,
                                            top=None, wspace=0.38, hspace=0.23)
                        axis.yaxis.set_major_formatter(
                            FormatStrFormatter('%0.2E'))
                        axes.append(axis)
                    if plt_idx_pg == 0:  # if first plot of page, add page title
                        title = grp_titles[grp_idx]
                        if pgs_per_grp > 1:
                            title += f'\n(pg. {pg_idx + 1} of {pgs_per_grp})'
                        fig.suptitle(title, y=0.99, fontsize=16)

                add_headers(exp_set, mech_names)
                add_footers(xlabel, ylabel)
                figs_axes.append([fig, axes])

        return figs_axes

    # Load and check the sorting method
    group_by = options['group_by']

    # Get the titles and axis labels
    cond_titles, xlabel = util.get_cond_titles(exp_set, options)
    target_titles, ylabel = util.get_target_titles(exp_set)

    # Set the titles according to measurement type and grouping method
    if group_by == 'cond':
        pg_titles = cond_titles
        plt_titles = target_titles
    else:  # 'target'
        pg_titles = target_titles
        plt_titles = cond_titles

    # Create the figures and axes
    figs_axes = create_figs(pg_titles, plt_titles, xlabel, ylabel,
                            mech_names)

    return figs_axes


def add_headers(exp_set, mech_names):

    # Make some text describing the legends
    header_x_positions = [0.06, 0.37, 0.68]
    header_y_positions = [0.9, 0.92]
    for mech_idx, mech_name in enumerate(mech_names):
        header = f'{COLORS[mech_idx % 4]} lines: {mech_name}'
        if mech_idx < 3:
            y_idx = 0
        else:
            y_idx = 1
        plt.figtext(header_x_positions[mech_idx % 3], header_y_positions[y_idx],
                    header, fontsize=12, color=COLORS[mech_idx % 4])

    # Make some text describing the experimental set
    source = exp_set['overall']['source']
    description = exp_set['overall']['description']
    reac_type = exp_set['overall']['reac_type']
    meas_type = exp_set['overall']['meas_type']
    plt.figtext(0.01, 0.98, f'Source: {source}', fontsize=10)
    plt.figtext(0.01, 0.96, f'Description: {description}', fontsize=10)
    plt.figtext(0.77, 0.98, f'Reac. type: {REAC[reac_type]}', fontsize=10)
    plt.figtext(0.77, 0.96, f'Meas. type: {MEAS[meas_type]}', fontsize=10)


def add_footers(x_units, y_units):

    footnote1 = f'Y-axis units: {x_units}\n'
    footnote2 = f'X-axis units: {y_units}\n'
    footnotes = footnote1 + footnote2
    plt.figtext(0.11, 0.06, footnotes, fontsize=8, va="top", ha="left")


def get_exp_data(exp_set):
    """ Searches all experiment objects in an exp_set and gets the experimental
        data

    :param exp_set:
    :return:
    """

    # Initialize
    shape = sim_util.get_mech_result_shape(exp_set, 'exps')
    exp_data = np.ndarray(shape, dtype=float)  # float makes NaNs work properly
    meas_type = exp_set['overall']['meas_type']
    uniform_times = sim_util.get_uniform_times(exp_set, 'exps')

    # Loop over each exp_obj and get the data based on the meas_type
    exp_objs = exp_set['exp_objs']
    for cond_idx, exp_obj in enumerate(exp_objs):  # each obj is a condition
        result = exp_obj['result']
        if meas_type == 'conc':
            exp_time = exp_obj['result']['time']
            # Get all the species concentrations
            for spc_idx, spc in enumerate(exp_set['spc'].keys()):
                if spc in result:
                    raw_ydata = result[spc]
                    ydata = sim_util.interp(raw_ydata, exp_time, uniform_times)
                    exp_data[cond_idx, spc_idx, :] = ydata
                else:  # set to NaN if the spcs isn't defined
                    exp_data[cond_idx, spc_idx, :] = np.nan

    return exp_data


def add_exp_data(figs_axes, exp_data, group_by):





