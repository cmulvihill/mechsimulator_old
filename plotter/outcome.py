import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from mechsimulator.plotter import util
from mechsimulator.sim import run_set
from mechsimulator.parser.exp import ALLOWED_UNITS
from mechsimulator.sim import util as sim_util

# Define some stuff for plot formatting
COLORS = ['Red', 'Blue', 'Green', 'Black']
LINESTYLES = ['-', '--', '-.', ':']
REAC = {
    'st':           'ST',
    'rcm':          'RCM',
    'jsr':          'JSR',
    'pfr':          'PFR',
    'const_t_p':    'Const. TP',
}
MEAS = {
    'abs':          'Absorption',
    'emis':         'Emission',
    'conc':         'Concentration',
    'ion':          'Ion',
    'pressure':     'Pressure',
    'idt':          'IDT',
    'outlet':       'Outlet',
}
DEFAULT_OPTIONS = {
    'plot_points':  False,
    'ignore_exps':  False,
    'xunit':        None,
    'xlim':         None,
    'yunit':        None,
    'ylim':         None,
    'omit_targets': None,
    'group_by':     'cond',
    'rows_cols':    (4, 3),
    'marker_size':  12,
    'exp_color':    'gray',
}
POINT_MEAS_TYPES = ('outlet', 'idt',)  # could replace with time-resolved??


def single_set(set_ydata, set_xdata, exp_set, mech_names=None,
               options=None):
    """ Plots results from a single set (i.e., any number of mechanisms). Also
        plots the experimental data if there are any

        :param set_result:
        :param set_xdata:
        :param exp_set:
        :param mech_names:
        :param options:
        :return:
    """

    # Initialize some variables
    options = update_options(options)
    num_mechs = len(set_ydata)
    mech_names = [f'mech {idx+1}' for idx in range(num_mechs)] or mech_names
    allowed_sort_methods = ('cond', 'target')
    group_by = options['group_by']
    assert group_by in allowed_sort_methods, (
        f"'group_by' should be in {allowed_sort_methods}, not {group_by}")
    meas_type = exp_set['overall']['meas_type']

    # Build the empty figures and axes
    figs_axes = build_figs_and_axes(exp_set, mech_names, options)

    # Loop over each mechanism and plot
    for mech_idx in range(num_mechs):
        mech_ydata = set_ydata[mech_idx]
        mech_format = get_mech_format(meas_type, options, mech_names,
                                    mech_idx=mech_idx)
        figs_axes = single_mech(mech_ydata, set_xdata, figs_axes, options,
                                mech_format)  # xdata is same for all mechs

    # Add experimental data if present
    exp_ydata = exp_set['overall']['exp_ydata']
    exp_xdata = exp_set['overall']['exp_xdata']
    mech_format = get_mech_format(meas_type, options, mech_names, mech_idx=None)
    figs_axes = single_mech(exp_ydata, exp_xdata, figs_axes, options, 
                            mech_format)

    return figs_axes


def single_mech(mech_ydata, mech_xdata, figs_axes, options, mech_format):
    """ PLots the results of a single mechanism on an existing set of figures

        :param mech_ydata: either the simulation results for a single mech or
            the experimental data
        :param mech_xdata: the uniform time array used as the x axis if the
            measurement type is time-resolved
        :param exp_set:
        :param figs_axes:
        :param options:
        :return:
    """

    # Load/check some things
    nrows, ncols = options['rows_cols']
    ndims = np.ndim(mech_ydata)
    assert ndims in (2, 3), (
        f"'mech_ydata' should have 2 or 3 dimensions, not {ndims}")
    group_by = options['group_by']
    label = mech_format['label']
    color = mech_format['color']
    linestyle = mech_format['linestyle']
    marker = mech_format['marker']

    # Get numbers of conditions and targets
    if ndims == 3:
        num_conds, num_targets, _ = np.shape(mech_ydata)
    elif ndims == 2:
        num_conds, num_targets = np.shape(mech_ydata)

    # Get information on how the plots are to be organized
    if ndims == 3:
        if group_by == 'cond':
            num_grps = num_conds  # number of plotting groups
            num_plts = num_targets  # number of plots per group
        else:  # 'target'
            num_grps = num_targets  # flipped relative to above
            num_plts = num_conds
    elif ndims == 2:
        num_grps = 1  # for ndims = 2, all targets are plotted as one group
        num_plts = num_targets

    plts_per_pg = nrows * ncols  # number of plots per page
    pgs_per_grp = int(num_plts / plts_per_pg) + 1  # number of pages per group

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
                    if ndims == 3:
                        if group_by == 'cond':
                            y_data = mech_ydata[grp_idx, plt_idx_grp]
                        else:  # 'target'
                            y_data = mech_ydata[plt_idx_grp, grp_idx]
                    elif ndims == 2:
                        y_data = mech_ydata[:, plt_idx_pg]
                    current_ax = axs[plt_idx_pg]
                    current_ax.plot(mech_xdata, y_data, label=label,
                                    color=color, linestyle=linestyle,
                                    marker=marker)

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
                        axis = fig.add_subplot(nrows, ncols, plt_idx_pg + 1)
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


def update_options(options):
    """ Fills in missing entries in the options dictionary with default values

        :param options: initial options dictionary
        :type options: dict
        :return updated_options: updated options dictionary
        :rtype: dict
    """

    # If options doesn't exist, just use the default dictionary
    if options is None:
        updated_options = copy.deepcopy(DEFAULT_OPTIONS)
    # If options does exist, fill in any missing values with defaults
    else:
        updated_options = copy.deepcopy(options)
        for parameter, value in DEFAULT_OPTIONS.items():
            if parameter not in updated_options:
                updated_options[parameter] = value

    return updated_options


def get_mech_format(meas_type, options, mech_names, mech_idx=None):
    """ Gets plot formatting info for a single mech (or the experimental data)

        :param meas_type:
        :param options:
        :param mech_names:
        :param mech_idx:
        :return:
    """

    mech_format = {}
    # If on a simulation
    if mech_idx is not None:
        mech_format['label'] = mech_names[mech_idx]
        mech_format['color'] = COLORS[mech_idx]
        mech_format['linestyle'] = LINESTYLES[mech_idx]
        mech_format['marker'] = ''

    # Otherwise, if on an experiment (i.e., if mech_idx is None)
    else:
        mech_format['label'] = 'Experiment'
        mech_format['color'] = options['exp_color']
        # If indicated in options or by the type of experiment, plot points
        if meas_type in POINT_MEAS_TYPES or options.get('plot_points'):
            mech_format['linestyle'] = ''
            mech_format['marker'] = '.'
        # Otherwise, plot lines
        else:
            mech_format['linestyle'] = '-'
            mech_format['marker'] = ''

    return mech_format
