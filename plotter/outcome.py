import copy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter as Formatter
from mechsimulator.plotter import util
from mechsimulator.parser.exp import ALLOWED_UNITS

# Define some stuff for plot formatting
COLORS = ['Red', 'Blue', 'Green', 'Black']
LINESTYLES = ['-', '--', '-.', ':']
REAC = {
    'st': 'ST',
    'rcm': 'RCM',
    'jsr': 'JSR',
    'pfr': 'PFR',
    'const_t_p': 'Const. TP',
}
MEAS = {
    'abs': 'Absorption',
    'emis': 'Emission',
    'conc': 'Concentration',
    'ion': 'Ion',
    'pressure': 'Pressure',
    'idt': 'IDT',
    'outlet': 'Outlet',
}
DEFAULT_OPTS = {
    'plot_points':  False,
    'xunit':        None,
    'xlim':         None,
    'yunit':        None,
    'ylim':         None,
    'omit_targets': None,
    'rows_cols':    (4, 3),
    'marker_size':  15,
    'group_by':     'target',
    'exp_color':    'black',
    'xscale':       'linear',
    'yscale':       'linear',
}
POINT_MEAS_TYPES = ('outlet', 'idt',)


# NOTE TO SELF: might delete this function entirely at some point in favor of
# the external call from plotter.main
def mult_sets(set_ydata_lst, set_xdata_lst, exp_sets, conds_source_lst,
              mech_names=None):
    pass


def single_set(set_ydata, set_xdata, exp_set, conds_source, mech_names=None):
    """ Plots results from a single set (i.e., any number of mechanisms). Also
        plots the experimental data if there are any

        :param set_ydata: the simulation results for a single set (i.e., for
            multiple mechanisms); shape=(num_mechs, num_conds, num_targets,
            num_times) (the last dimension is omitted for some meas_types)
        :type set_ydata: numpy.ndarray
        :param set_xdata: the uniform x data used for all mechanisms; 1-D
        :type set_xdata: numpy.ndarray
        :param exp_set: exp_set object
        :type exp_set: dict
        :param conds_source: the source of the conditions; 'plot' or 'exps'
        :type conds_source: str
        :param mech_names: mechanism nicknames
        :type mech_names: list
        :return figs_axes: list of figures with plots [(fig1, ax1), ...]
        :rtype: list
    """

    # Initialize some variables
    opts = read_options(exp_set)
    num_mechs = len(set_ydata)
    mech_names = mech_names or [f'mech {idx + 1}' for idx in range(num_mechs)]

    # Build the empty figures and axes
    figs_axes = build_figs_and_axes(exp_set, mech_names, opts, conds_source)

    # Loop over each mechanism and plot
    for mech_idx in range(num_mechs):
        mech_ydata = set_ydata[mech_idx]
        mech_frmt = get_mech_frmt(exp_set, opts, mech_names, conds_source,
                                  mech_idx=mech_idx)
        # Note: xdata is same for all mechs, so using set_xdata
        figs_axes = single_mech(mech_ydata, set_xdata, figs_axes, mech_frmt)

    # Add experimental data if present
    exp_ydata = exp_set['overall']['exp_ydata']
    exp_xdata = exp_set['overall']['exp_xdata']
    mech_frmt = get_mech_frmt(exp_set, opts, mech_names, conds_source)
    figs_axes = single_mech(exp_ydata, exp_xdata, figs_axes, mech_frmt)

    return figs_axes


def single_mech(mech_ydata, mech_xdata, figs_axes, mech_frmt):
    """ Plots the results of a single mechanism on an existing set of figures

        :param mech_ydata: either the simulation results for a single mech or
            the experimental data
        :type mech_ydata: numpy.ndarray
        :param mech_xdata: the uniform x data used for all mechanisms; 1-D
        :type mech_xdata: numpy.ndarray
        :param figs_axes: figures to which plots are to be added
        :type figs_axes: list
        :param mech_frmt: formatting instructions for a single mech (or exps)
        :type mech_frmt: dict
        :return figs_axes: list of figures with plots; [(fig1, ax1), ...]
        :rtype: list
    """

    # Load/check some things
    ndims = np.ndim(mech_ydata)
    assert ndims in (2, 3), (
        f"'mech_ydata' should have 2 or 3 dimensions, not {ndims}")
    nrows, ncols = mech_frmt['rows_cols']
    group_by = mech_frmt['group_by']
    label = mech_frmt['label']
    color = mech_frmt['color']
    xlim = mech_frmt['xlim']
    ylim = mech_frmt['ylim']
    xconv = mech_frmt['xconv']
    yconv = mech_frmt['yconv']
    xscale = mech_frmt['xscale']
    yscale = mech_frmt['yscale']
    linestyle = mech_frmt['linestyle']
    marker = mech_frmt['marker']

    # Check for the case of inverse temperature
    if xconv == 'inv':
        mech_xdata = 1000 / mech_xdata
        xconv = 1

    # Get numbers of conditions and targets
    if ndims == 3:
        num_conds, num_targets, _ = np.shape(mech_ydata)
    else:  # ndims = 2
        num_conds, num_targets = np.shape(mech_ydata)

    # Get information on how the plots are to be organized
    if ndims == 3:
        if group_by == 'cond':
            num_grps = num_conds  # number of plotting groups
            num_plts = num_targets  # number of plots per group
        else:  # 'target'
            num_grps = num_targets  # flipped relative to above
            num_plts = num_conds
    else:  # ndims = 2
        num_grps = 1  # for ndims = 2, all targets are plotted as one group
        num_plts = num_targets
    plts_per_pg = nrows * ncols  # maximum number of plots per page
    pgs_per_grp = int(num_plts / plts_per_pg) + 1  # number of pages per group

    # Loop over each group of pages
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
                            ydata = mech_ydata[grp_idx, plt_idx_grp]
                        else:  # 'target'
                            ydata = mech_ydata[plt_idx_grp, grp_idx]
                    else:  # ndims = 2
                        ydata = mech_ydata[:, plt_idx_pg]
                    current_ax = axs[plt_idx_pg]
                    current_ax.plot(mech_xdata * xconv, ydata * yconv,
                                    label=label, color=color,
                                    linestyle=linestyle, marker=marker)
                    if xlim is not None:
                        current_ax.set_xlim(xlim)
                    if ylim is not None:
                        current_ax.set_ylim(ylim)
                    current_ax.set_xscale(xscale)
                    current_ax.set_yscale(yscale)

    return figs_axes


def build_figs_and_axes(exp_set, mech_names, opts, conds_source):
    """ Initializes a set of figures with empty axes

        :param exp_set: exp_set object
        :type exp_set: dict
        :param mech_names: mechanism nicknames
        :type mech_names: list
        :param opts: plot options
        :type opts: dict
        :param conds_source: the source of the conditions; 'plot' or 'exps'
        :type conds_source: str
        :return figs_axes: list of figures with plots [(fig1, ax1), ...]
        :rtype: list
    """

    def create_figs(grp_titles, plt_titles, xlabel, ylabel, mech_names):
        """ Create the figures/axes given the group and plot titles

            :param grp_titles: titles for each group of plots
            :type grp_titles: list
            :param plt_titles: titles for the plots within each group
            :type plt_titles: list
            :param xlabel: label for the x axis
            :type xlabel: str
            :param ylabel: label for the y axis
            :type ylabel: str
        """

        # Load/calculate various data
        nrows, ncols = opts['rows_cols']
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
                        plt.subplots_adjust(wspace=0.38, hspace=0.23)
                        axis.yaxis.set_major_formatter(Formatter('%0.2E'))
                        axes.append(axis)
                    if plt_idx_pg == 0:  # if first plot of page, add page title
                        title = grp_titles[grp_idx]
                        title += f'\n(pg. {pg_idx + 1} of {pgs_per_grp})'
                        fig.suptitle(title, y=0.99, fontsize=16)

                add_headers(exp_set, mech_names)
                add_footers(xlabel, ylabel)
                figs_axes.append([fig, axes])

        return figs_axes

    # Load some information
    group_by = opts['group_by']
    xunit = opts['xunit']
    yunit = opts['yunit']

    # Get the titles and axis labels
    cond_titles, xlabel, _ = util.get_cond_titles(exp_set, conds_source,
                                                  xunit=xunit)
    target_titles, ylabel, _ = util.get_target_titles(exp_set, yunit=yunit)

    # Set the titles according to the grouping method
    if group_by == 'cond':
        grp_titles = cond_titles
        plt_titles = target_titles
    else:  # 'target'
        grp_titles = target_titles
        plt_titles = cond_titles

    # Create the figures and axes
    figs_axes = create_figs(grp_titles, plt_titles, xlabel, ylabel, mech_names)

    return figs_axes


def add_headers(exp_set, mech_names):
    """ Adds header text to a figure
    """

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


def add_footers(xlabel, ylabel):
    """ Adds footer text to a figure
    """

    footnote1 = f'Y-axis: {ylabel}\n'
    footnote2 = f'X-axis: {xlabel}\n'
    footnotes = footnote1 + footnote2
    plt.figtext(0.11, 0.06, footnotes, fontsize=10, va="top", ha="left")


def read_options(exp_set):
    """ Reads the plot options (i.e., plot_format) input from an exp_set; adds
        default values that are missing

        :param exp_set: exp_set object
        :type exp_set: dict
        :return opts: plot options dictionary
        :rtype: dict
    """

    meas_type = exp_set['overall']['meas_type']
    raw_opts = exp_set['plot_format']

    # If options doesn't exist, just use the default dictionary
    if raw_opts is {}:
        opts = copy.deepcopy(DEFAULT_OPTS)
    # If options does exist, fill in any missing values with defaults
    else:
        opts = copy.deepcopy(raw_opts)
        for parameter, value in DEFAULT_OPTS.items():
            if parameter not in opts:
                if meas_type == 'idt' and parameter == 'yscale':
                    opts[parameter] = 'log'  # the default for idt
                else:
                    opts[parameter] = value

    # If the meas_type should be plotted as points, set that
    if meas_type in POINT_MEAS_TYPES:
        opts['plot_points'] = True
        opts['group_by'] = 'cond'

    return opts


def get_mech_frmt(exp_set, opts, mech_names, conds_source, mech_idx=None):
    """ Gets plot formatting info for a single mech (or the experimental data)

        :param exp_set: exp_set object
        :type exp_set: dict
        :param opts: plot options
        :type opts: dict
        :param mech_names: mechanism nicknames
        :type mech_names: list
        :param conds_source: the source of the conditions; 'plot' or 'exps'
        :type conds_source: str
        :param mech_idx: index of the current mechanism, or None if experiment
        :type mech_idx: int
        :return mech_frmt: formatting instructions for a single mech (or exps)
        :rtype: dict
    """

    def get_conv_factors(units, quant):
        """ Get the factor to convert a quantity from the internal units to the
            desired units. Note: Celsius and Fahrenheit are not allowed (this is
            enforced in the parser)
        """

        if units is not None:
            allowed_units = ALLOWED_UNITS[quant][0]
            conv_factors = ALLOWED_UNITS[quant][1]
            assert units in allowed_units, (
                f"'{units}' are not allowed units for the quantity '{quant}'")
            idx = allowed_units.index(units)
            conv_factor = 1 / conv_factors[idx]  # conv_factor is the inverse
        else:
            conv_factor = 1

        return conv_factor

    # Initialize the dict with some general information
    mech_frmt = {'rows_cols': opts['rows_cols'],
                 'group_by': opts['group_by'],
                 'xlim': opts['xlim'],
                 'ylim': opts['ylim'],
                 'xscale': opts['xscale'],
                 'yscale': opts['yscale']}

    # Get the conversion factors
    xunit = opts['xunit']
    yunit = opts['yunit']
    _, _, xquant = util.get_cond_titles(exp_set, conds_source)
    _, _, yquant = util.get_target_titles(exp_set)
    mech_frmt['xconv'] = get_conv_factors(xunit, xquant)
    mech_frmt['yconv'] = get_conv_factors(yunit, yquant)
    if exp_set['overall']['meas_type'] == 'idt':
        mech_frmt['xconv'] = 'inv'  # so that inverse temp can be plotted

    # Get information that depends on whether a simulation or an experiment
    if mech_idx is not None:  # if mech_idx exists, it's a simulation
        mech_frmt['label'] = mech_names[mech_idx]
        mech_frmt['color'] = COLORS[mech_idx]
        mech_frmt['linestyle'] = LINESTYLES[mech_idx]
        mech_frmt['marker'] = ''
    else:  # if mech_idx is None, it's an experiment
        mech_frmt['label'] = 'Experiment'
        mech_frmt['color'] = opts['exp_color']
        # If indicated in options, plot points
        if opts['plot_points']:
            mech_frmt['linestyle'] = ''
            mech_frmt['marker'] = '.'
        # Otherwise, plot lines
        else:
            mech_frmt['linestyle'] = '-'
            mech_frmt['marker'] = ''

    return mech_frmt
