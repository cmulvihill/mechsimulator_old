import os
import matplotlib.backends.backend_pdf as plt_pdf
from mechsimulator.simulator import util as sim_util
from mechsimulator.parser.exp import ALLOWED_UNITS


def build_pdf(figs_axes, filename='output.pdf', path=None):
    """ Produce a PDF with one reaction per page

        :param figs_axes: list of pairs of figure and axis objects
        :type figs_axes: list [(fig1, axis1), ...]
        :param filename: filename for the output pdf; default is 'output.pdf'
        :type filename: str
        :param path: path for the output pdf; default is the current directory
        :type path: str
    """

    print('Producing PDF...')
    if path is not None:
        filename = os.path.join(path, filename)
    pdf = plt_pdf.PdfPages(filename)
    for fig, _ in figs_axes:  # don't need the axes
        pdf.savefig(fig)
    pdf.close()


def get_target_titles(exp_set, yunit=None):

    meas_type = exp_set['overall']['meas_type']
    if meas_type == 'abs':
        target_titles = exp_set['overall']['wavelength']
        yunit = yunit or '%'
        ylabel = f'Fractional absorption ({yunit})'
        yquant = 'abs'
    elif meas_type == 'emis':
        target_titles = exp_set['overall']['wavelength']
        ylabel = 'Normalized emission'
        yquant = None
    elif meas_type == 'ion':
        raise NotImplementedError("The 'ion' measurement type is not working.")
    elif meas_type in ('conc', 'outlet'):
        target_titles = list(exp_set['spc'].keys())
        yunit = yunit or ''
        ylabel = f'Mole fraction ({yunit})'
        yquant = 'conc'
    elif meas_type == 'pressure':
        target_titles = None
        yunit = yunit or 'atm'
        ylabel = f'Pressure ({yunit})'
        yquant = 'pressure'
    elif meas_type == 'idt':
        idt_targets = exp_set['plot']['idt_target']
        idt_methods = exp_set['plot']['idt_method']
        target_titles = []
        for idt_method in idt_methods:
            for idt_target in idt_targets:
                target_title = idt_target + ', ' + idt_method
                target_titles.append(target_title)
        yunit = yunit or 's'
        ylabel = f'Ignition delay time ({yunit})'
        yquant = 'time'

    return target_titles, ylabel, yquant


def get_cond_titles(exp_set, conds_source, xunit=None):

    sim_util.check_sources('plot', conds_source)  # dummy value for x_source

    # Load some data
    meas_type = exp_set['overall']['meas_type']
    plot_var = exp_set['plot']['variable']
    units = ALLOWED_UNITS[plot_var][0][0]

    # If on a time-resolved measurement
    if meas_type in ('abs', 'emis', 'conc', 'pressure'):
        cond_titles = []
        # If indicated, use the 'plot' conditions instead of the experiments
        if conds_source == 'plot':
            # Loop over each value in the 'plot'-defined values
            conds = sim_util.get_plot_conds(exp_set)
            for cond in conds:
                cond_title = f'{cond} {units}'
                cond_titles.append(cond_title)
        # Otherwise, use the experimental conditions
        elif conds_source == 'exps':
            # Loop over each experimental condition
            exp_objs = exp_set['exp_objs']
            for exp_obj in exp_objs:
                cond = exp_obj['conds'][plot_var][0]  # [0] omits uncertainties
                cond_title = f'{cond} {units}'
                cond_titles.append(cond_title)
        if xunit is not None:
            xlabel = f'Time ({xunit})'
        else:
            xlabel = 'Time (s)'
        xquant = 'time'

    # If on a non-time-resolved measurement
    elif meas_type in ('idt', 'outlet'):
        if meas_type == 'idt':
            cond_titles = ['Ignition delays']
        if meas_type == 'outlet':
            cond_titles = ['Outlet concentrations']
        if xunit is None:  # if no xunit given...
            xunit = ALLOWED_UNITS[plot_var][0][0]  # ...use the default unit
        xlabel = ALLOWED_UNITS[plot_var][2] + f' ({xunit})'
        # If on ignition delay time and using temperature, override the label
        if meas_type == 'idt' and plot_var == 'temp':
            xlabel = '1000/Temperature (K^-1)'
        xquant = plot_var

    # Haven't figured this one out yet
    elif meas_type == 'ion':
        raise NotImplementedError("The 'ion' measurement type is not working.")

    return cond_titles, xlabel, xquant
