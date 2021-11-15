import os
import numpy as np
import cantera as ct
import matplotlib.backends.backend_pdf as plt_pdf
from mechsimulator.sim import util as sim_util
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


def get_target_titles(exp_set):

    meas_type = exp_set['overall']['meas_type']
    if meas_type == 'abs':
        target_titles = exp_set['overall']['wavelength']
        ylabel = 'Fractional absorption (%)'
    elif meas_type == 'emis':
        target_titles = exp_set['overall']['wavelength']
        ylabel = 'Normalized emission'
    elif meas_type == 'ion':
        raise NotImplementedError("The 'ion' measurement type is not working.")
    elif meas_type in ('conc', 'outlet'):
        target_titles = list(exp_set['spc'].keys())
        ylabel = 'Mole fraction'
    elif meas_type == 'pressure':
        target_titles = None
        ylabel = 'Pressure (atm)'
    elif meas_type == 'idt':
        idt_targets = exp_set['plot']['idt_target']
        idt_methods = exp_set['plot']['idt_method']
        target_titles = []
        for idt_method in idt_methods:
            for idt_target in idt_targets:
                target_title = idt_target + ', ' + idt_method
                target_titles.append(target_title)
        ylabel = 'Ignition delay time (s)'

    return target_titles, ylabel


def get_cond_titles(exp_set, options):

    # Load some data
    ignore_exps = options['ignore_exps']
    meas_type = exp_set['overall']['meas_type']
    plot_var = sim_util.get_plot_variable(exp_set)
    units = ALLOWED_UNITS[plot_var][0][0]

    # Check that the ignore_exps option was correctly used
    if ignore_exps:
        allowed_ignore_types = ('pressure', 'conc')
        assert meas_type in allowed_ignore_types, (
        "For the 'ignore_exps' option to be used, the meas_type must be one"
        f" of the following: {allowed_ignore_types}, not {meas_type}.")

    # If on a time-resolved measurement
    if meas_type in ('abs', 'emis', 'conc', 'pressure'):
        cond_titles = []
        # If indicated, use the 'plot' conditions instead of the experiments
        if ignore_exps:
            # Loop over each value in the 'plot'-defined values
            conds = sim_util.get_plot_conds(exp_set)
            for cond in conds:
                cond_title = f'{cond} {units}'
                cond_titles.append(cond_title)
        # Otherwise, use the experimental conditions
        else:
            # Loop over each experimental condition
            exp_objs = exp_set['exp_objs']
            for exp_obj in exp_objs:
                cond = exp_obj['conds'][plot_var][0]  # [0] omits uncertainties
                cond_title = f'{cond} {units}'
                cond_titles.append(cond_title)
        xlabel = ALLOWED_UNITS['time'][2]

    # If on a non-time-resolved measurement
    elif meas_type in ('idt', 'outlet'):
        if meas_type == 'idt':
            cond_titles = ['Ignition delays']
        if meas_type == 'outlet':
            cond_titles = ['Outlet concentrations']
        xlabel = ALLOWED_UNITS[plot_var][2]

    # Haven't figured this one out yet
    elif meas_type == 'ion':
        raise NotImplementedError("The 'ion' measurement type is not working.")

    return cond_titles, xlabel
