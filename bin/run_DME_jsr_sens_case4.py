import cantera as ct
import numpy as np
from mechsimulator.sim import sens as sim
from mechsimulator.plotter import sens as plotter
from mechsimulator.writer import sens as writer

# Inputs

target_spcs = ['CH3OCH3',
               'CO',
               'CO2',
               'CH2O',
               'O2',
               'H2O',
               'CH3OCHO',
]

# This is for case 4
mech_filename = 'data/HPmechv3.cti'
temps = np.linspace(905, 495, 83)  # K
pressure = 100  # atm
res_time = None  # s
mdot = 7.86609e-5  # kg/s
vol = 0.4e-6  # m^3
# mix = {'CH3OCH3': 0.00644, 'O2': 0.1109, 'N2': 0.88266}  # real values

# To get the same fuel concentration as the stoich case but keeping the same phi
mix = {'CH3OCH3': 0.01203, 'O2': 0.2072, 'N2': 0.78077}
factor = 0.01
num_top_rxns = 35  # number of rxns to plot
file_prefix = 'case4_fake'

# Run calculations
mech = ct.Solution(mech_filename)
sens_coeffs, ref_results = sim.jsr(
    temps, pressure, mix, mech, target_spcs, res_time, vol, mdot=mdot, factor=factor)

# Get the unsorted rxn names
rxn_names = sim.get_rxn_names(mech)

# Sort the results
sorted_idxs = sim.get_sorted_idxs(sens_coeffs)

# Plot the results
plot_filename = file_prefix + '.pdf'
plotter.plot_sens(sens_coeffs, num_top_rxns, rxn_names, target_spcs, temps,
                  ref_results, filename=plot_filename)

# Write sensitivity results to file
excel_filename = '../bin/' + file_prefix + '.xlsx'
writer.write_file(sens_coeffs, target_spcs, temps, rxn_names, ref_results,
                   filename=excel_filename)
