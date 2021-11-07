from mechsimulator.plotter import outcome
from mechsimulator.plotter import util as plot_util
from mechsimulator.sim import run_set
from mechsimulator.sim import util as sim_util
from mechsimulator.parser import exp as exp_parser
from mechsimulator.parser import spc as spc_parser
import numpy as np


def test_st():
    """ Tests a shock tube simulation
    """

    exp_filenames = [
        'data/Mulvihill_2021_N2O.xlsx'
    ]
    mech_filenames = [
        'data/N2O_mechanism_Mulvihill.cti',
    ]
    spc_csv_filenames = [
        'data/Mulvihill_N2O_species.csv',
    ]
    calc_type = 'outcome'
    plot_or_exps = 'exps'

    # Load the various objects using the filenames
    exp_sets = exp_parser.load_exp_sets(exp_filenames)
    gases = sim_util.load_solution_objs(mech_filenames)
    mech_spc_dcts = spc_parser.mech_spc_dcts_from_filenames(spc_csv_filenames)

    # Run the simulation
    set_results, set_times = run_set.mult_sets(
        exp_sets, gases, mech_spc_dcts, calc_type, plot_or_exps)
    print(set_results[0][-1][0])
    dim_check = set_results[0]
    print('shape of mech_result: ', np.shape(dim_check))
    mech_nicknames = ['my_mech']
    # Plot the results...for now, just run for a single set

    figs_axes = outcome.single_set(set_results[0], set_times[0], exp_sets[0],
                                   mech_nicknames)


    plot_util.build_pdf(figs_axes, filename='test.pdf')




if __name__ == '__main__':
    # test_jsr()
    test_st()
