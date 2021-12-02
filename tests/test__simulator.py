from mechsimulator import simulator
from mechsimulator.parser import exp as exp_parser
from mechsimulator.parser import spc as spc_parser


def test_jsr():
    """ Tests a JSR simulation
    """
    exp_filename = 'data/Tran_JSR_mdot.xlsx'
    mech_filenames = [
        '../../dee/mechanisms/Danilack_thermal_stereo_v8.cti',
    ]

    spc_csv_filenames = [
        'data/Danilack_species.csv',
    ]
    calc_type = 'outcome'
    plot_or_exps = 'exps'

    simulator.main.single_set(exp_filename, mech_filenames, spc_csv_filenames,
                              calc_type, plot_or_exps)


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
    gases = simulator.util.load_solution_objs(mech_filenames)
    mech_spc_dcts = spc_parser.mech_spc_dcts_from_filenames(spc_csv_filenames)

    #
    set_results, set_times = simulator.main.mult_sets(
        exp_sets, gases, mech_spc_dcts, calc_type, plot_or_exps)
    print(set_results)


if __name__ == '__main__':
    # test_jsr()
    test_st()
