from os.path import join
from mechsimulator import plotter

# ------------------------------------------------------------------------------
# Path and filename inputs

EXP_PATH = '../lib/exp_data'
MECH_PATH = '../lib/mechanisms'
SPC_PATH = '../lib/mechanisms'
OUTPUT_FILENAME = 'outcome.pdf'  # filename for the generated PDF


# ------------------------------------------------------------------------------
# Experimental inputs (should all be the same length)

EXP_FILENAMES = [
    'moshammer_2016_dme.xlsx',
]
# Desired calculation types (should be 'outcome', 'sens', or 'pathways')
CALC_TYPES = [
    'outcome',
]
# Sources of x data; only used for experiments with x data (e.g., time, height
# above burner) (should be 'plot' or 'exps'). If unsure, just use 'plot'
X_SRCS = [
    'plot',
]
# Sources of experimental conditions (should be 'plot' or 'exps')
COND_SRCS = [
    'plot',
]


# ------------------------------------------------------------------------------
# Mechanism inputs (should all be the same length)

MECH_FILENAMES = [
    'dme_couch_v1.cti',
    'dme_couch_v2.cti',
    'dme_couch_v3.cti',
]
SPC_CSV_FILENAMES = [
    'dme_couch.csv',
    'dme_couch.csv',
    'dme_couch.csv',
]
MECH_NAMES = [
    'v1',
    'v2',
    'v3',
]


# ------------------------------------------------------------------------------
# DON'T CHANGE ANYTHING BELOW THIS LINE


def add_paths(filenames, path):
    """ Takes a list of filenames and adds the path to each filename
    """

    # If a list of integers, just return the integers
    if all([isinstance(filename, int) for filename in filenames]):
        new_filenames = filenames
    # Otherwise, add the path to each filename
    else:
        new_filenames = []
        for filename in filenames:
            new_filenames.append(join(path, filename))

    return new_filenames


# Add paths to filenames
EXP_FILENAMES = add_paths(EXP_FILENAMES, EXP_PATH)
MECH_FILENAMES = add_paths(MECH_FILENAMES, MECH_PATH)
SPC_CSV_FILENAMES = add_paths(SPC_CSV_FILENAMES, SPC_PATH)

# Plot the results
FIGS_AXES = plotter.main.mult_sets_filenames(
    EXP_FILENAMES, MECH_FILENAMES, SPC_CSV_FILENAMES, CALC_TYPES,
    X_SRCS, COND_SRCS, mech_names=MECH_NAMES)
plotter.util.build_pdf(FIGS_AXES, filename=OUTPUT_FILENAME)
