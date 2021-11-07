"""
Parses spc.csv files containing species information to obtain mech_spc_dcts
"""

import csv
import copy
from mechsimulator.parser import rdkit_

ALLOWED_COLUMN_NAMES = (
    'name',
    'smiles',
    'inchi',
    'inchikey',
    'mult',
    'charge',
    'exc_flag',
    'sens'
)


def mech_spc_dcts_from_filenames(filenames, quotechar="'"):
    """ Obtains a single mech_spc_dct given a spc.csv filename

        :param filenames: filenames of the spc.csv files to be read
        :type filenames: list [filename1, filename 2, ...]
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return mech_spc_dcts: mech_spc_dcts for all filenames given
        :rtype: list
    """

    mech_spc_dcts = []
    for filename in filenames:
        mech_spc_dct = mech_spc_dct_from_filename(filename, quotechar=quotechar)
        mech_spc_dcts.append(mech_spc_dct)

    return mech_spc_dcts


def mech_spc_dct_from_filename(filename, quotechar="'"):
    """ Obtains a single mech_spc_dct given a spc.csv filename

        :param filename: the filename of the spc.csv file to be read
        :type filename: str
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return mech_spc_dct: identifying information on species in a mech
        :rtype: dct {spc1: spc_dct1, spc2: ...}
    """

    with open(filename, 'r') as file_obj:
        file_str = file_obj.read()

    mech_spc_dct = mech_spc_dct_from_str(file_str, quotechar=quotechar)

    return mech_spc_dct


def mech_spc_dct_from_str(file_str, quotechar="'"):
    """ Obtains a single mech_spc_dct given a string parsed from a spc.csv file

        :param file_str: the string that was read directly from the .csv file
        :type file_str: str
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return mech_spc_dct: identifying information on species in a mech
        :rtype: dct {spc1: spc_dct1, spc2: ...}
    """

    # Check for incorrect quotechar usage
    if quotechar == '"':
        assert "'" not in file_str, (
            f'A quotechar input of double quotation marks was selected, but at '
            f'least one instance of a single quotation mark exists in the '
            f'csv file. Use only double quotation marks.')
    elif quotechar == "'":
        assert '"' not in file_str, (
            f'A quotechar input of single quotation marks was selected, but at '
            f'least one instance of a double quotation mark exists in the '
            f'csv file. Use only single quotation marks.')
    else:
        raise NotImplementedError(
            f'The quotechar {quotechar} is not allowed. Options are either '
            f'single or double quotation marks.')

    # Split into lines
    lines = file_str.split('\n')

    # Build the mech_spc_dct line by line
    mech_spc_dct = {}
    for idx, line in enumerate(lines):
        if idx == 0:
            col_headers, num_cols = parse_first_line(line, quotechar=quotechar)
        else:
            cols = parse_line(line, idx, num_cols, quotechar=quotechar)
            spc_name, spc_dct = get_spc_dct(cols, col_headers)

            # Check that the species name was not already defined
            assert spc_name not in mech_spc_dct.keys(), (
                f'The species name {spc_name} appears in the csv file more than'
                f' once! The second time is on line {idx}, {line}.')

            # Fill in the spc_dct and then add it to the mech_spc_dct
            spc_dct = fill_in_spc_dct(spc_dct)  # add some information
            mech_spc_dct[spc_name] = spc_dct

    return mech_spc_dct


def get_spc_dct(cols, col_headers):
    """ Reads the parsed contents of a single line and converts to a spc_dct

        :param cols: a list of the entries in the line
        :type cols: list
        :param col_headers: a list of the entries in the first line
        :type col_headers: list
        :return spc_name: the mechanism name for the species
        :rtype: str
        :return spc_dct: identifying information for a single species
        :rtype: dct
    """

    # Build the spc_dct for this species, column by column
    spc_dct = {}
    for idx, col in enumerate(cols):
        if col_headers[idx] == 'name':
            spc_name = col
        # If charge or mult or exc_flag, convert to int if not empty
        elif col_headers[idx] in ('charge', 'mult', 'exc_flag') and col != '':
            spc_dct[col_headers[idx]] = int(col)
        else:
            spc_dct[col_headers[idx]] = col

    return spc_name, spc_dct


def fill_in_spc_dct(spc_dct):
    """ Looks at various identifiers in a spc_dct and fills in certain values

        :param spc_dct: identifying information for a single species
        :type: dct
        :return filled_spc_dct: identifying information for a single species
            that may now be supplemented with some information that was missing
        :rtype: dct
    """

    filled_spc_dct = copy.deepcopy(spc_dct)

    # Add SMILES or InChI if missing
    if 'inchi' not in filled_spc_dct:
        filled_spc_dct['inchi'] = rdkit_.to_inchi(rdkit_.from_smiles(
            filled_spc_dct['smiles']))

    elif 'smiles' not in filled_spc_dct:
        filled_spc_dct['smiles'] = rdkit_.to_smiles(rdkit_.from_inchi(
            filled_spc_dct['inchi']))

    # Add charge and exc_flag if missing; assume 0 for both
    if 'charge' not in filled_spc_dct or filled_spc_dct['charge'] == '':
        filled_spc_dct['charge'] = 0
    if 'exc_flag' not in filled_spc_dct or filled_spc_dct['exc_flag'] == '':
        filled_spc_dct['exc_flag'] = 0

    return filled_spc_dct


def parse_first_line(first_line, quotechar="'"):
    """ Parse the first line in the spc.csv file

        :param first_line: a string with the contents of the first line in the
            csv file
        :type first_line: str
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return headers: a list of the entries in the first line
        :rtype: list
        :return num_cols: the number of columns indicated by the first line
        :rtype: int
    """

    # Remove the UTF encoding if present
    if '\ufeff' in first_line:
        first_line = first_line.replace('\ufeff', '')

    # Parse the line
    headers = next(csv.reader([first_line], quotechar=quotechar))
    num_cols = len(headers)

    # Check that all column names are allowed
    for col_header in headers:
        assert col_header in ALLOWED_COLUMN_NAMES, (
            f'The column header {col_header} is not allowed. Options are'
            f' {ALLOWED_COLUMN_NAMES}.')
    # Check that at least one of the required chemical identifiers was provided
    assert 'inchi' in headers or 'smiles' in headers, (
        "At least one of the following chemical identifiers must be included in"
        " the csv file header: 'inchi' or 'smiles'.")
    # Check that the multiplicity was provided
    assert 'mult' in headers, (
        "The multiplicity, 'mult', must be included in the csv file header.")
    # Check that the name was provided
    assert 'name' in headers, (
        "The species name, 'name', must be included in the csv file header.")

    return headers, num_cols


def parse_line(line, idx, num_cols, quotechar="'"):
    """ Parse a line in the spc.csv file (other than the first line)

        :param line: a string with the contents of a single line in the csv file
        :type line: str
        :param idx: the index of the line
        :type idx: int
        :param num_cols: the number of columns in the csv file
        :type num_cols: int
        :param quotechar: the quotechar used to optionally ignore commas
        :type quotechar: str
        :return cols: a list of the entries in the line
        :rtype: list
    """

    cols = next(csv.reader([line], quotechar=quotechar))
    assert len(cols) == num_cols, (
        f'Line {idx}, {line}, has {len(cols)} columns. The first line in the'
        f'csv file indicates {num_cols} columns are needed.')

    return cols
