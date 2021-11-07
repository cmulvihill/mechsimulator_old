""" Provides common functions for parsing files
"""

import numpy as np


def chk_entry(entry):
    """ Converts a value to None if it is either '-' or NaN, the latter of which
        means it was blank in the Excel file

        :param entry: entry from a Pandas dataframe that was read from Excel
        :type entry: either a str or a float
        :return chkd_entry: either the exact same as entry or None
        :rtype: either the same type as entry or None if it was '-' or NaN
    """

    if isinstance(entry, str) and entry == '-':
        chkd_entry = None
    elif isinstance(entry, float) and np.isnan(entry):
        chkd_entry = None
    else:
        chkd_entry = entry

    return chkd_entry
