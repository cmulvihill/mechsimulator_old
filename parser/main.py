from mechsimulator.parser import exp
from mechsimulator.parser import spc
from mechsimulator.parser import mech

ALLOWED_FILETYPES = ('exp', 'csv', 'cti')  # could add correlation for opt.


def mult_files(filenames, filetype, opts=None):
    """ Parses multiple files of the same type

        :param filenames: filenames
        :type: list
        :param filetype: file type; should be 'exp', 'spc', or 'mech'
        :type filetype: str
        :param opts: options for parsing
        :type opts: dict
        :return objs: a list of objects of the type indicated by filetype
        :rtype: list
    """

    objs = []
    for filename in filenames:
        obj = single_file(filename, filetype, opts=opts)
        objs.append(obj)

    return objs


def single_file(filename, filetype, opts=None):
    """ Loads a single object given a single filename

        :param filename: filename
        :type filename: str
        :param filetype: file type
        :type filetype: str
        :param opts: options for parsing
        :type opts: dict
        :return obj: object created from parsing the file
    """

    opts = opts or {}  # set to empty dict if None

    if filetype == 'exp':
        obj = exp.load_exp_set(filename)
    elif filetype == 'spc':
        quotechar = opts.get('quotechar')
        obj = spc.load_mech_spc_dct(filename, quotechar=quotechar)
    elif filetype == 'mech':
        obj = mech.load_solution_obj(filename)
    else:
        raise NotImplementedError(f"filetype '{filetype}' is not allowed")

    return obj
