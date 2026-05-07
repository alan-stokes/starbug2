import os
from starbug2.utils import p_error


def usage(docstring, verbose=0):
    """
    outputs the usage.
    :param docstring: the doc string to output
    :param verbose: if to do so in verbose mode
    :return: 1 when complete.
    """
    if verbose:
        p_error(docstring)
    else:
        p_error("%s\n" % docstring.split('\n')[1])
    return 1

def parse_cmd(args):
    """
    parses an args command.
    :param args: the args array.
    :return: tuple of the command and the rest of the args array.
    :rtype: (str, array[str])
    """
    cmd = os.path.basename(args[0])
    return cmd, args[1:]

