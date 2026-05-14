"""StarbugII Plotting Scripts
usage: starbug2-plot [-vhX] [-I CN000] [-o outfile] images.fits
    -h  --help           : show help screen
    -o  --output   fname : output filename
    -v  --verbose        : verbose mode

    -I  --inspect  CN000 : inspect a source in an array of images
    -X  --test           : plot a test image

        --style    fname : load a custom pyplot style sheet
        --dark           : plot in dark mode
"""
import os, sys, getopt

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

import starbug2.bin as scr
import starbug2
from starbug2.constants import (
    EXIT_EARLY, EXIT_FAIL, EXIT_SUCCESS, CAT_NUM, FILTER, EXT, IMAGE,
    BIN_TABLE, OUTPUT)
from starbug2.plot import load_style, plot_test, plot_inspect_source
from starbug2.utils import p_error, warn

VERBOSE = 0x01
SHOW_HELP = 0x02
STOP_PROC = 0x04
KILL_PROC = 0x08
DARK_MODE = 0x10

PTEST = 0x1000
PINSPECT = 0x2000


def plot_parse_argv(argv):
    options = 0
    set_opt = {}
    cmd, argv = scr.parse_cmd(argv)
    opts, args = getopt.gnu_getopt(
        argv, "hvXI:d:o:",
        ["help", "verbose", "test", "inspect=", "output=", "style=", "dark"]
    )

    for opt, opt_arg in opts:
        match opt:
            case "-h" | "--help":
                options |= (SHOW_HELP | STOP_PROC)
            case "-v" | "--verbose":
                options |= VERBOSE
            case "-o" | "--output":
                set_opt[OUTPUT] = opt_arg
            case "-d" | "--apfile":
                set_opt["APFILE"] = opt_arg

            case "-I"|"--inspect":
                options |= PINSPECT
                set_opt["INSPECT"] = opt_arg
            case "-X" | "--test":
                options |= PTEST
            case "--style":
                set_opt["STYLESHEET"] = opt_arg
            case "--dark":
                options |= DARK_MODE

    return options, set_opt, args


def plot_one_time_runs(options, set_opt, args):
    """
    runs plot one time

    :param options: the plot options
    :param set_opt: the options
    :param args: args
    :return: end state
    """

    if options & SHOW_HELP:
        scr.usage(__doc__, verbose=options&VERBOSE)

        if options & PINSPECT:
            p_error(fn_pinspect.__doc__)

        return EXIT_EARLY

    if _file_name := set_opt.get("STYLESHEET"):
        load_style(_file_name)

    if options & DARK_MODE:
        load_style("%s/extras/dark.style" % starbug2.__path__[0])
    
    if options & STOP_PROC:
        return EXIT_EARLY
    if options & KILL_PROC:
        p_error("..killing process\n")
        return EXIT_FAIL

    return EXIT_SUCCESS

def fn_pinspect(options, set_opt, images=None, tables=None):
    """
    Plot at a source position cutouts in a range of images.
    This requires a source list to be loaded, a list of image
    file and the source catalogue number to be given. This will
    take the form::

        $~ starbug2-plot -I CN123 source list.fits image*.fits

    :param options: The starbug2.bin.plot options integer
    :type options: int
    :param set_opt: The starbug2.bin.plot set opt dictionary
    :type set_opt: dict
    :param images: The list of fits image HDUs to cut out from
    :type images: list [HDU]
    :param tables: The source list to pull the source from. Must have a column
        with the name "Catalogue_Number"
    :type tables: list of astropy.Table
    :return: The output figure
    :rtype: plt.figure
    """

    fig = None
    if (cn := set_opt.get("INSPECT")) and images and tables:
        if (CAT_NUM in tables[0].col_names
                and cn in tables[0][CAT_NUM]):
            i = np.where(tables[0][CAT_NUM] == cn)[0]
            fig = plot_inspect_source(tables[0][i], images)
    else: p_error(
        "Must include the source {}, "
        "a list of images and a source list \n".format(CAT_NUM))
    return fig

def plot_main(argv):
    """
    plot main

    :param argv: the arguments for the plot.
    :return: None
    """
    warn("Still in development\n\n")
    options, set_opt, args = plot_parse_argv(argv)
    load_style("%s/extras/starbug.style" % starbug2.__path__[0])

    if options or set_opt:
        if exit_code := plot_one_time_runs(options, set_opt, args):
            return exit_code

    images = []
    tables = []
    for arg in args:
        if _file_name := os.path.exists(arg):
            fp = fits.open(arg)

            # THIS IS A HACK
            _filter = fp[0].header.get(FILTER)
            hdu = None
            for hdu in fp:
                if hdu.header.get(EXT) == IMAGE:
                    images.append(hdu)
                    break
                if hdu.header.get(EXT) == BIN_TABLE:
                    tables.append(Table(hdu.data))
                    break
            hdu.header[FILTER] = _filter

    fig = None
    if options & PTEST:
        fig, ax = plt.subplots(1, figsize=(3,2.5))
        plot_test(ax)

    if options & PINSPECT:
        fig = fn_pinspect(options, set_opt, images=images, tables=tables)

    if fig is not None:
        fig.tight_layout()
        if output := set_opt.get(OUTPUT):
            fig.savefig(output, dpi=300)
        else:
            plt.show()

def plot_main_entry():
    """Command Line entry point"""
    return plot_main(sys.argv)
