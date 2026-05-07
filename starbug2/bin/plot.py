"""StarbugII Plotting Scripts
usage: starbug2-plot [-vhX] [-I CN000] [-o outfile] images.fits ..
    -h  --help           : show help screen
    -o  --output   fname : output filename
    -v  --verbose        : verbose mode

    -I  --inspect  CN000 : inspect a source in an array of images
    -X  --test           : plot a test image

        --style    fname : load a custom pyplot style sheet
        --dark           : plot in dark mode
"""
import os, sys, getopt
from collections.abc import set_iterator

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

import starbug2.bin as scr
import starbug2
from starbug2.constants import EXIT_EARLY, EXIT_FAIL, EXIT_SUCCESS
from starbug2.plot import load_style, plot_test, plot_inspectsource
from starbug2.utils import p_error, warn

VERBOSE =0x01
SHOWHELP=0x02
STOPPROC=0x04
KILLPROC=0x08

DARKMODE=0x10

PTEST=   0x1000
PINSPECT=0x2000


def plot_parse_argv(argv):
    options=0
    set_opt={}
    cmd, argv = scr.parse_cmd(argv)
    opts, args = getopt.gnu_getopt(
        argv, "hvXI:d:o:",
        ["help", "verbose", "test", "inspect=", "output=", "style=", "dark"]
    )

    for opt, optarg in opts:
        match(opt):
            case "-h"|"--help":     options|=(SHOWHELP|STOPPROC)
            case "-v"|"--verbose":  options|=VERBOSE
            case "-o"|"--output":   set_opt["OUTPUT"]=optarg
            case "-d"|"--apfile":   set_opt["APFILE"]=optarg

            case "-I"|"--inspect":
                options|=PINSPECT
                set_opt["INSPECT"]=optarg
            case "-X"|"--test": options|=PTEST

            case "--style": set_opt["STYLESHEET"]=optarg
            case "--dark":  options|=DARKMODE

    return options, set_opt, args


def plot_one_time_runs(options, set_opt, args):
    if options & SHOWHELP:
        scr.usage(__doc__,verbose=options&VERBOSE)

        if options & PINSPECT: p_error(fn_pinspect.__doc__)

        return EXIT_EARLY

    if _file_name := set_opt.get("STYLESHEET"):
        load_style(_file_name)

    if options & DARKMODE:
        load_style("%s/extras/dark.style"%starbug2.__path__[0])
    
    if options & STOPPROC: return EXIT_EARLY
    if options & KILLPROC:
        p_error("..killing process\n")
        return EXIT_FAIL

    return EXIT_SUCCESS

def fn_pinspect(options, set_opt, images=None, tables=None):
    """
    Inspect Source
    --------------

    Plot at a source position cutouts in a range of images. 
    This requires a source list to be loaded, a list of image 
    file and the source catalogue number to be given. This will 
    take the form::

        $~ starbug2-plot -I CN123 sourcelist.fits image*.fits
    """
    """
    Parameters
    ----------
    options : int
        The starbug2.bin.plot options integar

    setopt : dict
        The starbug2.bin.plot setopt dictionary

    images : list
        The list of fits image HDUs to cut out from 

    tables : list
        The source list to pull the source from. Must have a column
        with the name "Catalogue_Number"

    Returns
    -------
    fig : plt.figure
        The output figure
    """
    fig=None
    if (cn:=set_opt.get("INSPECT")) and images and tables:
        if ("Catalogue_Number" in tables[0].col_names
                and cn in tables[0]["Catalogue_Number"]):
            i = np.where(tables[0]["Catalogue_Number"]==cn)[0]
            fig = plot_inspectsource(tables[0][i], images)

    else: p_error(
        "Must include the source Catalogue_Number,"
        " a list of images and a source list \n")
    return fig

def plot_main(argv):
    warn("Still in development\n\n")
    options,setopt,args=plot_parse_argv(argv)
    load_style("%s/extras/starbug.style"%starbug2.__path__[0])

    if options or setopt:
        if exit_code := plot_one_time_runs(options, setopt, args):
            return exit_code

    images=[]
    tables=[]
    for arg in args:
        if _file_name:=os.path.exists(arg):
            fp=fits.open(arg)
            _filter=fp[0].header.get("FILTER") # THIS IS A HACK
            for hdu in fp:
                if hdu.header.get("XTENSION") == "IMAGE":
                    images.append(hdu)
                    break
                if hdu.header.get("XTENSION") == "BINTABLE":
                    tables.append(Table(hdu.data))
                    break
            hdu.header["FILTER"]=_filter

    fig = None
    if options& PTEST:
        fig,ax=plt.subplots(1,figsize=(3,2.5))
        ax=plot_test(ax)

    if options& PINSPECT: fig=fn_pinspect(options, setopt, images=images, tables=tables)

    if fig is not None:
        fig.tight_layout()
        if output:=setopt.get("OUTPUT"):
            fig.savefig(output, dpi=300)
        else:
            plt.show()

def plot_main_entry():
    """Command Line entry point"""
    return plot_main(sys.argv)
