"""Copyright (C) 2026 UKATC

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>."""

# noinspection SpellCheckingInspection
"""StarbugII Plotting Scripts
usage: starbug2-plot [-vhX] [-I CN000] [-o outfile] images.fits
    -h  --help           : show help screen
    -o  --output   f_name : output filename
    -v  --verbose        : verbose mode

    -I  --inspect  CN000 : inspect a source in an array of images
    -X  --test           : plot a test image

        --style    f_name : load a custom pyplot style sheet
        --dark           : plot in dark mode
    -apfile              : ?????
"""
import os, sys, getopt
from typing import Final

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import starbug2
from starbug2.constants import (
    EXIT_EARLY, EXIT_FAIL, EXIT_SUCCESS, CAT_NUM, FILTER, EXT, IMAGE,
    BIN_TABLE, OUTPUT, INSPECT, STYLESHEET, AP_FILE_SET_OPT)
from starbug2.plot import load_style, plot_test, plot_inspect_source
from starbug2.utils import p_error, warn, parse_cmd, usage
from astropy.io.fits import PrimaryHDU, ImageHDU, BinTableHDU

VERBOSE: Final[int] = 0x01
SHOW_HELP: Final[int] = 0x02
STOP_PROC: Final[int] = 0x04
KILL_PROC: Final[int] = 0x08
DARK_MODE: Final[int] = 0x10

PTEST: Final[int] = 0x1000
PINSPECT: Final[int] = 0x2000


def plot_parse_argv(
        argv: list[str]) -> tuple[int,
                                  dict[str, float | int | str], list[str]]:
    """
    Parses configuration flags and image/catalogue parameters from the
     command line string.
    """
    options: int = 0
    set_opt: dict[str, float | int | str] = {}

    cmd: list[str]
    cmd, argv = parse_cmd(argv)

    # noinspection SpellCheckingInspection
    opts: list[tuple[str, str]]
    args: list[str]
    opts, args = getopt.gnu_getopt(
        argv, "hvXI:d:o:",
        ["help", "verbose", "test", "inspect=", "output=", "style=", "apfile",
         "dark"]
    )

    for opt, opt_arg in opts:
        # noinspection SpellCheckingInspection
        match opt:
            case "-h" | "--help":
                options |= (SHOW_HELP | STOP_PROC)
            case "-v" | "--verbose":
                options |= VERBOSE
            case "-o" | "--output":
                set_opt[OUTPUT] = opt_arg
            case "-d" | "--apfile":
                set_opt[AP_FILE_SET_OPT] = opt_arg
            case "-I" | "--inspect":
                options |= PINSPECT
                set_opt[INSPECT] = opt_arg
            case "-X" | "--test":
                options |= PTEST
            case "--style":
                set_opt[STYLESHEET] = opt_arg
            case "--dark":
                options |= DARK_MODE

    return options, set_opt, args


def plot_one_time_runs(
        options: int, set_opt: dict[str, Any], args: list[str]) -> int:
    """
    Handles initialization routines such as style sheet distribution or
     help menu warnings.
    """
    if options & SHOW_HELP:
        usage(__doc__, verbose=bool(options & VERBOSE))
        if options & PINSPECT:
            p_error(str(fn_pinspect.__doc__))
        return EXIT_EARLY

    # Only throw an error if files are missing when they are explicitly
    # required
    if not (options & PTEST) and len(args) == 0:
        p_error(
            "Error: Image or catalogue argument targets must be provided.\n")
        return EXIT_EARLY

    if _file_name := set_opt.get(STYLESHEET):
        load_style(str(_file_name))

    if options & DARK_MODE:
        load_style(f"{starbug2.__path__[0]}/extras/dark.style")

    if options & STOP_PROC:
        return EXIT_EARLY
    if options & KILL_PROC:
        p_error("..killing process\n")
        return EXIT_FAIL

    return EXIT_SUCCESS


def fn_pinspect(set_opt: dict[str, float | str | int],
                images: list[fits.HDUList | fits.ImageHDU] | None = None,
                tables: list[Table] | None = None) -> plt.Figure | None:
    """
    Plot cutouts at a source position across a range of images.

    This requires a source list to be loaded, a list of image
    files, and the source catalogue number to be given. This will
    take the form::

        $~ starbug2-plot -I CN123 source_list.fits image*.fits

    :param set_opt: The starbug2 config options dictionary
    :type set_opt: dict
    :param images: The list of FITS image HDUs to cut out from
    :type images: list
    :param tables: The source list containing coordinates matching a
                   'Catalogue_Number'
    :type tables: list of astropy.Table
    :return: The output figure object containing rendered cutouts
    :rtype: matplotlib.pyplot.Figure or None
    """
    fig: plt.Figure | None = None
    cn: str | None = set_opt.get(INSPECT)

    if cn and images and tables and len(tables) > 0:
        if CAT_NUM in tables[0].colnames and cn in tables[0][CAT_NUM]:
            i: np.ndarray = np.where(tables[0][CAT_NUM] == cn)[0]
            fig = plot_inspect_source(tables[0][i], images)
    else:
        p_error(
            f"Must include the source {CAT_NUM}, "
            f"a list of images and a source list \n"
        )
    return fig


def plot_main(argv: list[str]) -> int | None:
    """
    Main runtime entry path configuration structure for
    data visualization parsing loops.
    """
    warn("Still in development\n\n")
    options: int
    set_opt: dict[str, float | int | str]
    args: list[str]
    options, set_opt, args = plot_parse_argv(argv)

    load_style(f"{starbug2.__path__[0]}/extras/starbug.style")

    if options or set_opt:
        if (exit_code := plot_one_time_runs(
                options, set_opt, args)) != EXIT_SUCCESS:
            return exit_code

    images: list[PrimaryHDU | ImageHDU | BinTableHDU | None] = []
    tables: list[Table] = []

    for arg in args:
        if os.path.exists(arg):
            fp: fits.HDUList = fits.open(arg)
            _filter: str = fp[0].header.get(FILTER)

            # Use type tracking alias explicitly during extraction loop blocks
            hdu: PrimaryHDU | ImageHDU | BinTableHDU | None = None
            for hdu in fp:

                if hdu.header.get(EXT) == IMAGE:
                    images.append(hdu)
                    break
                if hdu.header.get(EXT) == BIN_TABLE:
                    tables.append(Table(hdu.data))
                    break
            if hdu is not None:
                hdu.header[FILTER] = _filter

    fig: plt.Figure | None = None

    if options & PTEST:
        ax: plt.Axes
        fig, ax = plt.subplots(1, figsize=(3, 2.5))
        plot_test(ax)

    if options & PINSPECT:
        fig = fn_pinspect(set_opt, images=images, tables=tables)

    if fig is not None:
        fig.tight_layout()
        if output := set_opt.get(OUTPUT):
            fig.savefig(str(output), dpi=300)
        else:
            plt.show()

    return EXIT_SUCCESS


def plot_main_entry() -> int | None:
    """Command Line package gateway binary endpoint entry pointer mapper."""
    return plot_main(sys.argv)