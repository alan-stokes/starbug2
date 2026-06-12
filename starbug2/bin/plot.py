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
import os, sys
from typing import Final

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table
import starbug2
from starbug2.constants import (
    EXIT_EARLY, EXIT_SUCCESS, CAT_NUM, FILTER, EXT, IMAGE, BIN_TABLE)
from starbug2.plot import load_style, plot_test, plot_inspect_source
from starbug2.star_bug_config import StarBugMainConfig
from starbug2.utils import p_error, warn, parse_cmd, usage
from astropy.io.fits import PrimaryHDU, ImageHDU, BinTableHDU

VERBOSE: Final[int] = 0x01
SHOW_HELP: Final[int] = 0x02
STOP_PROC: Final[int] = 0x04
KILL_PROC: Final[int] = 0x08
DARK_MODE: Final[int] = 0x10

PTEST: Final[int] = 0x1000
PINSPECT: Final[int] = 0x2000


def starbug_parse_argv(argv: list[str]) -> StarBugMainConfig:
    """
    Organise the sys argv line into options, values and arguments

    :param argv: the arguments
    :return: the config class
    :rtype: StarBugMainConfig
    """
    config: StarBugMainConfig = StarBugMainConfig()
    short_definition: str
    long_definition: list[str]
    short_definition, long_definition = (
        config.generate_plot_get_opt_definitions())

    _, argv = parse_cmd(argv)
    config.populate_params(
        argv, short_definition, long_definition, config.PLOT_FLAG_MAP)
    return config


def plot_one_time_runs(config: StarBugMainConfig) -> int:
    """
    Handles initialisation routines such as style sheet distribution or
     help menu warnings.
    """
    if config.show_plot_help:
        usage(__doc__, verbose=bool(config.verbose_logs))
        if config.inspect_parameter:
            p_error(str(fn_pinspect.__doc__))
        return EXIT_EARLY

    # Only throw an error if files are missing when they are explicitly
    # required
    if not config.test_mode and len(config.fits_images) == 0:
        p_error(
            "Error: Image or catalogue argument targets must be provided.\n")
        return EXIT_EARLY

    if config.plot_style is not None:
        load_style(str(config.plot_style))

    if config.dark_mode:
        load_style(f"{starbug2.__path__[0]}/extras/dark.style")

    return EXIT_SUCCESS


def fn_pinspect(
        inspect_string: str,
        images: list[fits.PrimaryHDU | fits.ImageHDU |
                     fits.BinTableHDU | None] | None = None,
        tables: list[Table] | None = None) -> plt.Figure | None:
    """
    Plot cutouts at a source position across a range of images.

    This requires a source list to be loaded, a list of image
    files, and the source catalogue number to be given. This will
    take the form::

        $~ starbug2-plot -I CN123 source_list.fits image*.fits

    :param inspect_string: the string to inspect
    :type inspect_string: str
    :param images: The list of FITS image HDUs to cut out from
    :type images: list
    :param tables: The source list containing coordinates matching a
                   'Catalogue_Number'
    :type tables: list of astropy.Table
    :return: The output figure object containing rendered cutouts
    :rtype: matplotlib.pyplot.Figure or None
    """
    fig: plt.Figure | None = None

    if inspect_string and images and tables and len(tables) > 0:
        if (CAT_NUM in tables[0].colnames
                and inspect_string in tables[0][CAT_NUM]):
            i: np.ndarray = np.where(tables[0][CAT_NUM] == inspect_string)[0]
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
    data visualisation parsing loops.
    """
    warn("Still in development\n\n")

    config: StarBugMainConfig = starbug_parse_argv(argv)

    load_style(f"{starbug2.__path__[0]}/extras/starbug.style")

    if exit_code := plot_one_time_runs(config) != EXIT_SUCCESS:
        return exit_code

    images: list[PrimaryHDU | ImageHDU | BinTableHDU | None] = []
    tables: list[Table] = []

    for arg in config.fits_images:
        if os.path.exists(arg):
            fp: fits.HDUList = fits.open(arg)
            _filter: str = fp[0].header.get(FILTER)

            # Use type tracking alias explicitly during extraction loop blocks
            hdu: PrimaryHDU | ImageHDU | BinTableHDU | None = None
            for hdu in fp:
                if hdu is None:
                    continue

                if hdu.header.get(EXT) == IMAGE:
                    images.append(hdu)
                    break
                if hdu.header.get(EXT) == BIN_TABLE:
                    tables.append(Table(hdu.data))
                    break
            if hdu is not None:
                hdu.header[FILTER] = _filter

    fig: plt.Figure | None = None

    if config.test_mode:
        ax: plt.Axes
        fig, ax = plt.subplots(1, figsize=(3, 2.5))
        plot_test(ax)

    inspect_string: str | None = config.inspect_parameter
    if inspect_string is not None:
        fig = fn_pinspect(inspect_string, images=images, tables=tables)

    if fig is not None:
        fig.tight_layout()
        if output := config.output_file:
            fig.savefig(str(output), dpi=300)
        else:
            plt.show()

    return EXIT_SUCCESS


def plot_main_entry() -> int | None:
    """Command Line package gateway binary endpoint entry pointer mapper."""
    return plot_main(sys.argv)