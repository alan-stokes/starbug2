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
"""StarbugII Matching
usage: starbug2-match [-BCGfhvX] [-e column] [-m mask] [-o output]
                      [-p file.param] [-s KEY=VAL] table.fits ...
    -B  --band               : match in "BAND" mode (does not preserve a
                               column for every frame)
    -C  --cascade            : match in "CASCADE" mode (left justify columns)
    -G  --generic            : match in "GENERIC" mode
    -X  --exact              : match in "EXACTVALUE" mode

    -e  --error   column     : photometric error column ("eflux" or "stdflux")
    -f  --full               : export full catalogue
    -h  --help               : show help message
    -m  --mask    eval       : column evaluation to mask out of matching e.g.
                               -m"~np.isnan(F444W)"
    -o  --output  file.fits  : output matched catalogue
    -p  --param   file.param : load starbug parameter file
    -s  --set     option     : set value in parameter file at runtime
                               (-s MATCH_THRESH=1)
    -v  --verbose            : display verbose outputs

        --band-depr          : match in "old" band mode

    --> typical runs
       $~ starbug2-match -Gfo outfile.fits tab1.fits tab2.fits
       $~ starbug2-match -sMATCH_THRESH=0.2 -sBRIDGE_COL=F444W -Bo out.fits
                         F*W.fits
"""
import os, sys, getopt
from typing import Final, Any

import numpy as np
from astropy.table import Table, vstack
from starbug2 import utils, param
from starbug2.constants import (
    PARAM_FILE_TAG, EXIT_EARLY, EXIT_SUCCESS, EXIT_FAIL, VERBOSE_TAG,
    CAT_NUM, MATCH_THRESH, FILTER, NEXP_THRESH, OUTPUT, STAR_BUG_MIRI, NIRCAM,
    match_cols, RA, DEC, FLAG, BRIDGE_COL, NUM, E_FLUX)
from starbug2.filters import STAR_BUG_FILTERS
from starbug2.matching.band_match import BandMatch
from starbug2.matching.cascade_match import CascadeMatch
from starbug2.matching.exact_value_match import ExactValueMatch
from starbug2.matching.generic_match import GenericMatch
from starbug2.misc import parse_mask
from starbug2.utils import parse_cmd, usage

# Random bit trackers.
VERBOSE: Final[int] = 0x01
KILL_PROC: Final[int] = 0x02
STOP_PROC: Final[int] = 0x04
SHOW_HELP: Final[int] = 0x08

BAND_MATCH: Final[int] = 0x10
BAND_DEPR: Final[int] = 0x20
GENERIC_MATCH: Final[int] = 0x40
CASCADE_MATCH: Final[int] = 0x80
EXACT_MATCH: Final[int] = 0x100

EXP_FULL: Final[int] = 0x1000

# Unique param tags
# noinspection SpellCheckingInspection
ERR_COL: Final[str] = "ERRORCOLUMN"
# noinspection SpellCheckingInspection
MASK_EVAL: Final[str] = "MASKEVAL"
MATCH_COLS: Final[str] = "MATCH_COLS"


def match_parse_m_argv(
        argv: list[str]) -> tuple[int, dict[str, Any], list[str]]:
    """
    Parses CLI command arguments for catalogue matching operations.

    :param argv: List of system arguments passed via the terminal execution
                 string
    :type argv: list of str
    :return: A parsed sequence consisting of options_bit_mask,
            configuration_dict, and raw arguments
    :rtype: tuple
    """
    options: int = 0
    set_opt: dict[str, float | int | str] = {}

    cmd: list[str]
    cmd, argv = parse_cmd(argv)

    opts: list[tuple[str, str]]
    args: list[str]
    opts, args = getopt.gnu_getopt(
        argv, "BCfGhvXe:m:o:p:s:",
        ["band", "cascade", "dither", "exact", "full", "generic", "help",
         VERBOSE_TAG.lower(), "error=", "mask=", "output=", "param=", "set=",
         "band-depr"]
    )

    for opt, opt_arg in opts:
        if opt in ("-h", "--help"):
            options |= (SHOW_HELP | STOP_PROC)
        if opt in ("-v", "--verbose"):
            options |= VERBOSE
        if opt in ("-o", "--output"):
            set_opt[OUTPUT] = opt_arg
        if opt in ("-p", "--param"):
            set_opt[PARAM_FILE_TAG] = opt_arg

        if opt in ("-e", "--error"):
            set_opt[ERR_COL] = opt_arg
        if opt in ("-f", "--full"):
            options |= EXP_FULL
        if opt in ("-m", "--mask"):
            set_opt[MASK_EVAL] = opt_arg
        if opt in ("-s", "--set"):
            if '=' in opt_arg:
                key, val_raw = opt_arg.split('=')
                key: str
                val_raw: str
                val: float | str
                val = val_raw
                try:
                    val = float(val_raw)
                except (ValueError, AttributeError, NameError):
                    pass
                set_opt[key] = val
            else:
                utils.p_error(
                    "unable to set parameter, use syntax -s KEY=VALUE\n")

        if opt in ("-B", "--band"):
            options |= BAND_MATCH
        if opt in ("-C", "--cascade"):
            options |= CASCADE_MATCH
        if opt in ("-G", "--generic"):
            options |= GENERIC_MATCH
        if opt in ("-X", "--exact"):
            options |= EXACT_MATCH
        if opt == "--band-depr":
            options |= BAND_DEPR
    return options, set_opt, args


def match_one_time_runs(
        options: int, set_opt: dict[str, float | int | str]) -> int:
    """
    Executes immediate utility operations such as version checks or help menus.
    """
    if options & VERBOSE:
        set_opt[VERBOSE_TAG] = 1
    if options & SHOW_HELP:
        usage(__doc__, verbose=bool(options & VERBOSE))
        return EXIT_EARLY
    return EXIT_SUCCESS


def match_full_band_match(
        tables: list[Table],
        parameters: dict[str, float | int | str]) -> Table:
    """
    Handles fallback deprecated band-matching configurations across diverse detectors.
    """
    utils.p_error("THIS NEEDS A TEST\n")
    to_match: dict[int, list[Table]] = {
        NIRCAM: [],
        STAR_BUG_MIRI: []
    }
    _col_names: list[str] = [RA, DEC, FLAG]
    d_threshold: float = float(parameters.get(MATCH_THRESH))
    band_matcher: BandMatch = BandMatch(threshold=d_threshold)

    for tab in tables:
        tab: Table
        filter_string: str = str(tab.meta.get(FILTER))
        to_match[STAR_BUG_FILTERS[filter_string].instr].append(tab)
        _col_names += [filter_string, f"e{filter_string}"]

    matched: Table
    if to_match[NIRCAM] and to_match[STAR_BUG_MIRI]:
        utils.printf("Detected NIRCam to MIRI matching\n")
        nir_cam_matched: Table = band_matcher.band_match(
            to_match[NIRCAM], col_names=_col_names)
        miri_matched: Table = band_matcher.band_match(
            to_match[STAR_BUG_MIRI], col_names=_col_names)

        # noinspection SpellCheckingInspection
        load: utils.Loading = utils.Loading(
            len(miri_matched),
            msg="Combining NIRCAM-MIRI(%.2g\")" % d_threshold
        )

        mask: np.ndarray
        if bridge_col := parameters.get(BRIDGE_COL):
            mask = np.isnan(nir_cam_matched[bridge_col])
            utils.printf("-> bridging catalogues with %s\n" % bridge_col)
        else:
            mask = np.full(len(nir_cam_matched), False)

        m: GenericMatch = GenericMatch(threshold=d_threshold, load=load)
        full: Any = m((nir_cam_matched[~mask], miri_matched))
        matched = m.finish_matching(full)
        matched.remove_column(NUM)
        matched = vstack((matched, nir_cam_matched[mask]))
    else:
        matched = band_matcher.band_match(tables, col_names=_col_names)

    return matched


def match_main(argv: list[str]) -> int:
    """
    Main runtime processing loop for executing cross-catalogue astronomical
    source coordinate matching.
    """
    options: int
    set_opt: dict[str, float | int | str]
    args: list[str]
    options, set_opt, args = match_parse_m_argv(argv)

    if options or set_opt:
        if (exit_code := match_one_time_runs(options, set_opt)) != EXIT_SUCCESS:
            return exit_code

    p_file: str | None = set_opt.get(PARAM_FILE_TAG)
    if not p_file:
        p_file = (
            "./starbug.param" if os.path.exists("./starbug.param") else None)

    parameters: dict[str, float | int | str] = param.load_params(p_file)
    parameters.update(set_opt)

    tables: list[Table] = []
    for f_name in args:
        t: Table | None = utils.import_table(f_name, verbose=True)
        if t is not None:
            tables.append(t)

    masks: list[np.ndarray] | None = None
    if raw := parameters.get(MASK_EVAL):
        masks = [parse_mask(raw, t) for t in tables]
        for m in masks:
            try:
                print(m, sum(m), len(m))
            except (TypeError, NameError, ImportError):
                print(m)

    if len(tables) > 1:
        col_names: list[str] = list(match_cols)
        col_names += [
            name for name in str(parameters[MATCH_COLS]).split()
            if name not in col_names
        ]
        d_threshold: float | int | str | np.ndarray = parameters[MATCH_THRESH]
        error_column: str = (
            str(set_opt.get(ERR_COL) if set_opt.get(ERR_COL) else E_FLUX))

        av: Table
        full: Table| None = None
        matcher: GenericMatch

        if options & BAND_DEPR:
            av = match_full_band_match(tables, parameters)
            options &= ~EXP_FULL
        else:
            if options & BAND_MATCH:
                if isinstance(d_threshold, str):
                    d_threshold = np.array(
                        parameters[MATCH_THRESH].split(','), float)

                filter_string: list[str]
                if parameters[FILTER] != "":
                    filter_string = str(parameters[FILTER]).split(',')
                else:
                    filter_string = utils.remove_duplicates(
                        [utils.find_filter(t) for t in tables])

                matcher = BandMatch(
                    threshold=d_threshold, fltr=filter_string,
                    verbose=parameters[VERBOSE_TAG]
                )
            elif options & CASCADE_MATCH:
                matcher = CascadeMatch(
                    threshold=d_threshold, colnames=col_names,
                    verbose=parameters[VERBOSE_TAG]
                )
            elif options & GENERIC_MATCH:
                matcher = GenericMatch(
                    threshold=d_threshold, col_names=col_names,
                    verbose=parameters[VERBOSE_TAG]
                )
            elif options & EXACT_MATCH:
                matcher = ExactValueMatch(
                    value=CAT_NUM, colnames=None,
                    verbose=parameters[VERBOSE_TAG]
                )
            else:
                matcher = GenericMatch(
                    threshold=d_threshold, verbose=parameters[VERBOSE_TAG]
                )
                options |= EXP_FULL

            if options & VERBOSE:
                print("\n%s" % matcher)

            full = matcher.match(tables, join_type="or", mask=masks)
            av = matcher.finish_matching(
                full,
                num_thresh=parameters[NEXP_THRESH],
                zp_mag=parameters["ZP_MAG"],
                error_column=error_column
            )

        output: str | None = parameters.get(OUTPUT)
        if output is None or output == '.':
            output = utils.combine_file_names(
                [name for name in args], n_mismatch=100)

        d_name: str
        f_name: str
        ext: str
        d_name, f_name, ext = utils.split_file_name(output)

        suffix: str = ""
        if options & EXP_FULL:
            utils.export_table(
                full, f_name="%s/%sfull.fits" % (d_name, f_name))
            utils.printf("-> %s/%sfull.fits\n" % (d_name, f_name))
            suffix = "match"

        if av:
            utils.export_table(av, "%s/%s%s.fits" % (d_name, f_name, suffix))
            utils.printf("-> %s/%s%s.fits\n" % (d_name, f_name, suffix))

        return EXIT_SUCCESS

    elif len(tables) == 1:
        return EXIT_EARLY
    else:
        utils.p_error("No tables loaded for matching.\n")
        return EXIT_FAIL


def match_main_entry() -> int:
    """StarbugII-match entry path map setup routing wrapper."""
    return match_main(sys.argv)