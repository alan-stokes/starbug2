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
import numpy as np
from astropy.table import vstack
from starbug2 import utils
from starbug2.constants import (
    PARAM_FILE_TAG, EXIT_EARLY, EXIT_SUCCESS, EXIT_FAIL, VERBOSE_TAG,
    CAT_NUM, MATCH_THRESH, FILTER, NEXP_THRESH, OUTPUT, MIRI, NIRCAM,
    match_cols)
from starbug2.filters import filters
from starbug2.matching import (
    GenericMatch, CascadeMatch, BandMatch, ExactValueMatch, band_match,
    parse_mask)
from starbug2 import param
import starbug2.bin as scr

# random bit trackers.
VERBOSE = 0x01
KILL_PROC = 0x02
STOP_PROC = 0x04
SHOW_HELP = 0x08

BAND_MATCH = 0x10
BAND_DEPR = 0x20
GENERIC_MATCH = 0x40
CASCADE_MATCH = 0x80
EXACT_MATCH = 0x100

EXP_FULL = 0x1000

# unique param tags
ERR_COL = "ERRORCOLUMN"
MASK_EVAL = "MASKEVAL"
MATCH_COLS = "MATCH_COLS"


def match_parse_m_argv(argv):
    """
    parses some arguments.

    :param argv:  the arg to parse.
    :return:
    """
    options = 0
    set_opt = {}

    cmd, argv = scr.parse_cmd(argv)
    opts, args = getopt.gnu_getopt(
        argv, "BCfGhvXe:m:o:p:s:",
        ["band", "cascade", "dither", "exact", "full", "generic", "help",
         VERBOSE_TAG.lower(), "error=", "mask=", "output=", "param=", "set=",
         "band-depr"]
    )
    for opt, optarg in opts:
        if opt in ("-h", "--help"):
            options |= (SHOW_HELP | STOP_PROC)
        if opt in ("-v", "--verbose"):
            options |= VERBOSE
        if opt in ("-o", "--output"):
            set_opt[OUTPUT] = optarg
        if opt in ("-p", "--param"):
            set_opt[PARAM_FILE_TAG] = optarg

        if opt in ("-e", "--error"):
            set_opt[ERR_COL] = optarg
        if opt in ("-f", "--full"):
            options |= EXP_FULL
        if opt in ("-m", "--mask"):
            set_opt[MASK_EVAL] = optarg
        if opt in ("-s", "--set"):
            if '=' in optarg:
                key, val = optarg.split('=')
                try:
                    val = float(val)
                except (ValueError, AttributeError, NameError):
                    pass
                set_opt[key] = val
            else: utils.p_error(
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

def match_one_time_runs(options, set_opt):
    """
    Options set, one time runs
    """
    if options & VERBOSE:
        set_opt[VERBOSE_TAG] = 1
    if options & SHOW_HELP:
        scr.usage(__doc__, verbose=options&VERBOSE)
        return EXIT_EARLY
    return EXIT_SUCCESS

def match_full_band_match(tables, parameters):
    utils.p_error("THIS NEEDS A TEST\n")
    to_match = { NIRCAM: [],
                 MIRI: [] }
    _col_names = ["RA","DEC","flag"]
    d_threshold = parameters.get(MATCH_THRESH)

    for i,tab in enumerate(tables):
        filter_string = tab.meta.get(FILTER)
        to_match[filters[filter_string].instr].append(tab)
        _col_names += ([filter_string, "e%s" % filter_string])
    
    if to_match[NIRCAM] and to_match[MIRI]:
        utils.printf("Detected NIRCam to MIRI matching\n")
        nircam_matched = band_match(
            to_match[NIRCAM], col_names=_col_names)
        miri_matched = band_match(
            to_match[MIRI], col_names=_col_names)

        load = utils.Loading(
            len(miri_matched),
            msg="Combining NIRCAM-MIRI(%.2g\")" % d_threshold)
        if bridge_col := parameters.get("BRIDGE_COL"):
            mask = np.isnan(nircam_matched[bridge_col])
            utils.printf("-> bridging catalogues with %s\n" % bridge_col)
        else:
            mask = np.full(len(nircam_matched), False)

        m = GenericMatch(threshold=d_threshold, load=load)
        full = m((nircam_matched[~mask], miri_matched))
        matched = m.finish_matching(full)
        matched.remove_column("NUM")
        matched = vstack((matched, nircam_matched[mask]))
    else:
        matched = band_match(tables, col_names=_col_names)
    return matched


def match_main(argv):
    """"""
    options, set_opt, args = match_parse_m_argv(argv)
    if options or set_opt:
        if ((exit_code := match_one_time_runs(options, set_opt)) !=
                EXIT_SUCCESS):
            return exit_code

    ##########
    # PARAMS #
    ##########
    if not (p_file := set_opt.get(PARAM_FILE_TAG)):
        if os.path.exists("./starbug.param"):
            p_file = "./starbug.param"
        else:
            p_file = None
    parameters = param.load_params(p_file)
    parameters.update(set_opt)

    #################
    # MAIN ROUTINES #
    #################

    tables=[ ]
    for f_name in args:
        t = utils.import_table(f_name, verbose=True)
        if t is not None:
            tables.append(t)
    if raw := parameters.get(MASK_EVAL):
        masks = [ parse_mask(raw,t) for t in tables ]
        for m in masks:
            try:
                print(m, sum(m), len(m))
            except (TypeError, NameError, ImportError):
                print( m )
    else:
        masks = None


    if len(tables) > 1:
        col_names = match_cols
        col_names += [ name for name in parameters[MATCH_COLS].split()
                       if name not in col_names]
        d_threshold = parameters[MATCH_THRESH]
        error_column = (
            set_opt.get(ERR_COL) if set_opt.get(ERR_COL) else "eflux")

        if options & BAND_DEPR:
            av = match_full_band_match(tables, parameters)
            full = None
            options &= ~EXP_FULL

        else:
            if options & BAND_MATCH:
                if isinstance(d_threshold, str):
                    d_threshold = np.array(
                        parameters[MATCH_THRESH].split(','), float)
                if parameters[FILTER] != "":
                    filter_string = parameters[FILTER].split(',')
                else:
                    filter_string = utils.remove_duplicates(
                        [utils.find_filter(t) for t in tables])
                matcher = BandMatch(
                    threshold=d_threshold, fltr=filter_string,
                    verbose=parameters[VERBOSE_TAG])

            elif options & CASCADE_MATCH:
                matcher = CascadeMatch(
                    threshold=d_threshold, colnames=col_names,
                    verbose=parameters[VERBOSE_TAG])
            elif options & GENERIC_MATCH:
                matcher = GenericMatch(
                    threshold=d_threshold, col_names=col_names,
                    verbose=parameters[VERBOSE_TAG])
            elif options & EXACT_MATCH:
                matcher = ExactValueMatch(
                    value=CAT_NUM, colnames=None,
                    verbose=parameters[VERBOSE_TAG])
            else: 
                matcher=GenericMatch(
                    threshold=d_threshold, verbose=parameters[VERBOSE_TAG])
                options |= EXP_FULL

            if options & VERBOSE:
                print("\n%s"%matcher)

            full = matcher.match( tables, join_type="or", mask=masks )
            av = matcher.finish_matching(
                full,
                num_thresh=parameters[NEXP_THRESH],
                zp_mag=parameters["ZP_MAG"],
                error_column=error_column)


        output = parameters.get(OUTPUT)
        if output is None or output == '.':
            output = utils.combine_file_names(
                [name for name in args], n_mismatch=100)
        d_name, f_name, ext = utils.split_file_name(output)

        suffix = ""
        if options & EXP_FULL:
            utils.export_table(
                full, f_name="%s/%sfull.fits" % (d_name, f_name))
            utils.printf("-> %s/%sfull.fits\n" % (d_name,f_name))
            suffix = "match"
        if av:
            utils.export_table(av,"%s/%s%s.fits" % (d_name,f_name,suffix))
            utils.printf("-> %s/%s%s.fits\n" % (d_name,f_name,suffix))

        return EXIT_SUCCESS
    
    elif len(tables) == 1:
        return EXIT_EARLY
    else:
        utils.p_error("No tables loaded for matching.\n")
        return EXIT_FAIL

def match_main_entry():
    """StarbugII-match entry"""
    return match_main(sys.argv)
