# noinspection SpellCheckingInspection
"""StarbugII - JWST PSF photometry
usage: starbug2 [-ABDfGhMPSv] [-b bgdfile] [-d apfile] [-n ncores] [-o ouput]
                [-p file.param] [-s opt=val] image.fits ...
   -A  --apphot          : run aperture photometry on a source list
   -B  --background      : run background estimation
   -D  --detect          : run source detection
   -G  --geom            : calculate geometric stats on source list
   -M  --match           : match outputs from all input image files
   -P  --psf             : run psf photometry
   -S  --subbgd          : subtract background from image

   -b  --bgdfile         : load background (-bgd.fits) file
   -d  --apfile  ap.fits : load a source detection (-ap.fits) file to skip the
                           source detection step
   -f  --find            : attempt to find associated -ap -bgd files
   -h  --help            : display uasage information
   -n  --ncores      num : number of CPU cores to split process between
   -o  --output      dir : output directory
   -p  --param   a.param : load parameter file
   -s  --set      option : set value in parameter file at runtime (-s SIGSKY=3)
   -v  --verbose         : display verbose outputs

   --> Single run commands
       --init                     : Initialise Starbug (post install)
       --local-param              : Make a local copy of the default
                                    parameter file
       --update-param             : Update an out-of-date local parameter file
       --generate-psf             : Generate a single PSF. Set FILTER,
                                    DET_NAME, PSF_SIZE with -s
       --generate-region   a.fits : Make a ds9 region file with a detection
                                    file
       --generate-run      *.fits : Generate a simple run script
       --version                  : Print starbug2 version

   --> typical runs
      //Source detect on image with a parameter file
      $~ starbug2 -vD -p file.param image.fits
      //Source detect and match outputs of a list of images
      $~ starbug2 -vDM -n4 images*.fits
      //PSF photometry on an image with a source file (image-ap.fits)
      $~ starbug2 -vd image-ap.fits -BP image.fits

To see more detailed information on an option, run [OPTION] --help:
    $~ starbug2 -D --help

See https://starbug2.readthedocs.io for full documentation.

"""
import os, sys, getopt

import numpy as np

from starbug2.constants import (
    SHOWHELP, STOPPROC, VERBOSE, PARAM_FILE_TAG, DOAPPHOT, DOBGDEST, DODETECT,
    DOGEOM, DOMATCH, DOPHOTOM, DOBGDSUB, DOARTIFL, FINDFILE, KILLPROC, INITSB,
    GENRATPSF, UPDATEPRM, GENRATRUN, GENRATREG, REGION_TAB, DETECTION,
    BACKGROUND, APP_HOT, PSFP_HOT, MATCH_OUTPUTS, OUTPUT, APPLYZP, CALCINSTZP,
    LOGO, HELP_STRINGS, N_CORES, EXIT_EARLY, EXIT_SUCCESS, EXIT_FAIL,
    EXIT_MIXED, READ_THE_DOCS_URL, FILTER, DET_NAME, PSF_SIZE, MATCH_THRESH, NEXP_THRESH, ZP_MAG, AP_FILE, BGD_FILE,
    FITS_EXTENSION, REGION_COL, REGION_SCAL, REGION_RAD, REGION_X_COL, REGION_Y_COL, REGION_WCS)
from starbug2.utils import (
    p_error, printf, get_version, warn, split_file_name, export_region,
    combine_file_names, export_table, puts)
from starbug2 import param
import starbug2.bin as scr
from astropy.table import Table

# noinspection SpellCheckingInspection
sys.stdout.write("\x1b[1mlaunching \x1b[36mstarbug\x1b[0m\n")

# noinspection SpellCheckingInspection
def starbug_parse_argv(argv):
    """
    Organise the sys argv line into options, values and arguments

    :param argv: the arguments
    :return: tuple containing (options, set_opt, args)
    :rtype: tuple int, dict of string, string, list of str
    """
    options = 0
    set_opt = {}

    cmd, argv = scr.parse_cmd(argv)

    opts, args = getopt.gnu_getopt(
        argv,
        "ABDfGhMPSvb:d:n:o:p:s:",
        [
            "apphot","background", "detect", "find", "geom", "help",
            "match", "psf", "subbgd", "verbose", "xtest",
            "bgdfile=", "apfile=", "ncores=", "output=", "param=", "set=",
            "init", "generate-psf", "local-param", "generate-region=",
            "version", "generate-run", "update-param", "debug", "dev"
        ]
    )

    for opt, opt_arg in opts:
        if opt in ("-h", "--help"):
            options |= (SHOWHELP | STOPPROC)
        if opt in ("-p", "--param"):
            set_opt[PARAM_FILE_TAG] = opt_arg
        if opt in ("-v", "--verbose"):
            options |= VERBOSE

        if opt in ("-A", "--apphot"):
            options |= DOAPPHOT
        if opt in ("-B", "--background"):
            options |= DOBGDEST
        if opt in ("-D", "--detect"):
            options |= DODETECT
        if opt in ("-G", "--geom"):
            options |= DOGEOM
        if opt in ("-M", "--match"):
            options |= DOMATCH
        if opt in ("-P", "--psf"):
            options |= DOPHOTOM
        if opt in ("-S", "--subbgd"):
            options |= DOBGDSUB

        if opt == "--dev":
            options |= DOARTIFL

        if opt in ("-d", "--apfile"):
            if os.path.exists(opt_arg):
                set_opt["AP_FILE"] = opt_arg
            else:
                p_error("AP_FILE \"%s\" does not exist\n" % opt_arg)

        if opt in ("-b", "--bgdfile"):
            if os.path.exists(opt_arg):
                set_opt["BGD_FILE"]=opt_arg
            else:
                p_error("BGD_FILE \"%s\" does not exist\n" % opt_arg)

        if opt in ("-f", "--find"):
            options |= FINDFILE
        if opt in ("-n", "--ncores"):
            set_opt["NCORES"] = max(1,int(opt_arg))

        if opt in ("-o", "--output"):
            set_opt["OUTPUT"] = opt_arg

        if opt in ("-s", "--set"):
            if '=' in opt_arg:
                key, val = opt_arg.split('=')
                try:
                    val = float(val)
                except ValueError:
                    pass
                set_opt[key] = val
            else:
                p_error("unable to set parameter, use syntax -s KEY=VALUE\n")
                options |= KILLPROC

        if opt == "--init":
            options |= ( INITSB | STOPPROC)
        if opt == "--generate-psf":
            options |= (GENRATPSF | STOPPROC)
        if opt == "--update-param":
            options |= (UPDATEPRM | STOPPROC)
        if opt == "--generate-run":
            options |= (GENRATRUN | STOPPROC)
        if opt == "--generate-region":
            set_opt[REGION_TAB] = opt_arg
            options |= (GENRATREG | STOPPROC)

        if opt == "--local-param":
            param.local_param()
            printf("--> generating starbug.param\n")
            options |= STOPPROC

        if opt == "--version":
            printf(LOGO % ("starbug2-v%s" % get_version()))
            options |= STOPPROC
    return options, set_opt, args

def starbug_one_time_runs(options, set_opt, args):
    """
    Options set, verify/run one time functions
    """

    # ABS why are we only importing these here?
    from starbug2.misc import init_starbug, generate_psf, generate_runscript

    if options & SHOWHELP:
        scr.usage(__doc__, verbose=options & VERBOSE)

        if options & DODETECT:
            p_error(HELP_STRINGS[DETECTION])
        if options & DOBGDEST:
            p_error(HELP_STRINGS[BACKGROUND])
        if options & DOAPPHOT:
            p_error(HELP_STRINGS[APP_HOT])
        if options & DOPHOTOM:
            p_error(HELP_STRINGS[PSFP_HOT])
        if options & DOMATCH:
            p_error(HELP_STRINGS[MATCH_OUTPUTS])
        return EXIT_EARLY

    ## Load parameter files for onetime runs
    if (p_file := set_opt.get(PARAM_FILE_TAG)) is None:
        if os.path.exists("./starbug.param"):
            p_file = "starbug.param"
        else:
            p_file = None

    init_parameters = param.load_params(p_file)

    if options & UPDATEPRM:
        param.update_param_file(p_file)
        return EXIT_EARLY

    tmp = param.load_default_params()
    if (set(tmp.keys()) - set(init_parameters.keys())
            | set(init_parameters.keys()) - set(tmp.keys())):
        warn("Parameter file version mismatch. "
             "Run starbug2 --update-param to update\nquitting :(\n")
        return EXIT_FAIL

    init_parameters.update(set_opt)
    if _output := init_parameters.get(OUTPUT):
        output = _output
    else:
        output = '.'

    #########################
    # One time run commands #
    #########################

    ## Initialise or update starbug
    if options & INITSB:
        init_starbug()

    ## Generate a single PSF
    if options & GENRATPSF:
        if filter_string := init_parameters.get(FILTER):
            detector = init_parameters.get(DET_NAME)
            psf_size = init_parameters.get(PSF_SIZE)
            printf(
                "Generating PSF: %s %s (%d)\n" %
                (filter_string, detector, psf_size))
            psf = generate_psf(
                filter_string, detector=detector, fov_pixels=psf_size)
            if psf: 
                name = (
                    "%s%s.fits" %
                      (filter_string, "" if detector is None else detector))
                printf("--> %s\n" % name)
                psf.writeto(name, overwrite=True)
            else:
                p_error("PSF Generation failed :(\n")
        else:
            # noinspection SpellCheckingInspection
            p_error(
                "Unable to generate PSF. Set filter with '-s FILTER=FXXX'\n")

    ## Generate a run script
    if options & GENRATRUN:
        generate_runscript(args, "starbug2 ")
        if not args:
            p_error("no files included to create runscript with\n")

    ## Generate a region from a table
    if options & GENRATREG:
        file_name = set_opt.get("REGION_TAB")
        if file_name and os.path.exists(file_name):
            table = Table.read(file_name, format="fits")
            _, name, _ = split_file_name(file_name)
            export_region(
                table, colour=init_parameters[REGION_COL],
                scale_radius=init_parameters[REGION_SCAL],
                region_radius=init_parameters[REGION_RAD],
                x_col=init_parameters[REGION_X_COL],
                y_col=init_parameters[REGION_Y_COL],
                wcs=init_parameters[REGION_WCS],
                f_name="%s/%s.reg" % (output, name))
            printf("generating region --> %s/%s.reg\n"%(output,name))

    ###########################
    # instrumental zero point #
    ###########################
    if options & (APPLYZP | CALCINSTZP):
        p_error("instrumental zero point application deprecated\n")

    if options & STOPPROC:
        return EXIT_EARLY ## quiet ending the process if required

    if options & KILLPROC:
        p_error("..quitting :(\n\n")
        return scr.usage(__doc__, verbose=options&VERBOSE)

    return EXIT_SUCCESS


def starbug_match_outputs(starbugs, options, set_opt):
    """
    Matching output catalogues

    :param starbugs: star bug instances
    :param options: options dict
    :param set_opt: other options.
    :return: None
    """
    from starbug2.matching import GenericMatch
    if options & VERBOSE:
        printf("Matching outputs\n")
    params = param.load_params(set_opt.get(PARAM_FILE_TAG))
    params.update(set_opt)

    if f_name := combine_file_names([sb.fname for sb in starbugs]):
        _, name ,_ = split_file_name(os.path.basename(f_name))
        f_name = "%s/%s"%(starbugs[0].out_dir, name)
    else:
        f_name = "out"

    header = starbugs[0].header

    match = GenericMatch(
        threshold = params[MATCH_THRESH],
        col_names = None,
        p_file = set_opt.get(PARAM_FILE_TAG))

    if options & (DODETECT | DOAPPHOT):
        full = match( [sb.detections for sb in starbugs], join_type="or")
        av = match.finish_matching(
            full, num_thresh=params[NEXP_THRESH], zpmag=params[ZP_MAG])

        printf("-> %s-ap*...\n" % f_name)

        # noinspection SpellCheckingInspection
        export_table(full, f_name="%s-apfull.fits" % f_name, header=header)

        # noinspection SpellCheckingInspection
        export_table(av, f_name="%s-apmatch.fits" % f_name, header=header)

    if options & DOPHOTOM:
        full = match( [sb.psfcatalogue for sb in starbugs], join_type="or")
        av = match.finish_matching(
            full, num_thresh=params[NEXP_THRESH], zpmag=params[ZP_MAG])

        printf("-> %s-psf*...\n" % f_name)

        # noinspection SpellCheckingInspection
        export_table(full, f_name="%s-psffull.fits" % f_name, header=header)

        # noinspection SpellCheckingInspection
        export_table(av, f_name="%s-psfmatch.fits" % f_name, header=header)


def fn(args):
    """
    ??????
    :param args: the args
    :return: the star bug instance for the function
    :rtype starbug2.StarBugBase
    """
    ## I've put this here because it takes some time
    from starbug2.starbug import StarbugBase
    star_bug_base = None
    f_name, options, set_opt = args
    if os.path.exists(f_name):
        folder, file_name, ext = split_file_name(f_name)

        if options & FINDFILE:
            ap = "%s/%s-ap.fits" % (folder,file_name)
            bgd = "%s/%s-bgd.fits" % (folder,file_name)
            if os.path.exists(ap)  and not set_opt.get(AP_FILE):
                set_opt[AP_FILE] = ap
            if os.path.exists(bgd) and not set_opt.get(BGD_FILE):
                set_opt[BGD_FILE] = bgd

        ## Sorting out the stdout
        if options & VERBOSE:
            printf("-> showing starbug stdout for \"%s\"\n" % f_name)
            set_opt[VERBOSE] = 1
        elif set_opt.get(N_CORES) > 1:
            printf("-> hiding starbug stdout for \"%s\"\n" % f_name)
        else: printf("-> %s\n" % f_name)

        if ext == FITS_EXTENSION:
            star_bug_base = StarbugBase(
                f_name, p_file=set_opt.get(PARAM_FILE_TAG), options=set_opt)
            if star_bug_base.verify():
                warn("System verification failed\n")
                return None

            if options & DODETECT:
                star_bug_base.detect()
            if options & DOBGDEST:
                star_bug_base.bgd_estimate()
            if options & DOBGDSUB:
                star_bug_base.bgd_subtraction()
            if options & DOGEOM:
                star_bug_base.source_geometry()

            if options & DOAPPHOT:
                star_bug_base.aperture_photometry()
            if options & DOPHOTOM:
                star_bug_base.photometry()

            if options & DOARTIFL:
                star_bug_base.artificial_stars()
        else:
            p_error("file must be type '.fits' not %s\n" % ext)
    else:
        p_error("can't access %s\n" % f_name)
    return star_bug_base


def starbug_main(argv):
    """Command entry"""
    options, set_opt, args= starbug_parse_argv(argv)

    if options or set_opt:

        if exit_code := starbug_one_time_runs(options, set_opt, args):
            return exit_code

    if args:
        # why import here
        import starbug2
        from multiprocessing import Pool
        from itertools import repeat

        puts(LOGO % READ_THE_DOCS_URL)
        exit_code = EXIT_SUCCESS

        if ((n_cores := set_opt.get(N_CORES)) is None
                or n_cores == 1 or len(args) == 1):
            set_opt[N_CORES] = 1
            starbugs = (
                [fn((file_name, options, set_opt)) for file_name in args])
        else:

            zip_options = np.full(len(args),options, dtype=int)
            for n in range(len(args)):
                if n > 0:
                    zip_options[n] &= ~VERBOSE

            pool = Pool(processes=n_cores)
            starbugs = pool.map(fn,zip(args, zip_options, repeat(set_opt)))
            pool.close()

        for n,sb in enumerate(starbugs): 
            if not sb: 
                p_error("FAILED: %s\n" % args[n])
                starbugs.remove(sb)
                exit_code = EXIT_MIXED

        if not starbug2:
            exit_code = EXIT_FAIL

            
        if options & DOMATCH and len(starbugs) > 1:
            starbug_match_outputs(starbugs, options, set_opt)
        

    else:
        p_error("fits image file must be included\n")
        exit_code = EXIT_FAIL

    return exit_code

def starbug_main_entry():
    """Entry point"""
    return starbug_main(sys.argv)
