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
from astropy.io.ascii.cparser import AstropyWarning

from starbug2.misc import generate_runscript
import warnings
import os, sys

from astropy.io.fits import PrimaryHDU
from astropy.io.fits.header import Header

from starbug2.matching.generic_match import GenericMatch
from starbug2.star_bug_config import StarBugMainConfig
from starbug2.starbug import StarbugBase

# quietens astropy so that it doesn't flood the terminal with warnings.
# ABS this seems concerning, if they're producing warnings we should be
# exploring those.
warnings.simplefilter("ignore", category=AstropyWarning)
warnings.simplefilter("ignore", category=RuntimeWarning) ## bit dodge that

from starbug2.initialise_psf_data import init_starbug_for_jwst, generate_psf

from starbug2.constants import (
    DETECTION, BACKGROUND, APP_HOT, PSFP_HOT, MATCH_OUTPUTS, LOGO,
    HELP_STRINGS,  EXIT_EARLY, EXIT_SUCCESS, EXIT_FAIL, EXIT_MIXED,
    READ_THE_DOCS_URL, FITS_EXTENSION)
from starbug2.utils import (
    p_error, printf, warn, split_file_name, export_region,
    combine_file_names, export_table, puts, parse_cmd,
    usage, get_version)
from starbug2 import param
from astropy.table import Table

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

# noinspection SpellCheckingInspection
sys.stdout.write("\x1b[1mlaunching \x1b[36mstarbug\x1b[0m\n")

# noinspection SpellCheckingInspection
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
        config.generate_main_get_opt_definitions())

    _, argv = parse_cmd(argv)
    config.populate_params(
        argv, short_definition, long_definition, config.MAIN_FLAG_MAP)
    return config

def starbug_one_time_runs(config: StarBugMainConfig) -> int:
    """
    Options set, verify/run one time functions
    """

    if config.show_version:
        printf(get_version())

    if config.show_help:
        usage(__doc__, verbose=config.verbose_logs)

        if config.do_star_detection:
            p_error(HELP_STRINGS[DETECTION])
        if config.do_bgd_estimate:
            p_error(HELP_STRINGS[BACKGROUND])
        if config.do_aperture_photometry:
            p_error(HELP_STRINGS[APP_HOT])
        if config.do_photometry_routine:
            p_error(HELP_STRINGS[PSFP_HOT])
        if config.do_matching:
            p_error(HELP_STRINGS[MATCH_OUTPUTS])
        return EXIT_EARLY

    ## Load parameter files for onetime runs
    parameter_file: str | None
    if (parameter_file := config.param_file) is None:
        if os.path.exists("./starbug.param"):
            parameter_file = "starbug.param"
        else:
            parameter_file = None

    config.load_params(parameter_file)

    if config.update_param:
        param.update_param_file(parameter_file)
        return EXIT_SUCCESS

    output: int | float | str
    if _output := config.output_file:
        _output: int | float | str
        output = _output
    else:
        output = '.'

    #########################
    # One time run commands #
    #########################

    ## Initialise or update starbug
    if config.execute_jwst_initialisation:
        init_starbug_for_jwst()

    ## Generate a single PSF
    if config.generate_psf:
        if config.got_valid_psf_generation_params():
            filter_string: str | None = config.custom_filter
            assert filter_string is not None
            detector: str| None = config.detector_name
            psf_size: int = config.psf_fit_size
            printf(
                "Generating PSF: %s %s (%d)\n" %
                (filter_string, detector, psf_size))
            psf: PrimaryHDU | None = generate_psf(
                filter_string, detector=detector, fov_pixels=psf_size)
            if psf: 
                name: str = (
                    "%s%s.fits" %
                      (filter_string, "" if detector is None else detector))
                printf("--> %s\n" % name)
                psf.writeto(name, overwrite=True)
            else:
                p_error("PSF Generation failed :(\n")
        else:
            # noinspection SpellCheckingInspection
            p_error(
                "Unable to generate PSF. Set filter with '-s FILTER=FXXX and "
                "Set detector name with '-s DET_NAME=XXX and "
                "Set psf_fit_size with '-s PSF_SIZE=XXX'\n")

    ## Generate a run script
    if config.generate_run:
        generate_runscript(config.fits_images, "starbug2 ")
        if not config.fits_images:
            p_error("no files included to create runscript with\n")

    ## Generate a region from a table
    if config.generate_region:
        file_name: str | None = config.region_file
        if file_name and os.path.exists(file_name):
            table: Table = Table.read(file_name, format="fits")
            _, name, _ = split_file_name(file_name)
            name: str
            export_region(
                table, colour=config.region_colour,
                scale_radius=config.region_scale,
                region_radius=config.region_radius,
                x_col=config.region_x_column_name,
                y_col=config.region_y_column_name,
                wcs=config.region_uses_wcs,
                f_name="%s/%s.reg" % (output, name))
            printf("generating region --> %s/%s.reg\n"%(output,name))

    # generate local param file as requested
    if config.generate_local_param_file:
        config.do_generate_local_param_file()

    return EXIT_SUCCESS


def starbug_match_outputs(
        starbugs: list[StarbugBase | None], config: StarBugMainConfig) -> None:
    """
    Matching output catalogues

    :param starbugs: star bug instances
    :type starbugs: list of starbugBase or None
    :param config: the config object
    :type config: StarBugMainConfig
    :return: None
    """
    if config.verbose_logs:
        printf("Matching outputs\n")

    f_name: str | None

    # filter out any Nones.
    valid_bugs: list[StarbugBase] = [sb for sb in starbugs if sb is not None]

    # get file name
    if f_name := combine_file_names([sb.f_name for sb in valid_bugs]):
        _, name ,_ = split_file_name(os.path.basename(f_name))
        name: str
        f_name = "%s/%s"%(valid_bugs[0].out_dir, name)
    else:
        f_name = "out"

    header: Header = valid_bugs[0].header

    match: GenericMatch = GenericMatch(
        threshold = config.match_threshold_arc_sec_as_an_arc_sec,
        col_names = None,
        p_file = config.param_file)

    if config.do_star_detection or config.do_aperture_photometry:
        full: Table = match(
            [sb.detections for sb in valid_bugs], join_type="or")
        av: Table = match.finish_matching(
            full, num_thresh=config.exposure_count_threshold,
            zp_mag=config.zero_point_magnitude)

        printf("-> %s-ap*...\n" % f_name)

        # noinspection SpellCheckingInspection
        export_table(full, f_name="%s-apfull.fits" % f_name, header=header)

        # noinspection SpellCheckingInspection
        export_table(av, f_name="%s-apmatch.fits" % f_name, header=header)

    if config.do_photometry_routine:
        full: Table = match(
            [sb.psf_catalogue for sb in valid_bugs], join_type="or")
        av: Table = match.finish_matching(
            full, num_thresh=config.exposure_count_threshold,
            zp_mag=config.zero_point_magnitude)

        printf("-> %s-psf*...\n" % f_name)

        # noinspection SpellCheckingInspection
        export_table(full, f_name="%s-psffull.fits" % f_name, header=header)

        # noinspection SpellCheckingInspection
        export_table(av, f_name="%s-psfmatch.fits" % f_name, header=header)


def execute_star_bug(
        args: tuple[str, StarBugMainConfig, bool]) -> StarbugBase | None:
    """
    Worker function to initialise and run standard photometry processes on a
    single file.

    :param args: A tuple containing (file_name, config, use_verbose)
    :type args: tuple
    :return: The verified StarbugBase pipeline wrapper instance, or None
             if validation fails
    :rtype: starbug2.StarbugBase or None
    """
    ## I've put this here because it takes some time
    from starbug2.starbug import StarbugBase
    star_bug_base: StarbugBase | None = None
    f_name: str
    config: StarBugMainConfig
    f_name, config, use_verbose = args
    if os.path.exists(f_name):
        folder, file_name, ext = split_file_name(f_name)

        ap_file: str | None = config.ap_file
        background_file: str | None = config.background_file

        if config.find_file:
            ap: str = "%s/%s-ap.fits" % (folder, file_name)
            bgd: str = "%s/%s-bgd.fits" % (folder, file_name)
            if os.path.exists(ap)  and config.ap_file is None:
                ap_file = ap
            if os.path.exists(bgd) and config.background_file is None:
                background_file = bgd

        ## Sorting out the stdout
        if use_verbose:
            printf("-> showing starbug stdout for \"%s\"\n" % f_name)
        elif config.n_cores > 1:
            printf("-> hiding starbug stdout for \"%s\"\n" % f_name)
        else:
            printf("-> %s\n" % f_name)

        if ext == FITS_EXTENSION:
            star_bug_base: StarbugBase = StarbugBase(
                f_name, config=config, ap_file=ap_file,
                bkg_file=background_file, verbose=use_verbose)
            if star_bug_base.verify():
                warn("System verification failed\n")
                return None

            if config.do_star_detection:
                star_bug_base.detect()
            if config.do_bgd_estimate:
                star_bug_base.bgd_estimate()
            if config.do_bgd_subtraction:
                star_bug_base.bgd_subtraction()
            if config.do_source_geometry:
                star_bug_base.source_geometry()
            if config.do_aperture_photometry:
                star_bug_base.aperture_photometry()
            if config.do_photometry_routine or config.generate_residual_image:
                star_bug_base.photometry_routine()

        else:
            p_error("file must be type '.fits' not %s\n" % ext)
    else:
        p_error("can't access %s\n" % f_name)
    return star_bug_base


def starbug_main(argv: list[str]) -> int:
    """
    Command-line execution orchestrator for processing astronomical image
    datasets.

    :param argv: System arguments mapping configurations and input filenames
    :type argv: list of str
    :return: System operational termination exit code status matrix
    :rtype: int
    """
    config: StarBugMainConfig = starbug_parse_argv(argv)

    if config.use_main_one_time_runs():
       return starbug_one_time_runs(config)

    if config.fits_images:
        # why import here
        import starbug2
        from multiprocessing import Pool
        from multiprocessing.pool import Pool as PoolType

        # freeze the config now to avoid writers
        config.freeze()

        puts(LOGO % READ_THE_DOCS_URL)
        exit_code: int = EXIT_SUCCESS
        starbugs: list[StarbugBase | None]

        if ((n_cores := config.n_cores) is None
                or n_cores == 1 or len(config.fits_images) == 1):

            config.unfreeze()
            config.n_cores = 1
            config.freeze()

            starbugs = (
                [execute_star_bug(
                    (file_name, config, config.verbose_logs))
                    for file_name in config.fits_images])
        else:
            pool: PoolType = Pool(processes=n_cores)

            # this ensures only the first worker executes verbose.
            worker_tasks = [
                (file_name, config, index == 0)
                for index, file_name in enumerate(config.fits_images)
            ]
            starbugs = pool.map(execute_star_bug, worker_tasks)
            pool.close()
            pool.join()

        to_remove: list[StarbugBase | None] = []
        sb: StarbugBase | None
        for n, sb in enumerate(starbugs):
            if not sb: 
                p_error("FAILED: %s\n" % config.fits_images[n])
                to_remove.append(sb)
                exit_code = EXIT_MIXED
        for sb in to_remove:
            starbugs.remove(sb)

        if not starbug2:
            exit_code = EXIT_FAIL

            
        if config.do_matching and len(starbugs) > 1:
            starbug_match_outputs(starbugs, config)
        

    else:
        p_error("fits image file must be included\n")
        exit_code = EXIT_FAIL

    return exit_code

def starbug_main_entry() -> int:
    """
    System binary path gateway routing console script entries.
    """
    return starbug_main(sys.argv)
