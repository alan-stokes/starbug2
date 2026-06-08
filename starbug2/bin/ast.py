# noinspection SpellCheckingInspection

"""
StarbugII Artificial Star Testing
usage: starbug2-ast [-vhR] [-N ntests] [-n ncores] [-p file.param] [-S nstars]
                    [-s opt=val] image.fits ...
    -h  --help          : show help screen
    -N  --ntests    num : number of tests to run
    -n  --ncores  cores : number of cores to split the tests over
    -o  --output output : output directory or filename to export results to
    -p  --param    file : load a parameter file
    -R  --recover       : recover incomplete test autosave files
    -S  --nstars    num : number of stars to inject per test
    -s  --set    option : set parameter at runtime with syntax "-s KEY=VALUE"
    -v  --verbose       : show verbose stdout output

        --autosave freq : frequency of quick save outputs
        --no-background : turn off background estimation routine
        --no-psfphot    : turn off psf photometry routine
"""

import os,sys,getopt
from multiprocessing.shared_memory import SharedMemory
from typing import Final, Dict, Tuple

import numpy as np
import glob
from multiprocessing import Pool, Process, shared_memory
from itertools import repeat
from time import sleep
from astropy.table import Table
from astropy.io.fits import HDUList

from starbug2.constants import (
    PARAM_FILE_TAG, N_CORES, OUTPUT, N_TESTS, N_STARS, AUTO_SAVE, QUIETMODE,
    EXIT_EARLY, EXIT_FAIL, EXIT_SUCCESS, MAX_MAG, MIN_MAG, FLUX, FLUX_DET,
    PLOTAST)
from starbug2.starbug import StarbugBase
from starbug2.artificialstars import ArtificialStars, compile_results
from starbug2.utils import (
    printf, p_error, combine_tables, fill_nan, translate_param_float,
    parse_cmd, usage)
from starbug2.param import load_params

# random bit markers
VERBOSE: Final[int] = 0x01
SHOW_HELP: Final[int] = 0x02
STOP_PROC: Final[int] = 0x04
KILL_PROC: Final[int] = 0x08
NO_BGD: Final[int] = 0x10
NO_PHOT: Final[int] = 0x20
RECOVER: Final[int] = 0x40

# globals
c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
share: SharedMemory = shared_memory.SharedMemory(create=True, size=c.nbytes)
buffer: np.ndarray = np.ndarray(c.shape, dtype=c.dtype, buffer=share.buf)


def load() -> None:
    """
    A loading bar that should be run in a subprocess
    It sits and watches the shared memory buffer and periodically
    prints out a progress bar
    """
    global buffer
    while buffer[0] < buffer[1]:
        sleep(1)
        p: np.ndarray = buffer[0] / buffer[1]
        msg: str = f"recovering:{buffer[2]}%"
        s: str = "\x1b[2K%s|%-40s|%d/%d\r" % (
            msg, int(p*40)*'=', int(buffer[0]), int(buffer[1]))
        printf(s)
        sys.stdout.flush()
    printf("\n")


def ast_parse_argv(argv: list[str]) -> (
        Tuple[int, Dict[str, int | str | float], list[str]]):
    """ Organise the argv line into options, values and arguments """
    options: int = 0
    set_opt: Dict[str, int | str | float] = {QUIETMODE : 1, AUTO_SAVE : 100}
    cmd, argv = parse_cmd(argv)
    cmd: str
    argv: list[str]

    # noinspection SpellCheckingInspection
    opts, args = getopt.gnu_getopt(
        argv, "hvN:n:p:R:S:s:o:",
        ["help", "verbose", "ncores=", "param=", "set=", "output=",
         "ntests=", "nstars=", "autosave=", "no-background", "no-psfphot",
         "recover"])
    opts: list[tuple[str, str]]
    args: list[str]

    for opt, opt_arg in opts:
        if opt in ("-h","--help"):
            options |= (SHOW_HELP | STOP_PROC)
        if opt in ("-v","--verbose"):
            options |= VERBOSE
        if opt in ("-p","--param"):
            set_opt[PARAM_FILE_TAG] = opt_arg
        # noinspection SpellCheckingInspection
        if opt in ("-n","--ncores"):
            set_opt[N_CORES] = int(opt_arg)
        if opt in ("-o","--output"):
            set_opt[OUTPUT] = opt_arg
        # noinspection SpellCheckingInspection
        if opt in ("-N","--ntests"):
            set_opt[N_TESTS] = int(opt_arg)
        # noinspection SpellCheckingInspection
        if opt in ("-S","--nstars"):
            set_opt[N_STARS] = int(opt_arg)

        # set options
        if opt in ("-R","--recover"):
            options |= (RECOVER | STOP_PROC)
        if opt == "--autosave":
            set_opt[AUTO_SAVE] = int(opt_arg)
        if opt == "--no-background" :
            options |= NO_BGD
        # noinspection SpellCheckingInspection
        if opt == "--no-psfphot"    :
            options |= NO_PHOT

        options, set_opt = translate_param_float(
            opt, opt_arg, set_opt, options, KILL_PROC)

    return options, set_opt, args

def ast_one_time_runs(options: int, args: list[str]) -> int:
    """
    Set options, verify run and execute one time functions
    """

    if options & SHOW_HELP:
        usage(__doc__, verbose=options & VERBOSE)
        return EXIT_EARLY

    if options & RECOVER:
        f_names: list[str] | None
        if not args:
            # noinspection SpellCheckingInspection
            f_names = glob.glob("sbast-autosave*.tmp")
        else:
            f_names = [a for a in args if os.path.exists(a)]
        if f_names:
            printf("Recovery Mode:\n-> %s\n"%("\n-> ".join(f_names)))
            raw: Table = Table()
            for f_name in f_names:
                f_name: str
                raw = combine_tables(raw, Table.read(f_name))
            results: HDUList
            if (results := compile_results(
                    fill_nan(raw), plot_ast="recovered.pdf")):
                printf("-> successful recovery!\n--> %s\n" % (
                    f_name := "recovered.fits"))
                results.writeto(f_name, overwrite=True)
            else:
                p_error("something went wrong\n")
        else:
            p_error("No files found to recover\n")



    if options & STOP_PROC:
        return EXIT_EARLY
    if options & KILL_PROC:
        p_error("..killing process\n")
        return EXIT_FAIL

    return EXIT_SUCCESS

def execute_artificial_stars(
        args: tuple[str, int, dict[str, int | str | float], int]) -> (
            Table | None):
    """
    Multiprocessing worker function to run artificial star tests on a given
    file.

    :param args: A tuple containing (f_name, options_flags,
                configuration_dict, worker_index)
    :type args: tuple
    :return: The generated artificial stars recovery catalogue table, or
             None if the file doesn't exist
    :rtype: astropy.table.Table or None
    """
    f_name: str
    options: int
    set_opt: dict[str, int | str | float]
    index: int
    f_name, options, set_opt, index = args

    global buffer
    out: Table | None = None
    if os.path.exists(f_name):
        star_bug_base: StarbugBase = StarbugBase(
            f_name, set_opt.get(PARAM_FILE_TAG), options=set_opt)
        opt: dict[str, int | str | float] = star_bug_base.options
        ast: ArtificialStars = ArtificialStars(star_bug_base, index=index)
        out = ast(
            opt.get(N_TESTS), stars_per_test=opt.get(N_STARS),
            mag_range=(opt.get(MAX_MAG),opt.get(MIN_MAG)),
            loading_buffer=buffer, autosave=opt.get(AUTO_SAVE),
            skip_phot=options & NO_PHOT, skip_background=options & NO_BGD)
    return out

def ast_main(argv: list[str]) -> int:
    global buffer, share

    options: int
    set_opt: dict[str, int | str | float]
    args: list[str]
    options, set_opt, args = ast_parse_argv(argv)

    exit_code: int = EXIT_SUCCESS

    if options or set_opt:
        if exit_code := ast_one_time_runs(options, args):
            share.unlink()
            return exit_code

    if params := load_params(set_opt.get(PARAM_FILE_TAG)):
        params.update(set_opt)
    else: 
        p_error("Failed to load parameters from file\n")
        return EXIT_FAIL

    if args:
        f_name: str = args[0]
        n_tests: int = int(params.get(N_TESTS))
        if options & VERBOSE:
            printf("Artificial Stars\n----------------\n")
            printf("-> loading %s\n"%f_name)
            if set_opt.get(PARAM_FILE_TAG):
                printf("-> parameters: %s\n" % set_opt.get(PARAM_FILE_TAG))
            printf("-> running %d tests with %d injections per test\n" % (
                n_tests, params.get(N_STARS)))
            printf("-> magnitude range: %.1f - %.1f\n" % (
                params.get(MAX_MAG), params.get(MIN_MAG)))
            if options & NO_PHOT:
                printf("-> skipping PSF photometry step\n")
            if options & NO_BGD:
                printf("-> skipping background estimation step\n")

        buffer[0] = 0
        buffer[1] = n_tests
        loading: Process = Process(target=load, args=())
        loading.start()

        # Initialise output container tracking tables
        outs: list[Table | None]

        if (n_cores := params.get(N_CORES)) is None or n_cores == 1:
            params[N_CORES] = 1
            outs = [execute_artificial_stars(
                (f_name, options, params, 0)) for f_name in args]
        else:
            n_cores: int = int(min(n_cores, n_tests))
            zip_options: np.ndarray = np.full(n_cores, options, dtype=int)
            for n in range(n_cores):
                if n > 0:
                    zip_options[n] &= ~VERBOSE
            params[N_TESTS] = int(np.ceil(n_tests / n_cores))
            params[AUTO_SAVE] = int(np.ceil(set_opt.get(AUTO_SAVE) / n_cores))


            pool: Pool = Pool(processes=n_cores)
            outs = pool.map(
                execute_artificial_stars, zip(
                    repeat(f_name), zip_options, repeat(params),
                    range(1, n_cores + 1)))
            pool.close()
            pool.join()

        #force finish
        buffer[0] = buffer[1]
        loading.join()
        
        #############################
        # COMPILING ALL THE RESULTS #
        #############################

        raw: Table = outs[0]
        for res in outs[1:]:
            raw = combine_tables(raw, res)
        star_bug_base: StarbugBase = StarbugBase(
            f_name, set_opt.get(PARAM_FILE_TAG), options=set_opt)
        if options & VERBOSE:
            printf("-> compiling results\n")
            printf("-> flux recovery: %.2g\n" % (
                np.nanmean(raw[FLUX] / raw[FLUX_DET])))

        results: HDUList
        if (results := compile_results(
                raw, image=star_bug_base.main_image.data,
                filter_string=star_bug_base.filter,
                plot_ast=set_opt.get(PLOTAST))):
            out_dir: str
            b_name: str
            out_dir, b_name, _= StarbugBase.sort_output_names(
                f_name, param_output=set_opt.get(OUTPUT))
            if options & VERBOSE:
                printf("--> %s/%s-ast.fits\n" % (out_dir, b_name))
            results.writeto("%s/%s-ast.fits"%(out_dir, b_name), overwrite=True)

            ## autosave cleanup
            # noinspection SpellCheckingInspection
            for _f_name in glob.glob("sbast-autosave*.tmp"):
                _f_name: str
                os.remove(_f_name)

        else:
            p_error("results compilation failed\n")

    else:
        p_error("must include a fits image to work on\n")
        exit_code = EXIT_FAIL

    # Wrapped fix to handle rapid multiprocess teardowns safely
    try:
        share.unlink()
    except FileNotFoundError:
        # The memory handle was already unlinked safely by another thread
        pass
    return exit_code

def ast_main_entry() -> int:
    """Command line entry point"""
    return ast_main(sys.argv)
