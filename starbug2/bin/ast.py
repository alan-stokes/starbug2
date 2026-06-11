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
import numpy as np
import glob
from multiprocessing import Pool, Process, shared_memory
from itertools import repeat
from time import sleep
from astropy.table import Table

from starbug2.constants import (
    PARAM_FILE_TAG, N_CORES, OUTPUT, N_TESTS, N_STARS, AUTO_SAVE, QUIETMODE,
    EXIT_EARLY, EXIT_FAIL, EXIT_SUCCESS, MAX_MAG, MIN_MAG, FLUX, FLUX_DET,
    PLOTAST)
from starbug2.starbug import StarbugBase
from starbug2.artificialstars import ArtificialStarsIII, compile_results
from starbug2.utils import (
    printf, p_error, combine_tables, fill_nan, translate_param_float,
    parse_cmd, usage)
from starbug2.param import load_params

# random bit markers
VERBOSE = 0x01
SHOW_HELP = 0x02
STOP_PROC = 0x04
KILL_PROC = 0x08
NO_BGD = 0x10
NO_PHOT = 0x20
RECOVER = 0x40

# globals
c = np.array([0, 0, 0], dtype=np.int64)
share = shared_memory.SharedMemory(create=True, size=c.nbytes)
buffer = np.ndarray(c.shape, dtype=c.dtype, buffer=share.buf)


def load():
    """
    A loading bar that should be run in a subprocess
    It sits and watches the shared memory buffer and periodically
    prints out a progress bar
    """
    global buffer
    while buffer[0] < buffer[1]:
        sleep(1)
        p = buffer[0] / buffer[1]
        msg = f"recovering:{buffer[2]}%"
        s = "\x1b[2K%s|%-40s|%d/%d\r" % (
            msg, int(p*40)*'=', int(buffer[0]), int(buffer[1]))
        printf(s)
        sys.stdout.flush()
    printf("\n")

def ast_parse_argv(argv):
    """ Organise the argv line into options, values and arguments """
    options = 0
    set_opt = {QUIETMODE : 1, AUTO_SAVE : 100}
    cmd, argv = parse_cmd(argv)
    # noinspection SpellCheckingInspection
    opts, args = getopt.gnu_getopt(
        argv, "hvN:n:p:R:S:s:o:",
        ["help", "verbose", "ncores=", "param=", "set=", "output=",
         "ntests=", "nstars=", "autosave=", "no-background", "no-psfphot",
         "recover"])

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
            opt, opt_arg, set_opt, options, KILL_PROC
        )

    return options, set_opt, args

def ast_one_time_runs(options, args):
    """
    Set options, verify run and execute one time functions
    """

    if options & SHOW_HELP:
        usage(__doc__, verbose=options & VERBOSE)
        return EXIT_EARLY

    if options & RECOVER:
        if not args:
            # noinspection SpellCheckingInspection
            f_names = glob.glob("sbast-autosave*.tmp")
        else:
            f_names = [a for a in args if os.path.exists(a)]
        if f_names:
            printf("Recovery Mode:\n-> %s\n"%("\n-> ".join(f_names)))
            raw = Table()
            for f_name in f_names:
                raw = combine_tables(raw, Table.read(f_name))
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

def fn(args):
    f_name, options, set_opt, index = args
    global buffer
    out = None
    if os.path.exists(f_name):
        star_bug_base = StarbugBase(
            f_name, set_opt.get(PARAM_FILE_TAG), options=set_opt)
        opt = star_bug_base.options
        ast = ArtificialStarsIII(star_bug_base, index=index)
        out = ast.auto_run(
            opt.get(N_TESTS), stars_per_test=opt.get(N_STARS),
            mag_range=(opt.get(MAX_MAG),opt.get(MIN_MAG)),
            loading_buffer=buffer, autosave=opt.get(AUTO_SAVE),
            skip_phot=options & NO_PHOT, skip_background=options & NO_BGD)
    return out

def ast_main(argv):
    global buffer, share
    options, set_opt, args = ast_parse_argv(argv)
    exit_code = EXIT_SUCCESS

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
        f_name = args[0]
        n_tests = params.get(N_TESTS)
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
        loading = Process(target=load, args=())
        loading.start()

        if (n_cores := params.get(N_CORES)) is None or n_cores == 1:
            params[N_CORES] = 1
            outs = [fn((f_name, options, params, 0)) for f_name in args]
        else:
            n_cores = min(n_cores, n_tests)
            zip_options = np.full(n_cores, options, dtype=int)
            for n in range(n_cores):
                if n > 0:
                    zip_options[n] &= ~VERBOSE
            params[N_TESTS] = int(np.ceil(n_tests / n_cores))
            params[AUTO_SAVE] = int(np.ceil(set_opt.get(AUTO_SAVE) / n_cores))


            pool = Pool(processes=n_cores)
            outs = pool.map(
                fn, zip(
                    repeat(f_name), zip_options, repeat(params),
                    range(1, n_cores + 1)))
            pool.close()

        #force finish
        buffer[0] = buffer[1]
        loading.join()
        
        #############################
        # COMPILING ALL THE RESULTS #
        #############################

        raw = outs[0]
        for res in outs[1:]:
            raw = combine_tables(raw, res)
        star_bug_base = StarbugBase(
            f_name, set_opt.get(PARAM_FILE_TAG), options=set_opt)
        if options & VERBOSE:
            printf("-> compiling results\n")
            printf("-> flux recovery: %.2g\n" % (
                np.nanmean(raw[FLUX] / raw[FLUX_DET])))

        if (results := compile_results(
                raw, image=star_bug_base.main_image,
                filter_string=star_bug_base.filter,
                plot_ast=set_opt.get(PLOTAST))):
            out_dir, b_name, _= StarbugBase.sort_output_names(
                f_name, param_output=set_opt.get(OUTPUT))
            if options & VERBOSE:
                printf("--> %s/%s-ast.fits\n" % (out_dir, b_name))
            results.writeto("%s/%s-ast.fits"%(out_dir, b_name), overwrite=True)

            ## autosave cleanup
            # noinspection SpellCheckingInspection
            for _f_name in glob.glob("sbast-autosave*.tmp"):
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

def ast_main_entry():
    """Command line entry point"""
    return ast_main(sys.argv)
