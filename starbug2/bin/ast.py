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
import os
import sys
from multiprocessing.shared_memory import SharedMemory
from multiprocessing import Pool, Process, shared_memory
from multiprocessing.pool import Pool as PoolType

import numpy as np
import glob
from time import sleep
from astropy.table import Table
from astropy.io.fits import HDUList

from starbug2.core.constants import ExitStates, TableColumn
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.core.artificialstars import ArtificialStars, compile_results
from starbug2.utilities.utils import (
    printf, p_error, combine_tables, fill_nan,  parse_cmd, usage)

import photutils

# Force photutils to strictly return standard QTables globally
photutils.future_column_names = True


def load(loading_buffer: np.ndarray) -> None:
    """
    A loading bar that should be run in a subprocess
    It sits and watches the shared memory buffer and periodically
    prints out a progress bar
    """
    while loading_buffer[0] < loading_buffer[1]:
        sleep(1)
        p: np.ndarray = loading_buffer[0] / loading_buffer[1]
        msg: str = f"recovering:{loading_buffer[2]}%"
        s: str = "\x1b[2K%s|%-40s|%d/%d\r" % (
            msg, int(p*40)*'=', int(loading_buffer[0]), int(loading_buffer[1]))
        printf(s)
        sys.stdout.flush()
    printf("\n")


def ast_parse_argv(argv: list[str]) -> StarBugMainConfig:
    """
    Organise the sys argv line into options, values and arguments

    :param argv: the arguments
    :return: the config class
    :rtype: StarBugMainConfig
    """
    cmd, argv = parse_cmd(argv)
    argv: list[str]
    short_definition: str
    long_definition: list[str]

    config: StarBugMainConfig = StarBugMainConfig()
    short_definition, long_definition = (
        config.generate_ast_get_opt_definitions())
    _, argv = parse_cmd(argv)
    config.populate_params(
        argv, short_definition, long_definition, config.AST_FLAG_MAP)
    return config


def ast_one_time_runs(config: StarBugMainConfig) -> ExitStates:
    """
    Set options, verify run and execute one time functions
    """

    if config.show_ast_help:
        usage(__doc__, verbose=config.verbose_logs)
        return ExitStates.EXIT_EARLY

    if config.ast_recover:
        f_names: list[str] | None
        if not config.fits_images:
            # noinspection SpellCheckingInspection
            f_names = glob.glob("sbast-autosave*.tmp")
        else:
            f_names = [a for a in config.fits_images if os.path.exists(a)]
        if f_names:
            printf("Recovery Mode:\n-> %s\n" % ("\n-> ".join(f_names)))
            raw: Table | None = Table()
            for f_name in f_names:
                f_name: str
                read_table: Table | None = Table.read(f_name)
                if read_table is None:
                    p_error(f"failed to read table at path {f_name}")
                    return ExitStates.EXIT_FAIL
                raw = combine_tables(raw, read_table)
            results: HDUList
            assert raw is not None
            if (results := compile_results(
                    fill_nan(raw), plot_ast="recovered.pdf")):
                printf("-> successful recovery!\n--> %s\n" % (
                    f_name := "recovered.fits"))
                results.writeto(f_name, overwrite=True)
            else:
                p_error("something went wrong\n")
        else:
            p_error("No files found to recover\n")
    return ExitStates.EXIT_SUCCESS


def execute_artificial_stars(
        f_name: str, config: StarBugMainConfig, verbose: bool,
        index: int, test_count: int, ast_auto_save: int,
        loading_buffer: np.ndarray) -> Table | None:
    """
    Multiprocessing worker function to run artificial star tests on a given
    file.
    :param f_name: the file to process
    :type f_name: str
    :param config: the config object
    :type config: StarBugMainConfig
    :param verbose: bool flag if to use verbose
    :type verbose: bool
    :param index: the index
    :type index: int
    :param test_count: the amount of tests
    :type test_count: int
    :param ast_auto_save: how many tests between saves
    :type ast_auto_save: int.
    :param loading_buffer: the loading buffer
    :type loading_buffer: np.ndarray
    :return: The generated artificial stars recovery catalogue table, or
             None if the file doesn't exist.
    :rtype: astropy.table.Table or None.
    """
    out: Table | None = None
    if os.path.exists(f_name):
        star_bug_base: StarbugBase = StarbugBase(
            f_name, config, ap_file=config.ap_file,
            bkg_file=config.background_file, verbose=verbose)
        ast: ArtificialStars = ArtificialStars(star_bug_base, index=index)
        out = ast.execute_ast(
            test_count,
            stars_per_test=config.stars_per_artificial_test,
            mag_range=(
                config.test_magnitude_bright_limit,
                config.test_magnitude_faint_limit),
            loading_buffer=loading_buffer,
            autosave=ast_auto_save,
            skip_phot=config.ast_no_psf_phot,
            skip_background=config.ast_no_background,
            zp_mag=config.zero_point_magnitude,
            sub_image_size=config.sub_image_crop_size,
            save_image=config.save_added_image,
            save_image_path=config.save_added_image_path,
            ast_seed=config.ast_seed)
    return out


def ast_main(
        argv: list[str], share_memory: SharedMemory,
        loading_buffer: np.ndarray) -> ExitStates:

    config: StarBugMainConfig = ast_parse_argv(argv)

    exit_code: ExitStates = ExitStates.EXIT_SUCCESS

    if config.use_ast_one_time_runs():
        if exit_code := ast_one_time_runs(config):
            share_memory.unlink()
            return exit_code
    config.freeze()

    print(f"{config.fits_images}")

    if config.fits_images:
        f_name: str = config.fits_images[0]
        n_tests: int = int(config.artificial_star_tests_count)
        if config.verbose_logs:
            printf("Artificial Stars\n----------------\n")
            printf("-> loading %s\n" % f_name)
            if config.param_file:
                printf("-> parameters: %s\n" % config.param_file)
            printf("-> running %d tests with %d injections per test\n" % (
                n_tests, config.stars_per_artificial_test))
            printf("-> magnitude range: %.1f - %.1f\n" % (
                config.test_magnitude_bright_limit,
                config.test_magnitude_faint_limit))
            if config.ast_no_psf_phot:
                printf("-> skipping PSF photometry step\n")
            if config.ast_no_background:
                printf("-> skipping background estimation step\n")

        loading_buffer[0] = 0
        loading_buffer[1] = n_tests
        loading: Process = Process(target=load, args=[loading_buffer])
        loading.start()

        # Initialise output container tracking tables
        outs: list[Table | None]

        if (n_cores := config.n_cores) is None or n_cores == 1:
            config.unfreeze()
            config.n_cores = 1
            config.freeze()
            outs = ([execute_artificial_stars(
                f_name, config, config.verbose_logs, index,
                config.artificial_star_tests_count, config.ast_auto_save,
                loading_buffer)
                    for index, f_name in enumerate(config.fits_images)])
        else:
            n_cores: int = int(min(n_cores, n_tests))
            per_process_n_test: int = int(np.ceil(n_tests / n_cores))
            per_process_tests_per_save: int = int(
                np.ceil(config.ast_auto_save / n_cores))

            worker_tasks = [
                (file_name, config, index == 0, index, per_process_n_test,
                 per_process_tests_per_save, loading_buffer)
                for index, file_name in enumerate(config.fits_images)
            ]

            pool: PoolType = Pool(processes=n_cores)
            outs = pool.starmap(execute_artificial_stars, worker_tasks)
            pool.close()
            pool.join()

        # force finish
        loading_buffer[0] = loading_buffer[1]
        loading.join()

        #############################
        # COMPILING ALL THE RESULTS #
        #############################

        raw: Table | None = outs[0]
        for res in outs[1:]:
            raw = combine_tables(raw, res)
        assert raw is not None
        star_bug_base: StarbugBase = StarbugBase(
            f_name, config, ap_file=config.ap_file,
            bkg_file=config.background_file, verbose=config.verbose_logs)
        if config.verbose_logs:
            printf("-> compiling results\n")
            printf("-> flux recovery: %.2g\n" % (
                np.nanmean(raw[TableColumn.FLUX] /
                           raw[TableColumn.FLUX_DET])))

        results: HDUList
        filter_string: str | None = star_bug_base.filter
        assert filter_string is not None
        assert raw is not None
        if (results := compile_results(
                raw, image=star_bug_base.main_image.data,
                filter_string=filter_string,
                plot_ast=config.ast_plot_filename)):
            out_dir: str
            b_name: str
            out_dir, b_name, _ = StarbugBase.sort_output_names(
                f_name, param_output=config.output_file)
            if config.verbose_logs:
                printf("--> %s/%s-ast.fits\n" % (out_dir, b_name))
            results.writeto("%s/%s-ast.fits" % (out_dir, b_name),
                            overwrite=True)

            # autosave clean-up
            # noinspection SpellCheckingInspection
            for _f_name in glob.glob("sbast-autosave*.tmp"):
                _f_name: str
                os.remove(_f_name)

        else:
            p_error("results compilation failed\n")

    else:
        p_error("must include a fits image to work on\n")
        exit_code = ExitStates.EXIT_FAIL

    # Wrapped fix to handle rapid multiprocess teardowns safely
    try:
        share_memory.unlink()
    except FileNotFoundError:
        # The memory handle was already unlinked safely by another thread
        pass
    return exit_code


def ast_main_entry() -> ExitStates:
    """Command line entry point"""
    # globals
    c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
    share_memory: SharedMemory = (
        shared_memory.SharedMemory(create=True, size=c.nbytes))
    loading_buffer: np.ndarray = np.ndarray(
        c.shape, dtype=c.dtype, buffer=share_memory.buf)
    return ast_main(sys.argv, share_memory, loading_buffer)
