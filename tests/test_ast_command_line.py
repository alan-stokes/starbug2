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
from astropy.table import Table
import subprocess
import os

import tests.generic as generic
from starbug2.command_line_interfaces.main import starbug_internal_main
from starbug2.constants import DEFAULT_PARAM_FILE_NAME, ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from tests.generic import TEST_PATH_STR, TEST_IMAGE_FITS, TEST_PSF_FITS


def create_working_param_file(
        config: StarBugMainConfig = StarBugMainConfig()) -> None:
    """
    generate the param file used for command line behaviour.
    :param config: the config, or uses a default
    :return: None
    """
    config.unfreeze()
    config.custom_filter = "F444W"
    config.do_artificial_star_test = True
    config.do_artificial_star_test_results = True
    config.generate_local_param_file = True
    config.output_file = TEST_PATH_STR
    config.psf_file_override = TEST_PSF_FITS
    config.freeze()
    starbug_internal_main(config)


def test_artificial_star_command_line() -> None:
    """
    basic test of ast.
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()
    create_working_param_file()
    result = subprocess.run(
        ["starbug2-ast",
         f"-p{os.path.join(TEST_PATH_STR, DEFAULT_PARAM_FILE_NAME)}",
         f"{TEST_IMAGE_FITS}"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == ExitStates.EXIT_SUCCESS

    ast_file: str = os.path.join(TEST_PATH_STR, "image-ast.fits")
    assert os.path.exists(ast_file)
    ast_psf_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ast.fits"), hdu=2).copy()
    assert len(ast_psf_file) == 1000

    generic.clean()


def test_artificial_star_command_line_change_n_tests() -> None:
    """
    test change to n tests
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()

    # create the param file for this command line test.
    config: StarBugMainConfig = StarBugMainConfig()
    config.unfreeze()
    config.artificial_star_tests_count = 2
    config.freeze()
    create_working_param_file(config)

    # execute the ast framework
    result = subprocess.run(
        ["starbug2-ast",
         f"-p{os.path.join(TEST_PATH_STR, DEFAULT_PARAM_FILE_NAME)}",
         f"{TEST_IMAGE_FITS}"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == ExitStates.EXIT_SUCCESS

    ast_file: str = os.path.join(TEST_PATH_STR, "image-ast.fits")
    assert os.path.exists(ast_file)
    ast_psf_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ast.fits"), hdu=2).copy()
    assert len(ast_psf_file) == 20
    generic.clean()


def test_artificial_star_command_line_change_n_stars() -> None:
    """
    change of n_stars
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()

    # create the param file for this command line test.
    config: StarBugMainConfig = StarBugMainConfig()
    config.unfreeze()
    config.stars_per_artificial_test = 2
    config.freeze()
    create_working_param_file(config)

    # execute the ast framework
    result = subprocess.run(
        ["starbug2-ast",
         f"-p{os.path.join(TEST_PATH_STR, DEFAULT_PARAM_FILE_NAME)}",
         f"{TEST_IMAGE_FITS}"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == ExitStates.EXIT_SUCCESS

    ast_file: str = os.path.join(TEST_PATH_STR, "image-ast.fits")
    assert os.path.exists(ast_file)
    ast_psf_file: Table = Table.read(
        os.path.join(TEST_PATH_STR, "image-ast.fits"), hdu=2).copy()
    assert len(ast_psf_file) == 200
    generic.clean()


def test_artificial_star_command_line_help() -> None:
    """
    checks help works.
    :return: None
    """
    generic.verify_test_data_exists()
    generic.clean()

    # execute the ast framework
    result = subprocess.run(
        ["starbug2-ast", "-h"],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == ExitStates.EXIT_EARLY

    output = result.stdout + result.stderr

    config: StarBugMainConfig = StarBugMainConfig()
    assert config.generate_ast_help_string() in output

    generic.clean()
