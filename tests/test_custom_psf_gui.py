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

import pytest

from generic import TEST_IMAGE_FITS
from starbug2.command_line_interfaces.main import starbug_internal_main
from starbug2.constants import ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from tests.generic import (
    TEST_PATH_STR, clean, verify_test_data_exists, TEST_JWST_FITS)


def create_config_file(
        config: StarBugMainConfig = StarBugMainConfig()) -> StarBugMainConfig:
    """
    generate the param file used for command line behaviour.
    :param config: the config, or uses a default
    :return: None
    """
    config.unfreeze()
    config.do_custom_psf_gui = True
    config.custom_psf_size_pixels = 51
    config.output_file = TEST_PATH_STR
    config.fits_images = [TEST_JWST_FITS]
    config.custom_filter = "F444W"
    config.full_width_half_max = 2
    config.sharp_cutoff_high = 1
    config.sharp_cutoff_low = 0
    config.freeze()
    return config


#@pytest.mark.skipif(
#    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") is None or
#    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") == "false",
#    reason="UI test not been fully implemented"
#)
def test_custom_psf_gui(qtbot) -> None:
    clean()
    verify_test_data_exists()
    config: StarBugMainConfig = create_config_file()
    exit_code: ExitStates

    exit_code = starbug_internal_main(config)
    qtbot.stopForInteraction()
    assert exit_code == ExitStates.EXIT_SUCCESS

    # verify files were made as expected.
    custom_stars_file: str = os.path.join(
        TEST_PATH_STR, "image_custom_fit_stars-ap.fits")
    custom_c_psf_file: str = os.path.join(
        TEST_PATH_STR, "image_custom-c-psf.fits")

    assert os.path.exists(custom_stars_file)
    assert os.path.exists(custom_c_psf_file)

    clean()
