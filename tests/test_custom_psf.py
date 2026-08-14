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
import pytest

from command_line_interfaces.main import starbug_internal_main
from constants import ExitStates
from core.star_bug_config import StarBugMainConfig
from generic import TEST_PATH_STR, TEST_IMAGE_FITS


def create_config_file(
    config: StarBugMainConfig = StarBugMainConfig()) -> StarBugMainConfig:
    """
    generate the param file used for command line behaviour.
    :param config: the config, or uses a default
    :return: None
    """
    config.unfreeze()
    config.do_custom_psf = True
    config.custom_psf_size_pixels = 51
    config.output_file = TEST_PATH_STR
    config.fits_images = [TEST_IMAGE_FITS]
    config.custom_filter = "F444W"
    config.full_width_half_max = 2
    config.freeze()
    return config


def test_custom_psf() -> None:
    config: StarBugMainConfig = create_config_file()
    exit_code: ExitStates
    exit_code = starbug_internal_main(config)
    assert exit_code == ExitStates.EXIT_SUCCESS


def test_custom_psf_even_fail() -> None:
    """
    ensures the system fails due to even pixels.
    :return: None
    """
    config: StarBugMainConfig = create_config_file()
    config.unfreeze()
    with pytest.raises(Exception):
        config.custom_psf_size_pixels = 50
