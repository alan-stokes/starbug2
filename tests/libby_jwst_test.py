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
from starbug2.constants import ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from tests.generic import (
    TEST_PATH_STR, verify_test_data_exists, clean, TEST_JWST_FITS)


def create_working_config_file(
        config: StarBugMainConfig = StarBugMainConfig()) -> StarBugMainConfig:
    """
    generate the param file used for command line behaviour.
    :param config: the config, or uses a default
    :return: None
    """
    config.unfreeze()
    config.do_artificial_star_test = True
    config.do_artificial_star_test_results = True
    config.artificial_star_tests_count = 1
    config.stars_per_artificial_test = 1
    config.generate_local_param_file = True
    config.output_file = TEST_PATH_STR
    config.fits_images = [TEST_JWST_FITS]
    config.freeze()
    return config


def test_jwst_image_for_filter_processing() -> None:
    """
    test filter detection and ap_corr usage.
    :return: None
    """
    clean()
    verify_test_data_exists()
    config: StarBugMainConfig = create_working_config_file()
    star_bug_base: StarbugBase | None = StarbugBase(
        f_name=TEST_JWST_FITS, config=config, ap_file=None,
        bkg_file=None)
    assert star_bug_base is not None
    exit_state: int = star_bug_base.run_starbug(config)
    assert exit_state == ExitStates.EXIT_SUCCESS
    clean()
