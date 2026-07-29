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
from astropy.table import Table

from core.main_components.one_time_runs import starbug_one_time_runs
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from tests.generic import TEST_PATH_STR, TEST_PSF_FITS
from starbug2.constants import TableColumn
from tests import generic


def update_config_for_fake_stars_into_blank(config: StarBugMainConfig) -> None:
    config.unfreeze()
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.fits_images = [generic.TEST_BLANK]
    config.psf_file_override = TEST_PSF_FITS
    config.do_star_detection = False
    config.do_aperture_photometry = False
    config.do_artificial_star_test = True
    config.ast_load_psf = True
    config.ast_seed = 42
    config.test_magnitude_bright_limit = 20
    config.test_magnitude_faint_limit = 22
    config.stars_per_artificial_test = 1
    config.artificial_star_tests_count = 1
    config.ast_save_added_image = True
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR
    config.freeze()


def test_ast_output_data():
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits()
    update_config_for_fake_stars_into_blank(config)
    config.unfreeze()
    config.verbose_logs = True
    config.sigma_sky = 4
    config.sigma_source = 10
    config.freeze()

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)

    # execute add stars and do test
    entrance.run_starbug()

    artificial_stars_detections: Table = entrance.detections
    fake_star_locations: Table = entrance.ast_star_source_list

    assert len(artificial_stars_detections) == 1
    assert len(fake_star_locations) == 1

    assert (artificial_stars_detections[0][TableColumn.X_CENTROID] ==
            pytest.approx(fake_star_locations[0][TableColumn.X_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.Y_CENTROID] ==
            pytest.approx(fake_star_locations[0][TableColumn.Y_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.FLUX] ==
            pytest.approx(fake_star_locations[0][TableColumn.FLUX], abs=0.1))

    # execute output generation
    config = StarBugMainConfig()
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.fits_images = [generic.TEST_BLANK]
    config.psf_file_override = TEST_PSF_FITS
    config.do_artificial_star_test_results = True
    config.plot_ast = os.path.join(TEST_PATH_STR, "plot")
    config.ast_plot_filename = os.path.join(TEST_PATH_STR, "plot")
    config.ast_out_tables = [entrance.ast_test_results]
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)
    entrance.run_starbug()

    # check output generated.
    output_file: str = os.path.join(TEST_PATH_STR, "blank-ast.fits")
    output_file2: str = os.path.join(TEST_PATH_STR, "plot.png")
    assert os.path.exists(output_file)
    assert os.path.exists(output_file2)

    config.generate_local_param_file = True
    starbug_one_time_runs(config)

    # clean setup
    #generic.clean()

@pytest.mark.skip
def test_ast_output_psf_photo_data():
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits()
    update_config_for_fake_stars_into_blank(config)
    config.unfreeze()
    config.verbose_logs = True
    config.sigma_sky = 4
    config.sigma_source = 10
    config.ast_no_background = False
    config.ast_no_psf_phot = False
    config.freeze()

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)

    # execute add stars and do test
    entrance.run_starbug()

    artificial_stars_detections: Table = entrance.detections
    fake_star_locations: Table = entrance.ast_star_source_list

    assert len(artificial_stars_detections) == 1
    assert len(fake_star_locations) == 1

    assert (artificial_stars_detections[0][TableColumn.X_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.X_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.Y_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.Y_0], abs=0.5))
    assert (artificial_stars_detections[0][TableColumn.FLUX_DET] ==
            pytest.approx(fake_star_locations[0][TableColumn.FLUX], abs=0.1))

    # execute output generation
    config = StarBugMainConfig()
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.fits_images = [generic.TEST_BLANK]
    config.psf_file_override = TEST_PSF_FITS
    config.do_artificial_star_test_results = True
    config.plot_ast = os.path.join(TEST_PATH_STR, "plot")
    config.ast_plot_filename = os.path.join(TEST_PATH_STR, "plot")
    config.ast_out_tables = [entrance.ast_test_results]
    config.ast_save_added_image_path = TEST_PATH_STR
    config.output_file = TEST_PATH_STR

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)
    entrance.run_starbug()

    # check output generated.
    output_file: str = os.path.join(TEST_PATH_STR, "blank-ast.fits")
    output_file2: str = os.path.join(TEST_PATH_STR, "plot.png")
    assert os.path.exists(output_file)
    assert os.path.exists(output_file2)

    config.generate_local_param_file = True
    starbug_one_time_runs(config)

    # clean setup
    #generic.clean()
