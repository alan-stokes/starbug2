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

import numpy as np
import pytest
from astropy.io.fits import ImageHDU, PrimaryHDU, Header
from astropy.nddata import NDData
from astropy.stats import sigma_clipped_stats
from astropy.table import Table

# needed for the clone of the code on the readme of:
# https://photutils.readthedocs.io/en/latest/user_guide/epsf_building.html
from photutils.datasets import load_simulated_hst_star_image
from photutils.datasets import make_noise_image
from photutils.detection import DAOStarFinder
from photutils.psf import (
    extract_stars, EPSFBuilder, EPSFStars, EPSFBuildResult, ImagePSF)

from constants import TableColumn
from custom_psf_gui import common_gui_code
from custom_psf_gui.psf_star_selector import find_stars_to_select
from generic import TEST_JWST_FITS, TEST_JWST_CUSTOM_FILTER
from main_components.custom_psf import CustomPSF
from main_components.one_time_runs import starbug_one_time_runs
from main_components.photometry import Photometry
from starbug2.command_line_interfaces.main import starbug_internal_main
from starbug2.constants import ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug_main import StarbugBase
from tests.generic import (
    TEST_PATH_STR, TEST_IMAGE_FITS, clean, verify_test_data_exists)
from utilities.utils import export_table


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
    config.custom_filter = TEST_JWST_CUSTOM_FILTER
    config.full_width_half_max = 2
    config.freeze()
    return config


def test_custom_psf() -> None:
    clean()
    verify_test_data_exists()
    config: StarBugMainConfig = create_config_file()
    exit_code: ExitStates
    exit_code = starbug_internal_main(config)
    assert exit_code == ExitStates.EXIT_SUCCESS

    # verify files were made as expected.
    custom_stars_file: str = os.path.join(
        TEST_PATH_STR, "image_custom_fit_stars-ap.fits")
    custom_c_psf_file: str = os.path.join(
        TEST_PATH_STR, "image_custom-c-psf.fits")

    assert os.path.exists(custom_stars_file)
    assert os.path.exists(custom_c_psf_file)

    clean()


def run_photutils_selector(
        data: np.ndarray, config: StarBugMainConfig) -> None:
    """
    tests the
    :param data:
    :param config:
    :return:
    """
    starbug_base, exit_state = common_gui_code.detect_stars(config)
    sources_before = starbug_base.detections
    assert sources_before is not None
    (sources, error) = find_stars_to_select(
        data, sources_before, config.psf_generator_stars_to_select,
        config.psf_generator_min_separation,
        config.psf_generator_saturation_limit, config.sharp_cutoff_low,
        config.sharp_cutoff_high, config.psf_generator_grid_bin_x,
        config.psf_generator_grid_bin_y, config.psf_generator_edge_buffer)
    assert sources is not None
    mean_val: float
    median_val: float
    std_val: float
    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.0)
    data -= median_val
    nd_data: NDData = NDData(data=data)

    size: int = 25
    hsize: float = (size - 1) / 2
    x = sources[TableColumn.X_CENTROID]
    y = sources[TableColumn.Y_CENTROID]
    mask: np.ndarray = (
        (x > hsize) & (x < (data.shape[1] - 1 - hsize)) &
        (y > hsize) & (y < (data.shape[0] - 1 - hsize)))
    stars_tbl = Table()
    stars_tbl[TableColumn.X] = x[mask]
    stars_tbl[TableColumn.Y] = y[mask]

    stars: EPSFStars = extract_stars(nd_data, stars_tbl, size=25)
    epsf_builder: EPSFBuilder = EPSFBuilder(
        oversampling=4, maxiters=3, progress_bar=False)
    result: EPSFBuildResult = epsf_builder(stars)
    epsf: ImagePSF = result.epsf
    fitted_stars: EPSFStars = result.fitted_stars

    config: StarBugMainConfig = create_config_file()
    output_dir: str | None = config.output_file
    assert output_dir is not None
    CustomPSF.write_files_to_disk(
        output_dir, epsf, fitted_stars, "plutUtilsTest")


def run_photutils(data: np.ndarray) -> None:
    finder: DAOStarFinder = DAOStarFinder(threshold=100.0, fwhm=1.5)
    sources: Table | None = finder(data)
    assert sources is not None
    mean_val: float
    median_val: float
    std_val: float
    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.0)
    data -= median_val
    nd_data: NDData = NDData(data=data)

    size: int = 25
    hsize: float = (size - 1) / 2
    x = sources[TableColumn.X_CENTROID]
    y = sources[TableColumn.Y_CENTROID]
    mask: np.ndarray = (
        (x > hsize) & (x < (data.shape[1] - 1 - hsize)) &
        (y > hsize) & (y < (data.shape[0] - 1 - hsize)))
    stars_tbl = Table()
    stars_tbl[TableColumn.X] = x[mask]
    stars_tbl[TableColumn.Y] = y[mask]

    stars: EPSFStars = extract_stars(nd_data, stars_tbl, size=25)
    epsf_builder: EPSFBuilder = EPSFBuilder(
        oversampling=4, maxiters=3, progress_bar=False)
    result: EPSFBuildResult = epsf_builder(stars)
    epsf: ImagePSF = result.epsf
    fitted_stars: EPSFStars = result.fitted_stars

    config: StarBugMainConfig = create_config_file()
    output_dir: str | None = config.output_file
    assert output_dir is not None
    CustomPSF.write_files_to_disk(
        output_dir, epsf, fitted_stars, "plutUtilsTest")


def test_custom_psf_using_epsf_building_example() -> None:
    # note this test assumes inspection manually of the output file to verify
    # it looks like a psf. And as of writing the test, its validated that
    # starbugs path does not work properly.
    clean()
    hdu: ImageHDU = load_simulated_hst_star_image()

    data: np.ndarray = hdu.data
    data += make_noise_image(
        data.shape, distribution='gaussian', mean=10.0, stddev=5.0, seed=0)

    run_photutils(data)
    clean()


def test_custom_psf_using_image_fits_and_building_example() -> None:
    # note this test assumes inspection manually of the output file to verify
    # it looks like a psf. And as of writing the test, its validated that
    # starbugs path does not work properly.
    clean()
    verify_test_data_exists()
    config: StarBugMainConfig = create_config_file()
    config.unfreeze()
    config.fits_images = [TEST_IMAGE_FITS]
    config.freeze()
    star_bug_base: StarbugBase | None = StarbugBase(
        TEST_IMAGE_FITS, config=config, ap_file=None,
        bkg_file=None)
    assert star_bug_base is not None
    main_image: ImageHDU | PrimaryHDU = star_bug_base.main_image()
    data: np.ndarray = main_image.data
    run_photutils(data)
    clean()


def test_custom_psf_using_jwst_fits_and_building_example() -> None:
    # note this test assumes inspection manually of the output file to verify
    # it looks like a psf. And as of writing the test, its validated that
    # starbugs path does not work properly.
    clean()
    verify_test_data_exists()
    config: StarBugMainConfig = create_config_file()
    config.unfreeze()
    config.fits_images = [TEST_JWST_FITS]
    config.freeze()
    star_bug_base: StarbugBase | None = StarbugBase(
        TEST_JWST_FITS, config=config, ap_file=None,
        bkg_file=None)
    assert star_bug_base is not None
    main_image: ImageHDU | PrimaryHDU = star_bug_base.main_image()
    data: np.ndarray = main_image.data
    run_photutils(data)
    clean()


def test_jwst_custom_psf_with_selector_from_gui() -> None:
    clean()
    verify_test_data_exists()
    config: StarBugMainConfig = create_config_file()
    config.unfreeze()
    config.fits_images = [TEST_JWST_FITS]
    config.freeze()
    star_bug_base: StarbugBase | None = StarbugBase(
        TEST_JWST_FITS, config=config, ap_file=None,
        bkg_file=None)
    assert star_bug_base is not None
    main_image: ImageHDU | PrimaryHDU = star_bug_base.main_image()
    data: np.ndarray = main_image.data
    run_photutils_selector(data, config)
    clean()


def test_custom_epsf_against_default_epsf() -> None:
    clean()
    verify_test_data_exists()
    config: StarBugMainConfig = create_config_file()

    # set up to run detections and aperture
    config.unfreeze()
    config.do_star_detection = True
    config.do_aperture_photometry = True
    config.fits_images = [TEST_IMAGE_FITS]
    config.freeze()

    # execute starbug to get detections and aperture results
    star_bug_base: StarbugBase = StarbugBase(
        config=config, ap_file=None, bkg_file=None, f_name=TEST_IMAGE_FITS)
    exit_state: ExitStates = star_bug_base.run_starbug(config)
    assert (exit_state == ExitStates.EXIT_SUCCESS)

    detections: Table | None = star_bug_base.detections
    assert detections is not None

    # generate psf aperture with default psf.
    exit_state: ExitStates = star_bug_base.psf_photometry_routine()
    assert (exit_state == ExitStates.EXIT_SUCCESS)
    normal_psf_catalogue: Table | None = star_bug_base.psf_catalogue
    assert normal_psf_catalogue is not None

    # generate custom psf.
    config = create_config_file()
    config.unfreeze()
    config.do_custom_psf = True
    config.ap_file = os.path.join(TEST_PATH_STR, "image-ap.fits")
    config.fits_images = [TEST_IMAGE_FITS]
    config.freeze()
    exit_state: ExitStates = starbug_one_time_runs(config)
    assert (exit_state == ExitStates.EXIT_SUCCESS)

    # verify the custom psf has been generated.
    #psf_file_name = "jw01234-c1003_t005_miri_f770w_i2d-psf.fits"
    psf_file_name = "image_custom-c-psf.fits"
    psf_file_path = os.path.join(TEST_PATH_STR, psf_file_name)
    assert(os.path.exists(psf_file_path))

    # set up to use custom psf and generate psf aperture.
    config = create_config_file()
    config.unfreeze()
    config.fits_images = [TEST_IMAGE_FITS]
    config.psf_file_override = psf_file_path
    config.ap_file = os.path.join(TEST_PATH_STR, "image-ap.fits")
    config.freeze()
    star_bug_base = StarbugBase(
        config=config, ap_file=config.ap_file, bkg_file=None,
        f_name=TEST_IMAGE_FITS)
    star_bug_photometry: Photometry = Photometry()
    (exit_state, custom_psf_catalogue, _) = (
        star_bug_photometry.photometry_routine(
            star_bug_base.filter, star_bug_base.wcs, config,
            star_bug_base.main_image(),
            star_bug_base.log, star_bug_base.image, star_bug_base.info,
            star_bug_base.background,
            star_bug_base.header(), star_bug_base.detections,
            star_bug_base.psf, star_bug_base.full_width_half_max,
            star_bug_base.ap_file, star_bug_base.background_file,
            star_bug_base._out_dir, star_bug_base.b_name))
    assert (exit_state == ExitStates.EXIT_SUCCESS)
    assert custom_psf_catalogue is not None

    # compare results.
    print(f"normal psf catalogue had {len(normal_psf_catalogue)} entries and "
          f"the custom psf catalogue had {len(custom_psf_catalogue)} entries")
    normal_path = os.path.join(TEST_PATH_STR, "normal_catalgogue.fits")
    custom_path = os.path.join(TEST_PATH_STR, "custom_catalgogue.fits")
    export_table(normal_psf_catalogue, normal_path, header=Header())
    export_table(custom_psf_catalogue, custom_path, header=Header())

    # wrap up
    clean()

def test_custom_psf_even_fail() -> None:
    """
    ensures the system fails due to even pixels.
    :return: None
    """
    config: StarBugMainConfig = create_config_file()
    config.unfreeze()
    with pytest.raises(Exception):
        config.custom_psf_size_pixels = 50
