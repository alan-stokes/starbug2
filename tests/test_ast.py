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
import math
import os
from multiprocessing import shared_memory
from multiprocessing.shared_memory import SharedMemory
from typing import Final

import numpy as np
import pytest
from astropy.table import Table
from photutils.psf import ImagePSF

from core.star_bug_config import StarBugMainConfig
from core.starbug_main import StarbugBase
from generic import TEST_PATH_STR, TEST_PSF_FITS
from starbug2.command_line_interfaces.ast import ast_main
from constants import ExitStates, TableColumn
from tests import generic

# main ast run
c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
share_memory: SharedMemory = (
    shared_memory.SharedMemory(create=True, size=c.nbytes))
loading_buffer: np.ndarray = np.ndarray(
    c.shape, dtype=c.dtype, buffer=share_memory.buf)
TEST_FILTER_STRING: Final[str] = "-s FILTER=F444W"


def run(s):
    return ast_main(
        s.split() + [generic.TEST_IMAGE_FITS], share_memory, loading_buffer
    )


def test_run_basic():
    generic.verify_test_data_exists()
    generic.clean()
    # run single core
    assert (run(
        f"starbug2-ast -N10 -S10 -sPSF_FILE={TEST_PSF_FITS} "
        f"--output={generic.TEST_PATH} {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)

    # run multi-core
    assert (run(
        f"starbug2-ast -N30 -S10 -n3 -sPSF_FILE={TEST_PSF_FITS} "
        f"--output={generic.TEST_PATH} {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)

    # run multi-core in /tmp
    assert run(
        f"starbug2-ast -N30 -S10 -n3 -sPSF_FILE={TEST_PSF_FITS} -o /tmp/"
        f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    generic.clean()


@pytest.mark.skipif(
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") is None or
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") == "false",
    reason="Harsh stress test locked out of normal development runs due to "
           "length of time to run, CPU resources required which nearly slags"
           " the machine."
)
def test_run_harsh_inputs():
    generic.clean()
    assert (run(
        f"starbug2-ast -N1 -S1000 {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-ast -N1000 -S1 {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-ast -N10 -S10 -n100 {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert run(
        f"starbug2-ast -N1000 -S1000 -n1000"
        f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    generic.clean()


def test_add_stars_logic():
    generic.verify_test_data_exists()
    generic.clean()

    # create blank fits file.
    config: StarBugMainConfig = StarBugMainConfig()
    generic.create_blank_fits()

    # create config
    config.unfreeze()
    config.ast_add_stars = True
    config.fits_images = [generic.TEST_BLANK]
    config.ast_load_psf = True
    config.psf_file_override = TEST_PSF_FITS
    config.ast_seed = 42
    config.stars_per_artificial_test = 1
    config.ast_save_added_image = True
    config.ast_save_added_image_path = TEST_PATH_STR
    config.custom_filter = generic.TEST_CUSTOM_FILTER
    config.freeze()

    entrance: StarbugBase = StarbugBase(
        config=config, f_name=generic.TEST_BLANK, ap_file=None, bkg_file=None)

    # extract image data before adding stars
    original_image_data: np.ndarray = entrance.main_image().data.copy()

    # execute add stars
    entrance.run_starbug()

    # extract locations and new image data.
    locations: Table = entrance.ast_star_source_list
    image_data: np.ndarray = entrance.main_image().data.copy()

    # set to 0 all the values in both images where these fake star has affected
    # the image.
    x_coord: int = int(locations[0][TableColumn.X_0])
    y_coord: int = int(locations[0][TableColumn.Y_0])
    half_width_x: int = math.ceil(ImagePSF(entrance.psf).origin[0])
    half_width_y: int = math.ceil(ImagePSF(entrance.psf).origin[1])

    # Slice the mask at the star location and flip it to False (
    # do not check these pixels)
    comparison_mask = np.ones_like(original_image_data, dtype=bool)
    comparison_mask[
        y_coord - half_width_y : y_coord + half_width_y,
        x_coord - half_width_x : x_coord + half_width_x
    ] = False

    # test outside fake star location it is identical
    np.testing.assert_array_equal(
        original_image_data[comparison_mask],
        image_data[comparison_mask]
    )

    # check inside the star its not identical
    star_box_original = original_image_data[~comparison_mask]
    star_box_new = image_data[~comparison_mask]
    assert not np.array_equal(star_box_original, star_box_new)


if __name__ == "__main__":
    # This allows you to run the harsh test directly.
    test_run_harsh_inputs()
