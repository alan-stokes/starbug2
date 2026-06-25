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
import glob
from urllib import request
from typing import Final

import numpy as np
import pytest
from astropy.io import fits

from starbug2.constants import (
    STAR_BUG_TEST_DAT_ENV, ImageHeaderTags, MIRI_STRING, MIRI_IMAGE)
from starbug2.star_bug_config import StarBugMainConfig

# paths to test files
TEST_PATH: Final[str | None] = os.getenv(STAR_BUG_TEST_DAT_ENV)
if TEST_PATH is None:
    raise Exception("cant find the test data environmental variable")
TEST_PATH_STR: Final[str] = str(TEST_PATH)
TEST_IMAGE_FITS: Final[str] = os.path.join(TEST_PATH, "image.fits")
TEST_PSF_FITS: Final[str] = os.path.join(TEST_PATH, "psf.fits")
TEST_NGC_FITS: Final[str] = os.path.join(TEST_PATH, "ngc6822_F770W_i2d.fits")
TEST_README: Final[str] = os.path.join(TEST_PATH, "readme.txt")
TEST_BLANK: Final[str] = str(os.path.join(str(TEST_PATH), "blank.fits"))
TEST_AST_FILLED: Final[str] = str(
    os.path.join(str(TEST_PATH), "inserted_image_for_test_1.fits"))
TEST_SEED = 42

# the filter string for tests to ensure they all use the same stuff
TEST_FILTER_STRING = "-s FILTER=F444W -G"
GITHUB_RELEASE_URL = (
    "https://github.com/alan-stokes/starbug2/releases/download/TEST_DATA/")
REQUIRED_FILES = ["image.fits", "psf.fits", "ngc6822_F770W_i2d.fits"]

def verify_test_data_exists() -> None:
    # Check if the specific FITS file is missing
    if not os.path.exists(TEST_IMAGE_FITS):
        print(
            f"\n⚠️ Test file missing due to merge. "
            f"Downloading all from GitHub Releases...")

        for filename in REQUIRED_FILES:
            file_path = os.path.join(str(TEST_PATH), filename)
            url = f"{GITHUB_RELEASE_URL}/{filename}"
            if not os.path.exists(file_path):
                try:
                    request.urlretrieve(url, file_path)
                except Exception as e:
                    pytest.fail(
                        f"Failed to download test asset from GitHub "
                        f"Release: {e}")


def create_default_config() -> StarBugMainConfig:
    """
    creates a default config where everything points to the test output dir
    :return: a config
    :rtype StarBugMainConfig
    """
    config: StarBugMainConfig = StarBugMainConfig()
    config.output_file = TEST_PATH
    return config


def clean() -> None:
    """
    cleans up the test data folder for new tests.
    :return: None
    """
    files = glob.glob(os.path.join(str(TEST_PATH), "*"))
    files.remove(TEST_IMAGE_FITS)
    files.remove(TEST_PSF_FITS)
    files.remove(TEST_NGC_FITS)
    files.remove(TEST_README)
    for file_name in files:
        os.remove(file_name)
    if os.path.exists("starbug.param"):
        os.remove("starbug.param")


def check_shape(c, out) -> None:
    """
    checks shape.
    :param c: array 1
    :param out: array 2
    :return: None
    """
    assert np.shape(c) == np.shape(out)
    for m in range(len(c)):
        for n in range(len(c[m])):
            a = c[m][n]
            b = out[m][n]
            assert np.isnan(a) == np.isnan(b)
            if not np.isnan(a) or not np.isnan(b):
                assert a == b


def create_blank_fits(size=(2048, 2048)):
    """
    creates a blank fits file.
    :param size: the size of the fits file.
    :return: None
    """
    print(f"Generating blank space image of size {size[0]}x{size[1]}...")

    # Create a 2D numpy array of zeros (using float32 for standard precision)
    blank_data = np.zeros(size, dtype=np.float32)

    # Create background noise: mean of 10.0 counts, standard deviation of 1.0
    rng = np.random.default_rng(seed=TEST_SEED)
    background_noise = rng.normal(loc=10.0, scale=1.0, size=blank_data.shape)

    # Wrap the data inside a Primary HDU
    primary_hdu = fits.PrimaryHDU(data=background_noise)

    # Add essential metadata headers so pipeline loaders don't choke
    header = primary_hdu.header
    header["EXTNAME"] = "PRIMARY"
    header["OBJECT"] = "BLANK_SPACE_CI"
    header["COMMENT"] = (
        "Artificial black space for starbug2 integration tests.")
    header[ImageHeaderTags.DETECTOR] = MIRI_IMAGE
    header[ImageHeaderTags.INSTRUMENT] = MIRI_STRING

    # Write the file out to disk
    # overwrite=True ensures test scripts can recreate this file on
    # every run
    primary_hdu.writeto(TEST_BLANK, overwrite=True)
    print(f"✅ Successfully saved to {TEST_BLANK}")
