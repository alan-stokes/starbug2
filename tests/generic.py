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

import os, glob

import numpy as np

from starbug2.constants import  STAR_BUG_TEST_DAT_ENV

# paths to test files
TEST_PATH = os.getenv(STAR_BUG_TEST_DAT_ENV)
TEST_IMAGE_FITS = os.path.join(TEST_PATH, "image.fits")
TEST_PSF_FITS = os.path.join(TEST_PATH, "psf.fits")

# the filter string for tests to ensure they all use the same stuff
TEST_FILTER_STRING = "-s FILTER=F444W -G"


def clean():
    files = glob.glob(os.path.join(TEST_PATH, "*"))
    files.remove(TEST_IMAGE_FITS)
    files.remove(TEST_PSF_FITS)
    for file_name in files:
        os.remove(file_name)
    if os.path.exists("starbug.param"):
        os.remove("starbug.param")

def check_shape(c, out):
    assert np.shape(c) == np.shape(out)
    for m in range(len(c)):
        for n in range(len(c[m])):
            a = c[m][n]
            b = out[m][n]
            assert np.isnan(a) == np.isnan(b)
            if not np.isnan(a) or not np.isnan(b):
                assert a == b