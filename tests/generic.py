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