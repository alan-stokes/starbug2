import glob
import os

from starbug2.filters import STAR_BUG_FILTERS
from starbug2.constants import STARBUG_DATA_DIR
from starbug2.initialise_psf_data import init_starbug_for_jwst


def test_init():
    os.environ[STARBUG_DATA_DIR] = "/tmp/starbug"
    d = os.getenv(STARBUG_DATA_DIR)
    init_starbug_for_jwst()

    for f in STAR_BUG_FILTERS:
        assert glob.glob("%s/*%s*" % (d, f))


