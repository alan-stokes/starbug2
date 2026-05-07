import glob
import os

from starbug2 import filters
from starbug2.constants import STARBUG_DATA_DIR
from starbug2.misc import init_starbug


def xtest_init():
    os.environ[STARBUG_DATA_DIR]="/tmp/starbug"
    d=os.getenv(STARBUG_DATA_DIR)
    init_starbug()

    for f in filters:
        assert glob.glob("%s/*%s*"%(d,f))


