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

import glob
import os

from starbug2.filters import STAR_BUG_FILTERS
from starbug2.constants import STARBUG_DATA_DIR
from starbug2.misc import init_starbug


def xtest_init():
    os.environ[STARBUG_DATA_DIR] = "/tmp/starbug"
    d = os.getenv(STARBUG_DATA_DIR)
    init_starbug()

    for f in STAR_BUG_FILTERS:
        assert glob.glob("%s/*%s*" % (d, f))


