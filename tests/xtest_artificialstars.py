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

from starbug2.bin.ast import ast_main
from starbug2.constants import ExitStates
from tests.generic import TEST_IMAGE_FITS

run = lambda s : (
    ast_main(["starbug2-afs"] + s.split()))

def _test_run():
    assert run(TEST_IMAGE_FITS) == ExitStates.EXIT_SUCCESS
    assert run("nope") == ExitStates.EXIT_FAIL
    
