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
from multiprocessing import shared_memory
from multiprocessing.shared_memory import SharedMemory
from typing import Final

import numpy as np
import pytest
from starbug2.command_line_interfaces.ast import ast_main
from starbug2.core.constants import ExitStates
from tests.generic import (
    TEST_IMAGE_FITS, clean, verify_test_data_exists, TEST_PATH)

# main ast run
c: np.ndarray = np.array([0, 0, 0], dtype=np.int64)
share_memory: SharedMemory = (
    shared_memory.SharedMemory(create=True, size=c.nbytes))
loading_buffer: np.ndarray = np.ndarray(
    c.shape, dtype=c.dtype, buffer=share_memory.buf)
TEST_FILTER_STRING: Final[str] = "-s FILTER=F444W"


def run(s):
    return ast_main(
        s.split() + [TEST_IMAGE_FITS], share_memory, loading_buffer
    )


def test_run_basic():
    verify_test_data_exists()
    clean()
    assert (run(
        f"starbug2-ast -N10 -S10 "
        f"--output={TEST_PATH} {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert (run(
        f"starbug2-ast -N30 -S10 -n3 "
        f"--output={TEST_PATH} {TEST_FILTER_STRING}") ==
            ExitStates.EXIT_SUCCESS)
    assert run(
        f"starbug2-ast -N30 -S10 -n3 -o /tmp/"
        f" {TEST_FILTER_STRING}") == ExitStates.EXIT_SUCCESS
    clean()


@pytest.mark.skipif(
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") is None or
    os.getenv("RUN_STAR_BUG_PRODUCTION_TESTS") == "false",
    reason="Harsh stress test locked out of normal development runs due to "
           "length of time to run, CPU resources required which nearly slags"
           " the machine."
)
def test_run_harsh_inputs():
    clean()
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
    clean()


if __name__ == "__main__":
    # This allows you to run the harsh test directly.
    test_run_harsh_inputs()
