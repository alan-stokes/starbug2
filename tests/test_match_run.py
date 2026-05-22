import os
import pytest
from starbug2.bin.main import starbug_main
from starbug2.bin.match import match_main
from starbug2.constants import EXIT_FAIL, EXIT_EARLY, EXIT_SUCCESS
from tests.generic import clean, TEST_IMAGE_FITS, TEST_PATH

run = lambda s:match_main(s.split())

def test_match_start():
    assert run("starbug2-match") == EXIT_FAIL
    assert run("starbug2-match -h") == EXIT_EARLY
    assert run("starbug2-match -vh") == EXIT_EARLY

def test_match_bad_input():
    assert run("starbug2-match ") == EXIT_FAIL
    assert run(f"starbug2-match {TEST_IMAGE_FITS}") == EXIT_EARLY
    assert run("starbug2-match badinput.fits") == EXIT_FAIL
    assert run("starbug2-match badinput.txt") == EXIT_FAIL
    starbug_main(f"starbug2 -D {TEST_IMAGE_FITS}".split())
    assert run(f"starbug2-match {
        os.path.join(TEST_PATH, "image-ap.fits")}") == EXIT_FAIL

def test_match_basic_run_through():
    starbug_main(
        f"starbug2 -Do {os.path.join(TEST_PATH, "out1.fits")}"
        f"  {TEST_IMAGE_FITS}".split())
    starbug_main(
        f"starbug2 -Do {os.path.join(TEST_PATH, "out2.fits")} "
        f" {TEST_IMAGE_FITS}".split())
    assert (run(
        f"starbug2-match {os.path.join(TEST_PATH, "out1-ap.fits")}"
        f" {os.path.join(TEST_PATH, "out2-ap.fits")}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match -G"
        f" {os.path.join(TEST_PATH, "out1-ap.fits")}"
        f" tests/dat/out2-ap.fits") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match -C"
        f" {os.path.join(TEST_PATH, "out1-ap.fits")} "
        f"{os.path.join(TEST_PATH, "out2-ap.fits")}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match -f "
        f"{os.path.join(TEST_PATH, "out1-ap.fits")} "
        f"{os.path.join(TEST_PATH, "out2-ap.fits")}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match -fG "
        f"{os.path.join(TEST_PATH, "out1-ap.fits")}"
        f"{os.path.join(TEST_PATH, "out2-ap.fits")}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match -fC "
        f"{os.path.join(TEST_PATH, "out1-ap.fits")}"
        f"{os.path.join(TEST_PATH, "out2-ap.fits")}") == EXIT_SUCCESS)

def test_mask():
    starbug_main(
        f"starbug2 -Do "
        f"{os.path.join(TEST_PATH, "out1.fits")}  {TEST_IMAGE_FITS}".split())
    starbug_main(
        f"starbug2 -Do "
        f"{os.path.join(TEST_PATH, "out2.fits")}  {TEST_IMAGE_FITS}".split())
    assert run(
        f"starbug2-match -vmF444W>20 "
        f"{os.path.join(TEST_PATH, "out1-ap.fits")}"
        f"{os.path.join(TEST_PATH, "out2-ap.fits")}") == EXIT_SUCCESS



@pytest.fixture(autouse=True)
def init():
    clean()


