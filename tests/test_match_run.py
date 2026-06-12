import os

import pytest

from starbug2.bin.main import starbug_main
from starbug2.bin.match import match_main
from starbug2.constants import EXIT_FAIL, EXIT_EARLY, EXIT_SUCCESS
from tests.generic import (
    clean, TEST_IMAGE_FITS, TEST_PATH, TEST_FILTER_STRING)

run = lambda s:match_main(s.split())

OUT_1_FITS = os.path.join(TEST_PATH, "out1.fits")
OUT_2_FITS = os.path.join(TEST_PATH, "out2.fits")
OUT_1_AP_FITS = os.path.join(TEST_PATH, "out1-ap.fits")
OUT_2_AP_FITS = os.path.join(TEST_PATH, "out2-ap.fits")
IMAGE_AP_FITS = os.path.join(TEST_PATH, "image-ap.fits")

def test_match_start():
    assert run("starbug2-match") == EXIT_FAIL
    assert run("starbug2-match -h") == EXIT_SUCCESS
    assert run("starbug2-match -vh") == EXIT_SUCCESS

def test_match_bad_input():
    assert run("starbug2-match ") == EXIT_FAIL
    assert run(f"starbug2-match {TEST_IMAGE_FITS}") == EXIT_EARLY
    assert run("starbug2-match badinput.fits") == EXIT_FAIL
    assert run("starbug2-match badinput.txt") == EXIT_FAIL
    starbug_main(f"starbug2 -D {TEST_IMAGE_FITS} {TEST_FILTER_STRING}".split())
    assert run(f"starbug2-match {IMAGE_AP_FITS}") == EXIT_EARLY

def test_match_basic_run_through():
    starbug_main(
        f"starbug2 -Do {OUT_1_FITS}"
        f"  {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}".split())
    starbug_main(
        f"starbug2 -Do {OUT_2_FITS} "
        f" {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}".split())
    assert (run(
        f"starbug2-match "
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(
        f"starbug2-match"
        f" {OUT_1_AP_FITS}"
        f" {OUT_2_AP_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)

def test_mask():
    starbug_main(
        f"starbug2 -Do "
        f"{OUT_1_FITS}  {TEST_IMAGE_FITS}"
        f" -s FILTER=F444W".split())
    starbug_main(
        f"starbug2 -Do"
        f"{OUT_2_FITS}  {TEST_IMAGE_FITS}"
        f" -s FILTER=F444W ".split())
    assert run(
        f"starbug2-match -vmF444W>20 "
        f"{OUT_1_AP_FITS} "
        f"{OUT_2_AP_FITS}") == EXIT_SUCCESS



@pytest.fixture(autouse=True)
def init():
    clean()


