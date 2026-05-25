import os
from starbug2.bin.main import starbug_main
from starbug2.constants import EXIT_EARLY, EXIT_FAIL, EXIT_SUCCESS, EXIT_MIXED
from tests.generic import (
    clean, TEST_IMAGE_FITS, TEST_PATH, TEST_FILTER_STRING)

run = lambda s:starbug_main(s.split())

# different fit files paths
TEST_IMAGE_AP_FITS = os.path.join(TEST_PATH, "image-ap.fits")
TEST_PSF_FITS =  os.path.join(TEST_PATH, "psf.fits")
TEST_IMAGE_BGD_FITS = os.path.join(TEST_PATH, "image-bgd.fits")
TEST_IMAGE_RES_FIT = os.path.join(TEST_PATH, "image-res.fits")
TEST_IMAGE_2_FITS = os.path.join(TEST_PATH, "image2.fits")

def test_start():
    clean()
    assert run("starbug2 -h") == EXIT_EARLY
    assert run("starbug2 -vh") == EXIT_EARLY
    assert run("starbug2 --version") == EXIT_EARLY
    assert run("starbug2 -vDABPh") == EXIT_EARLY
    assert run("starbug2") == EXIT_FAIL
    clean()

def test_param():
    clean()
    assert run("starbug2 --local-param") == EXIT_EARLY
    assert run("starbug2 --update-param") == EXIT_EARLY
    assert (run(
        f"starbug2 -p starbug.param {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    clean()

def test_detect():
    clean()
    assert run(
        f"starbug2 -v {TEST_IMAGE_FITS} {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2 -D {TEST_IMAGE_FITS} {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2 --detect {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert (run(f"starbug2 -D -sSIGSKY=3 -sSIGSRC=15 {TEST_IMAGE_FITS}"
                f" {TEST_FILTER_STRING}") ==
            EXIT_SUCCESS)
    clean()

def test_bgd():
    clean()
    assert run(
        f"starbug2 -D {TEST_IMAGE_FITS} {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert (run(f"starbug2 -d {TEST_IMAGE_AP_FITS}"
                f" -B {TEST_IMAGE_FITS} {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert run(f"starbug2 -d {TEST_IMAGE_AP_FITS} "
               f"--background {TEST_IMAGE_FITS} "
               f"{TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -vf -B {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    clean()

def test_psf():
    clean()
    assert run(f"starbug2 -DB {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -d {TEST_IMAGE_AP_FITS} -b "
               f"{TEST_IMAGE_BGD_FITS} -P {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -fP {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -d {TEST_IMAGE_AP_FITS}"
               f" -P {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -fBP {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -fPs GEN_RESIDUAL=1 {TEST_IMAGE_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    clean()

def test_residual():
    clean()
    assert run(f"starbug2 -DB {TEST_IMAGE_FITS} "
               f"{TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2 -fSs GEN_RESIDUAL=1 {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 {TEST_IMAGE_RES_FIT}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -D {TEST_IMAGE_RES_FIT}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -fB {TEST_IMAGE_RES_FIT}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -fP {TEST_IMAGE_RES_FIT} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -fPs GEN_RESIDUAL=1 {TEST_IMAGE_RES_FIT}"
               f" -sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS

    assert run(f"starbug2 -fSA {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    clean()


def test_n_cores():
    clean()
    os.system(f"cp {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}")
    assert run(
        f"starbug2 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
        f" {TEST_FILTER_STRING}")==EXIT_SUCCESS
    assert run(
        f"starbug2 -n2 {TEST_IMAGE_FITS} "
        f"{TEST_IMAGE_2_FITS} {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert (run(f"starbug2 -vD {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(f"starbug2 -Dn0 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(f"starbug2 -Dn1 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(f"starbug2 -Dn2 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(f"starbug2 -Dn4 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)

    assert run(f"starbug2 -DBP {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS} "
               f"-sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(f"starbug2 -vDBPn2 {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
               f" -sPSF_FILE={TEST_PSF_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert (run(f"starbug2 -DM {TEST_IMAGE_FITS} {TEST_IMAGE_2_FITS}"
                f" {TEST_FILTER_STRING}") == EXIT_SUCCESS)
    assert (run(f"starbug2 -DMn2 {TEST_IMAGE_FITS} "
                f"{TEST_IMAGE_2_FITS} {TEST_FILTER_STRING}") == EXIT_SUCCESS)

    assert run(f"starbug2 -D {TEST_IMAGE_AP_FITS} "
               f"{TEST_IMAGE_FITS} {TEST_FILTER_STRING}") == EXIT_MIXED
    assert run(f"starbug2 -D bad.fits {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING}") == EXIT_MIXED
    clean()

