from starbug2.bin.ast import ast_main
from starbug2.constants import EXIT_SUCCESS
from tests.generic import TEST_IMAGE_FITS, clean, TEST_FILTER_STRING_NO_G

# main ast run
run = lambda s:ast_main(s.split())

def test_run_basic():
    clean()
    assert run(f"starbug2-ast -N10 -S10"
               f" {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING_NO_G}") == EXIT_SUCCESS
    assert run(
        f"starbug2-ast -N30 -S10 -n3"
        f" {TEST_IMAGE_FITS} {TEST_FILTER_STRING_NO_G}") == EXIT_SUCCESS
    assert (run(f"starbug2-ast -N30 -S10 -n3 -o/tmp/ "
                f"{TEST_IMAGE_FITS}"
                f" {TEST_FILTER_STRING_NO_G}") == EXIT_SUCCESS)
    clean()

def test_run_harsh_inputs():
    clean()
    assert run(f"starbug2-ast -N1 -S1000"
               f" {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING_NO_G}") == EXIT_SUCCESS
    assert run(f"starbug2-ast -N1000 -S1"
               f" {TEST_IMAGE_FITS}"
               f" {TEST_FILTER_STRING_NO_G}") == EXIT_SUCCESS
    assert (run(f"starbug2-ast -N10 -S10 -n100"
                f" {TEST_IMAGE_FITS}"
                f" {TEST_FILTER_STRING_NO_G}") == EXIT_SUCCESS)
    assert (run(f"starbug2-ast -N1000 -S1000 -n1000"
                f" {TEST_IMAGE_FITS}"
                f" {TEST_FILTER_STRING_NO_G}") == EXIT_SUCCESS)
