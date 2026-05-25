from starbug2.bin.ast import ast_main
from starbug2.constants import EXIT_SUCCESS
from tests.generic import TEST_IMAGE_FITS, clean

# main ast run
run = lambda s: ast_main(s.split() + [TEST_IMAGE_FITS])
TEST_FILTER_STRING = "-s FILTER=F444W"

def test_run_basic():
    clean()
    assert run(f"starbug2-ast -N10 -S10 {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2-ast -N30 -S10 -n3 {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2-ast -N30 -S10 -n3 -o /tmp/"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    clean()

def test_run_harsh_inputs():
    clean()
    assert run(
        f"starbug2-ast -N1 -S1000 {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2-ast -N1000 -S1 {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2-ast -N10 -S10 -n100 {TEST_FILTER_STRING}") == EXIT_SUCCESS
    assert run(
        f"starbug2-ast -N1000 -S1000 -n1000"
        f" {TEST_FILTER_STRING}") == EXIT_SUCCESS
    clean()