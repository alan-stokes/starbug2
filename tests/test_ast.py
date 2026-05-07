import os, glob
import pytest
from starbug2.bin.ast import ast_main
from starbug2.constants import EXIT_SUCCESS

run = lambda s:ast_main(s.split())

def test_run_basic():
    clean()
    assert run("starbug2-ast -N10 -S10 tests/dat/image.fits") == EXIT_SUCCESS
    assert run("starbug2-ast -N30 -S10 -n3 tests/dat/image.fits") == EXIT_SUCCESS
    assert run("starbug2-ast -N30 -S10 -n3 -o/tmp/ tests/dat/image.fits") == EXIT_SUCCESS
    clean()

def test_run_harsh_inputs():
    clean()
    assert run("starbug2-ast -N1 -S1000 tests/dat/image.fits") == EXIT_SUCCESS
    assert run("starbug2-ast -N1000 -S1 tests/dat/image.fits") == EXIT_SUCCESS
    assert run("starbug2-ast -N10 -S10 -n100 tests/dat/image.fits") == EXIT_SUCCESS
    assert run("starbug2-ast -N1000 -S1000 -n1000 tests/dat/image.fits") == EXIT_SUCCESS

def clean():
    files=glob.glob("tests/dat/*")
    files.remove("tests/dat/image.fits")
    files.remove("tests/dat/psf.fits")
    for file_name in files:
        os.remove(file_name)
    if os.path.exists("starbug.param"): os.remove("starbug.param")

