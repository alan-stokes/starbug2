from starbug2.bin.ast import ast_main
from starbug2.constants import EXIT_SUCCESS, EXIT_FAIL

run = lambda s : ast_main(["starbug2-afs"] + s.split())

def _test_run():
    assert run("tests/dat/image.fits") == EXIT_SUCCESS
    assert run("nope") == EXIT_FAIL
    
