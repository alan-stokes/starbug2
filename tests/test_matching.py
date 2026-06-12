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

import os,numpy as np
from typing import Final

import pytest

from starbug2.constants import (
    CAT_NUM, E_FLUX, RA, DEC, FLUX, FILTER, NUM, FLAG, MAG)
from starbug2.matching.band_match import BandMatch
from starbug2.matching.cascade_match import CascadeMatch
from starbug2.matching.exact_value_match import ExactValueMatch
from starbug2.matching.generic_match import GenericMatch
from starbug2.star_bug_config import StarBugMainConfig
from starbug2.utils import import_table, fill_nan
from starbug2.bin.main import starbug_main
from astropy.table import Table
from astropy import units
from astropy.units import Quantity

from tests.generic import (
    TEST_IMAGE_FITS, check_shape, clean, TEST_FILTER_STRING, TEST_PATH_STR)

IMAGE_2_FITS: Final[str] = os.path.join(TEST_PATH_STR, "image2.fits")
IMAGE_AP_FITS: Final[str] = os.path.join(TEST_PATH_STR, "image-ap.fits")
IMAGE_2_AP_FITS: Final[str] = os.path.join(TEST_PATH_STR, "image2-ap.fits")



@pytest.fixture(autouse=True)
def init():

    clean()
    # noinspection SpellCheckingInspection
    starbug_main(
        f"starbug2 -Ds SIGSRC=10 {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}".split())
    # noinspection SpellCheckingInspection
    starbug_main(
        "starbug2 -Ds SIGSRC=3 -o "
        f"{IMAGE_2_FITS} {TEST_IMAGE_FITS} {TEST_FILTER_STRING}".split())
    starbug_main(
        f"starbug2 -d {IMAGE_AP_FITS} --background {TEST_IMAGE_FITS}"
        f" {TEST_FILTER_STRING}".split())
    os.system(f"cp {TEST_IMAGE_FITS} {IMAGE_2_FITS}")
    starbug_main(
        f"starbug2 -d {IMAGE_2_AP_FITS} --background {IMAGE_2_FITS}"
        f" {TEST_FILTER_STRING}".split())

def cats():
    t1 = [[ 0.0, 0.0, 1.0, 0.1],
          [ 0.1, 0.1, 1.0, 0.1],
          [ 0.2, 0.1, 2.0, 0.2],
          [ 0.1, 0.2, 2.0, 0.2],
          [ 1.0, 1.0, 100, 0.1],
          [ 1.1, 1.0, 100, 0.1],
        ]

    t2 = [[ 0.0, 0.0, 1.1, 0.1],
          [ 0.2, 0.1, 2.1, 0.2],
          [ 0.1, 0.2, 2.1, 0.2],
          [ 1.0, 1.0, 101, 0.1],
          [ 1.1, 1.0, 101, 0.1],
          [ 2.0, 2.0, 201, 0.1],
        ]

    cat1 = Table(
        np.array(t1), names=[RA, DEC, FLUX, E_FLUX], meta={FILTER:'a'})
    cat2 = Table(
        np.array(t2), names=[RA, DEC, FLUX, E_FLUX], meta={FILTER:'b'})
    return [cat1, cat2]

            

class TestGenericMatch:

    def test_initialing(self):
        config = StarBugMainConfig()

        m = GenericMatch()
        assert m.col_names is None
        assert not m.filter 
        assert m.threshold is None
        assert m.verbose == config.verbose_logs

        m = GenericMatch(
            filter_string=MAG, col_names=[RA],
            threshold=config.match_threshold_arc_sec_as_an_arc_sec,
            verbose=True)
        assert m.col_names == [RA]
        assert m.filter == MAG
        assert (m.threshold.value ==
                config.match_threshold_arc_sec_as_an_arc_sec.value)
        assert m.verbose == True

        assert isinstance(m.__str__(), str)

    def test_generic_match1(self):
        categories: list[Table | None] = [import_table(f) for f in (
            f"{IMAGE_AP_FITS}", f"{IMAGE_2_AP_FITS}")]
        assert categories is not None
        category1: Table | None = categories[0]
        category2: Table | None = categories[1]
        if category1 is None or category2 is None:
            raise Exception("failed to import tables")

        config = StarBugMainConfig()
        m : GenericMatch = GenericMatch(
            threshold=config.match_threshold_arc_sec_as_an_arc_sec)

        out: Table = m(categories)
        assert isinstance(out, Table)
        name : str
        for name in category1.colnames:
            if name != CAT_NUM:
                assert "%s_1" % name in out.colnames
                assert "%s_2" % name in out.colnames
        assert len(out) >= len(category1)
        assert len(out) >= len(category2)
        assert m.filter == "F444W"

        out = m(categories, join_type="and")
        print(out)
        assert len(out) <= len(category1)
        assert len(out) <= len(category2)

    def test_generic_match2(self):
        categories = [import_table(f) for f in (
            f"{IMAGE_AP_FITS}",
            f"{IMAGE_2_AP_FITS}")]
        config = StarBugMainConfig()
        m = GenericMatch(
            col_names=[RA],
            threshold=config.match_threshold_arc_sec_as_an_arc_sec)
        out = m(categories)

        assert out.colnames == ["RA_1", "RA_2"]

    def test_finish_matching(self):
        categories: list[Table | None] = [import_table(f) for f in (
            f"{IMAGE_AP_FITS}",
            f"{IMAGE_2_AP_FITS}")]
        config = StarBugMainConfig()
        m: GenericMatch = GenericMatch(
            threshold=config.match_threshold_arc_sec_as_an_arc_sec
        )
        out: Table = m(categories)
        m.finish_matching(out)

        category1: Table | None = categories[0]
        category2: Table | None = categories[1]
        if category1 is None or category2 is None:
            raise Exception("failed to import tables")

        m: GenericMatch = GenericMatch(
            col_names=[RA, DEC, FLUX], filter_string="F444W",
            threshold=config.match_threshold_arc_sec_as_an_arc_sec)
        av: Table = m.finish_matching(m.match([category1, category2]))
        # noinspection SpellCheckingInspection
        assert av.colnames == [
            "RA", "DEC", "flux", "stdflux", "flag", "F444W", "eF444W", "NUM"]

        m: GenericMatch = GenericMatch(
            col_names=[RA, DEC, FLUX],
            threshold=config.match_threshold_arc_sec_as_an_arc_sec)
        c: Table
        for c in categories:
            del c.meta[FILTER]
        av: Table = m.finish_matching(m.match([category1, category2]))
        # noinspection SpellCheckingInspection
        assert av.colnames == [
            "RA", "DEC", "flux", "stdflux", "flag", "MAG", "eMAG", "NUM"]


    def test_vals(self):
        config = StarBugMainConfig()
        m = GenericMatch(
            threshold=config.match_threshold_arc_sec_as_an_arc_sec)
        out = m(cats())
        t = [[ 0.0, 0.0, 1.0, 0.1,   0.0, 0.0, 1.1, 0.1],
             [ 0.1, 0.1, 1.0, 0.1,   np.nan, np.nan, np.nan, np.nan],
             [ 0.2, 0.1, 2.0, 0.2,   0.2, 0.1, 2.1, 0.2],
             [ 0.1, 0.2, 2.0, 0.2,   0.1, 0.2, 2.1, 0.2],
             [ 1.0, 1.0, 100, 0.1,   1.0, 1.0, 101, 0.1],
             [ 1.1, 1.0, 100, 0.1,   1.1, 1.0, 101, 0.1],
             [np.nan, np.nan, np.nan, np.nan, 2.0, 2.0, 201, 0.1] ]
        c = Table(
            np.array(t),
            names=[
                "RA_1", "DEC_1", "flux_1", "eflux_1", "RA_2", "DEC_2",
                "flux_2", "eflux_2"])
        check_shape(c, out)


class TestCascade:
    def test_cascade_match(self):
        [import_table(f) for f in (f"{IMAGE_AP_FITS}", f"{IMAGE_2_AP_FITS}")]
        config = StarBugMainConfig()
        CascadeMatch(threshold=config.match_threshold_arc_sec_as_an_arc_sec)
        
        
    def test_vals(self):
        t = [[ 0.0, 0.0, 1.0, 0.1,   0.0, 0.0, 1.1, 0.1],
             [ 0.1, 0.1, 1.0, 0.1,   np.nan, np.nan, np.nan, np.nan],
             [ 0.2, 0.1, 2.0, 0.2,   0.2, 0.1, 2.1, 0.2],
             [ 0.1, 0.2, 2.0, 0.2,   0.1, 0.2, 2.1, 0.2],
             [ 1.0, 1.0, 100, 0.1,   1.0, 1.0, 101, 0.1],
             [ 1.1, 1.0, 100, 0.1,   1.1, 1.0, 101, 0.1],
             [2.0, 2.0, 201, 0.1,    np.nan, np.nan, np.nan, np.nan] ]
        c = Table(
            np.array(t),
            names=[ "RA_1", "DEC_1", "flux_1", "eflux_1", "RA_2", "DEC_2",
                    "flux_2", "eflux_2"])
        m = CascadeMatch(threshold=Quantity(2, unit=units.arcsec))
        out = m.match(cats())
        print(out)

        check_shape(c, out)


class TestBandMatch:
    def test_init(self):
        filters = ["a", "b", "c"]
        m = BandMatch(fltr=filters)
        assert m.filter_list == ["a", "b", "c"]


    def test_order_catalogue_jwst_meta(self):
        a = Table(None, meta={"FILTER": 'F115W'})
        b = Table(None, meta={"FILTER": 'F187N'})
        c = Table(None, meta={"FILTER": 'F770W'})

        m = BandMatch()
        assert m.filter is None
        assert m.order_catalogues( [a,c,b] ) == [a,b,c]
        assert m.filter_list == ["F115W", "F187N", "F770W"]


    def test_order_catalogue_jwst_col_names(self):
        a = Table(None, names=['F115W'])
        b = Table(None, names=['F187N'])
        c = Table(None, names=['F770W'])

        m = BandMatch()
        assert m.filter is None
        assert m.order_catalogues( [a, c, b] ) == [a, b, c]
        assert m.filter_list == ["F115W", "F187N", "F770W"]


    def test_order_catalogue_filter_meta(self):
        a = Table(None, meta={"FILTER":'a'})
        b = Table(None, meta={"FILTER":'b'})
        c = Table(None, meta={"FILTER":'c'})

        m = BandMatch(fltr=["a", "b", "c"])
        assert m.order_catalogues( [a,c,b] ) == [a,b,c]


    def test_order_catalogue_filter_col_names(self):
        a = Table(None, names=['a'])
        b = Table(None, names=['b'])
        c = Table(None, names=['c'])

        m = BandMatch(fltr=["a", "b", "c"])
        assert m.order_catalogues( [a, c, b] ) == [a, b, c]


    def test_match(self):
        t1 = [[1., 1., 1, 1,0],
              [2., 2., 2, 2, 0],
              [3., 3., 3, 3, 0],
             ]
        t2 = [[1., 1., 1, 1, 0],
              [2., 2., 2, 2, 0],
              [4., 4., 4, 4, 1],
             ]
        t3 = [[1., 1., 1, 1, 0],
              [4., 4., 4, 4, 2],
             ]

        f = float
        categories = [
            Table(np.array(t1), names=[RA, DEC, "A", NUM, FLAG],
                  dtype=[f, f, f, f, np.uint16], meta={FILTER: "A"}),
            Table(np.array(t2), names=[RA, DEC, "B", NUM, FLAG],
                  dtype=[f, f, f, f, np.uint16], meta={FILTER: "B"}),
            Table(np.array(t3), names=[RA, DEC, "C", NUM, FLAG],
                  dtype=[f, f, f, f, np.uint16], meta={FILTER: "C"})]

        bm = BandMatch(fltr=["A", "B", "C"],
                       threshold=[0.1 * units.arcsec, 0.2 * units.arcsec])
        res = bm(categories)
        print(res)
        assert res.colnames == [RA, DEC, NUM, FLAG, "A", "B", "C"]
        bm(categories, method="bootstrap")


def test_parse_mask():
    import_table(f"{IMAGE_AP_FITS}")


def test_match_with_masks():
    t1 = [[0, 0, 1],
          [1, 1, 1],
          [2, 2, 1],
          [3, 3, 1]]
    t2 = [[0, 0, 1],
          [1, 1, 1],
          [2, 2, 0],
          [3, 3, 1]]
    t3 = [[0, 0, 1],
          [1, 1, 1],
          [2, 2, 0],
          [3, 3, 1]]
    cat1 = Table(np.array(t1, float), names=[RA, DEC, "a"])
    cat2 = Table(np.array(t2, float), names=[RA, DEC, "a"])
    cat3 = Table(np.array(t3, float), names=[RA, DEC, "a"])
    mask = [
        np.array([True, True, False, True]),
        None,
        np.array([True, True, True, False])]
    config = StarBugMainConfig()
    res = GenericMatch(
        threshold=config.match_threshold_arc_sec_as_an_arc_sec
    ).match([cat1, cat2, cat3], mask=mask)
    print(res)


def test_exact_match():
    cat1 = Table(
        np.array([["CN1", 1], ["CN2", 1], ["CN3", 1]]),
        names=["CN", "i"], dtype=(str, int))
    cat2 = Table(
        np.array([["CN2", 2], ["CN3", 2], ["CN4", 2]]),
        names=["CN", "i"], dtype=(str, int))
    cat3 = Table(
        np.array([["CN3", 3], ["CN4", 3], ["CN5", 3]]),
        names=["CN", "i"], dtype=(str, int))
    
    arr = [[1, np.nan, np.nan],
           [1,      2, np.nan],
           [1,      2,      3],
           [np.nan, 2,      3],
           [np.nan, np.nan, 3]]
    correct = Table(
        np.array(arr), dtype=[int, int, int], names=["i_1", "i_2", "i_3"], 
        masked=True)
    for i, col in enumerate(correct.columns.values()):
        col.mask = np.isnan(np.array(arr)[:, i])
    correct.add_column(["CN1", "CN2", "CN3", "CN4", "CN5"], name="CN", index=0)
    res = ExactValueMatch(value="CN").match([cat1, cat2, cat3])
    assert all(res==fill_nan(correct)) # type: ignore
        
if __name__=="__main__":
    test_exact_match()
