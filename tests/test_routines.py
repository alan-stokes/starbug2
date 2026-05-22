from astropy.io import fits
from astropy.table import Table

from starbug2.constants import X_CENTROID, Y_CENTROID, SCI
from starbug2.routines.background_estimate_routine import (
    BackGroundEstimateRoutine)
from starbug2.routines.detection_routines import DetectionRoutine
from starbug2.routines.source_properties import SourceProperties
from tests.generic import TEST_IMAGE_FITS


class TestDetection:
    im = fits.open(TEST_IMAGE_FITS)[SCI].data
    a = Table([[0, 10], [0, 10]], names=[X_CENTROID, Y_CENTROID])
    b = Table([[20, 10, 50], [20, 10, 0]], names=[X_CENTROID, Y_CENTROID])

    def test_detection_routine_none(self):
        dt = DetectionRoutine()
        assert dt.find_stars(None) is None

    def test_detection_routine_crashes(self):
        dt = DetectionRoutine()
        out = dt.find_stars(self.im.copy())
        assert out is not None

    def test_detection_match(self):
        dt = DetectionRoutine()
        _a = self.a.copy()
        _b = self.b.copy()
        c = dt.match(_a, _b)
        assert type(c) == Table
        assert len(_a) == len(self.a)
        assert len(_b) == len(self.b)
        assert len(c) == 4

    def test_bkg2d(self):
        b = DetectionRoutine().bkg2d(self.im.copy())
        assert type(b) == type(self.im)
        assert b.shape == self.im.shape


class TestBackground:
    def test_background_estimate_routine_none(self):
        bg = BackGroundEstimateRoutine(None)
        assert bg(None) is None


def test_source_properties_none():
    sp = SourceProperties(None, None)
    assert sp.calculate_crowding() is None
    assert sp.calculate_geometry(1) is None
