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

"""
Core routines for StarbugII.
"""
import time
import numpy as np
from astropy.table import Column, Table
from photutils.datasets import make_model_image, make_random_models_table

from starbug2.constants import (
    X_CENTROID, Y_CENTROID, OUT_FLUX, FLUX, X_0, Y_0, X_DET, Y_DET, FLUX_FIT,
    ID)
from starbug2.utils import Loading, warn, export_table


class ArtificialStarRoutine(object):
    def __init__(self, detector, psf_fitter, psf):
        """
        :param detector: Detection class that fits the StarFinder base class
        :type detector: photutils.detection.StarFinder
        :param psf_fitter: PSF fitting class that fits the
                           IterativelySubtractedPSFPhotometry base class
        :type psf_fitter: photutils.psf.IterativelySubtractedPSFPhotometry
        :param psf: Discrete Point Response Function
        :type psf: photutils.psf.DiscretePRF
        """
        self._detector = detector
        self._psf_fitter = psf_fitter
        self._psf = psf

        print("WARNING: THIS IS UNDER DEVELOPMENT")


    def run(self, image, n_tests=1000, sub_image_size=500, sources=None,
            full_width_half_max=1, flux_range=(0,1e5), separation_thresh=2,
            save_progress=1):
        # noinspection SpellCheckingInspection
        """
        Run artificial star testing on an image

        :param image: the image to run artifical star routine one.
        :type image: numpy.ndarray
        :param n_tests: Number of tests to conduct
        :type n_tests: int
        :param sub_image_size: Size of the cropped sub_image
        :type sub_image_size: int
        :param sources: Precalculated positions to test stars with x_0, y_0,
                        and flux columns
        :type sources: astropy.table.Table
        :param full_width_half_max: FWHM of the stars to be added (used for
                                    border safety checks)
        :type full_width_half_max: float
        :param flux_range: Range of fluxes to test
        :type flux_range: list or tuple
        :param separation_thresh: Number of pixels above which a detection is
                                  considered a failure
        :type separation_thresh: float
        :param save_progress: Periodically save the catalogue during the run
        :type save_progress: bool
        """
        shape = np.array(image.shape)
        if np.any(sub_image_size > shape):
            warn("sub image_size bigger than image dimensions\n")
            sub_image_size = min(shape)
        sub_image_size = int(sub_image_size)

        if not sources:
            x_range = [
                2.0 * full_width_half_max, shape[0]
                - (2.0 * full_width_half_max)]
            y_range = [
                2.0 * full_width_half_max, shape[1]
                - (2.0 * full_width_half_max)]
            sources = make_random_models_table(
                int(n_tests),
                {X_0: x_range, Y_0: y_range, FLUX: flux_range},
                seed=int(time.time()))

        # noinspection SpellCheckingInspection
        sources.add_column(Column(np.zeros(len(sources)), name="outflux"))
        sources.add_column(Column(np.zeros(len(sources)), name="x_det"))
        sources.add_column(Column(np.zeros(len(sources)), name="y_det"))
        sources.add_column(Column(np.zeros(len(sources)), name="status"))

        load = Loading(len(sources), msg="artificial star tests")
        load.show()
        for n, src in enumerate(sources):

            subx = 0
            suby = 0
            if sub_image_size > 0:
                ## !! I might change this to be PSF_SIZE not
                # 2_Full_width_1/2_max
                subx = np.random.randint(
                    max(0,
                        src[X_0] + (2 * full_width_half_max)
                        - sub_image_size),
                    np.min(shape[0] - sub_image_size,
                        src[X_0] - (2 * full_width_half_max)))
                suby = np.random.randint(
                    max(0,
                        src[Y_0] + (2 * full_width_half_max)
                        - sub_image_size),
                    np.min(shape[1] - sub_image_size,
                        src[Y_0] - (2 * full_width_half_max)))

            # src mod translates the position within the sub-image
            src_mod = Table(src)
            src_mod[X_0] -= subx
            src_mod[Y_0] -= suby
            sky = image[
                  subx : subx + sub_image_size,
                  suby : suby + sub_image_size]
            base = np.copy(sky) + make_model_image(
                2 * [sub_image_size], self._psf, src_mod)

            detections = self._detector(base)
            detections.rename_column(X_CENTROID, X_0)
            detections.rename_column(Y_CENTROID, Y_0)

            separations = (
                (src_mod[X_0] - detections[X_0]) ** 2
                + (src_mod[Y_0] - detections[Y_0]) ** 2)
            best_match = np.argmin(separations)
            if np.sqrt(separations[best_match]) <= separation_thresh:
                psf_tab = self._psf_fitter(base, init_guesses=detections)
                index = np.where(psf_tab[ID] == detections[best_match][ID])

                sources[n][OUT_FLUX] = psf_tab[index][FLUX_FIT]
                sources[n][X_DET] = psf_tab[index][X_0] + subx
                sources[n][Y_DET] = psf_tab[index][Y_0] + suby

                if (abs(sources[n][OUT_FLUX] - sources[n][FLUX])
                        < (sources[n][FLUX] / 100.0)):
                    # star matched
                    sources[n]["status"] = 1
            load()
            load.show()

            if save_progress and not n % 10:
                export_table(
                    sources[0:n], f_name="/tmp/artificial_stars.save")

        return sources
