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
from photutils.detection import StarFinder
from photutils.psf import IterativePSFPhotometry, ImagePSF

from starbug2.constants import (
    X_CENTROID, Y_CENTROID, OUT_FLUX, FLUX, X_0, Y_0, X_DET, Y_DET, FLUX_FIT,
    ID, STATUS)
from starbug2.utils import Loading, warn, export_table


class ArtificialStarRoutine:
    def __init__(
        self,
        detector: StarFinder,
        psf_fitter: IterativePSFPhotometry,
        psf: ImagePSF):
        """
        :param detector: Detection class that fits the StarFinder base class
        :type detector: photutils.detection.StarFinder
        :param psf_fitter: PSF fitting class that fits the
                           IterativePSFPhotometry base class
        :type psf_fitter: photutils.psf.IterativePSFPhotometry
        :param psf: Empirical Image Point Spread Function model
        :type psf: photutils.psf.ImagePSF
        """
        self._detector: StarFinder = detector
        self._psf_fitter: IterativePSFPhotometry = psf_fitter
        self._psf: ImagePSF = psf

        print("WARNING: THIS IS UNDER DEVELOPMENT")

    def run(self,
            image: np.ndarray,
            n_tests: int = 1000,
            sub_image_size: int = 500,
            sources: Table | None = None,
            full_width_half_max: float = 1.0,
            flux_range: tuple[float, float] | list[float] = (0.0, 1e5),
            separation_thresh: float = 2.0,
            save_progress: bool = True) -> Table:
        # noinspection SpellCheckingInspection
        """
        Run artificial star testing on an image.

        :param image: the image to run artificial star routine on.
        :type image: numpy.ndarray
        :param n_tests: Number of tests to conduct
        :type n_tests: int
        :param sub_image_size: Size of the cropped sub_image
        :type sub_image_size: int
        :param sources: Precalculated positions to test stars with x_0, y_0,
                        and flux columns
        :type sources: astropy.table.Table or None
        :param full_width_half_max: FWHM of the stars to be added (used for
                                    border safety checks).
        :type full_width_half_max: float.
        :param flux_range: Range of fluxes to test.
        :type flux_range: tuple or list.
        :param separation_thresh: Number of pixels above which a detection is
                                  considered a failure.
        :type separation_thresh: float.
        :param save_progress: Periodically save the catalogue during the run.
        :type save_progress: bool.
        :return: Table containing the injected parameters alongside recoveries.
        :rtype: astropy.table.Table.
        """
        shape: np.ndarray = np.array(image.shape)
        if np.any(sub_image_size > shape):
            warn("sub image_size bigger than image dimensions\n")
            sub_image_size = int(min(shape))
        sub_image_size = int(sub_image_size)

        if sources is None:
            x_range: list[float] = [
                2.0 * full_width_half_max,
                float(shape[0] - (2.0 * full_width_half_max))
            ]
            y_range: list[float] = [
                2.0 * full_width_half_max,
                float(shape[1] - (2.0 * full_width_half_max))
            ]
            sources = make_random_models_table(
                int(n_tests),
                {X_0: x_range, Y_0: y_range, FLUX: flux_range},
                seed=int(time.time()))

        # noinspection SpellCheckingInspection
        sources.add_column(Column(np.zeros(len(sources)), name=OUT_FLUX))
        sources.add_column(Column(np.zeros(len(sources)), name=X_DET))
        sources.add_column(Column(np.zeros(len(sources)), name=Y_DET))
        sources.add_column(Column(np.zeros(len(sources)), name=STATUS))

        load: Loading = Loading(len(sources), msg="artificial star tests")
        load.show()

        active_sources: Table = sources

        # noinspection PyTypeChecker
        for n, src in enumerate(iter(active_sources)):
            subx: int = 0
            suby: int = 0

            if sub_image_size > 0:
                subx = int(np.random.randint(
                    max(0,
                        src[X_0] + (2 * full_width_half_max) - sub_image_size),
                    np.min(shape[0] - sub_image_size,
                           src[X_0] - (2 * full_width_half_max))))
                suby = int(np.random.randint(
                    max(0,
                        src[Y_0] + (2 * full_width_half_max) - sub_image_size),
                    np.min(shape[1] - sub_image_size,
                           src[Y_0] - (2 * full_width_half_max))))

            # src mod translates the position within the sub-image
            src_mod: Table = Table(src)
            src_mod[X_0] -= subx
            src_mod[Y_0] -= suby

            sky: np.ndarray = image[
                subx : subx + sub_image_size, suby : suby + sub_image_size]

            # Dynamically extract matrix geometry from sky to support
            # rectangular crops safely
            base: np.ndarray = np.copy(sky) + make_model_image(
                shape=sky.shape, model=self._psf, params_table=src_mod
            )

            detections: Table = self._detector(base)
            detections.rename_column(X_CENTROID, X_0)
            detections.rename_column(Y_CENTROID, Y_0)

            # Check positional matches
            separations: np.ndarray = (
                (src_mod[X_0] - detections[X_0]) ** 2
                + (src_mod[Y_0] - detections[Y_0]) ** 2
            )
            best_match: int = int(np.argmin(separations))

            if np.sqrt(separations[best_match]) <= separation_thresh:
                psf_tab: Table = self._psf_fitter(
                    base, init_params=detections)
                index: np.ndarray = np.where(
                    psf_tab[ID] == detections[best_match][ID])[0]

                if len(index) > 0:
                    matched_idx: int = int(index[0])
                    sources[n][OUT_FLUX] = psf_tab[matched_idx][FLUX_FIT]
                    sources[n][X_DET] = psf_tab[matched_idx][X_0] + subx
                    sources[n][Y_DET] = psf_tab[matched_idx][Y_0] + suby

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