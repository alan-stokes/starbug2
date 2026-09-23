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
from typing import Tuple

import numpy
import os

import numpy as np
from astropy.io.fits import Header, ImageHDU
from astropy.nddata import NDData
from photutils.centroids import centroid_com
from photutils.psf import (
    EPSFBuilder, extract_stars, EPSFStars, EPSFBuildResult, ImagePSF)
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.table import Table
from scipy.ndimage import gaussian_filter

from routines.detection_routines import DetectionRoutine
from starbug2.constants import TableColumn, FileExtensions, ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.utilities.utils import (
    export_table, split_file_name, p_error, printf)


def _epsf_centering_wrapper(data, mask=None):
    """
    Wrapper matching photutils EPSFBuilder signature:
    callable(data, mask=None) -> (x, y)
    """
    clean_data: np.ndarray = numpy.nan_to_num(data, nan=0.0)
    if mask is not None:
        clean_data = clean_data.copy()
        clean_data[mask] = 0.0

    # centroid_com returns (x_center, y_center) cleanly
    return centroid_com(clean_data)


class CustomPSF:
    """
    utilises photutils to build a e-PSF using the image
    """

    @staticmethod
    def generate_epsf(
            sources: Table, data: numpy.ndarray,
            config: StarBugMainConfig) -> EPSFBuildResult:
        """
        generates the epsf.
        :param sources: the stars to use
        :type sources: Table
        :param data: the image data
        :type data: numpy.ndarray
        :param config: the main config
        :type config: StarBugMainConfig
        :return: the generated epsf.
        """
        h_size: float = float((config.custom_psf_size_pixels - 1) / 2)

        # collect stars from the sources
        assert sources is not None
        stars: EPSFStars = CustomPSF.extract_stars(
            sources, h_size, data, config)

        # Filter out invalid cutouts where stars have no flux or
        # has corrupted data.
        valid_stars = [
            star for star in stars
            if numpy.all(numpy.isfinite(star.data)) and star.flux > 0
        ]
        if len(valid_stars) != stars.n_stars:
            printf(f"Have lost {stars.n_stars - len(valid_stars)} stars due"
                   f" to loss of flux or the data is not finite.")

        # turn into the photutils star object.
        clean_stars = EPSFStars(valid_stars)

        # determines which pixels will be ignored. any pixels outside the
        # sigma will be ignored.
        sig_clip = SigmaClip(sigma=config.epsf_clipping_sigma,
                             maxiters=config.epsf_clipping_iterations)

        # Build e-PSF using photutils native sub-pixel alignment
        epsf_builder = EPSFBuilder(
            oversampling=config.epsf_oversampling,
            maxiters=config.epsf_iterations,
            progress_bar=config.verbose_logs,
            sigma_clip=sig_clip,
            smoothing_kernel='quartic',
            recentering_maxiters=config.epsf_centering_iterations
        )
        result: EPSFBuildResult = epsf_builder(clean_stars)
        psf_data = result.epsf.data

        # Post-process: smooth over empty sub-pixel grid cells if oversampled
        if config.epsf_oversampling > 1 and config.epsf_execute_post_smoothing:
            # Sigma scales with oversampling factor to bridge grid gaps
            smooth_sigma = 0.5 * config.epsf_oversampling
            psf_data = gaussian_filter(psf_data, sigma=smooth_sigma)

            # Re-normalise flux to 1.0 after smoothing.
            total_flux = numpy.sum(psf_data)
            if total_flux > 0:
                psf_data /= total_flux
            result.epsf.data = numpy.ascontiguousarray(
                psf_data, dtype=numpy.float64)
        return result

    @staticmethod
    def extract_stars(
            sources: Table, h_size: float, data: numpy.ndarray,
            config: StarBugMainConfig) -> EPSFStars:
        """
        extracts the stars from the image
        :param sources: the locations of the stars
        :type sources: Table
        :param h_size: the size of the extraction in pixels.
        :type h_size: int
        :param data: the image data.
        :type data: numpy.ndarray
        :param config: the main config.
        :type config: StarBugMainConfig
        :return: the stars in the format photutils requires.
        :rtype: EPSFStars
        """
        clean_data = numpy.ascontiguousarray(
            numpy.nan_to_num(data.copy(), nan=0.0, posinf=0.0, neginf=0.0),
            dtype=numpy.float64
        )
        nd_data = NDData(data=clean_data)

        # Convert FITS 1-based coordinates to Python 0-based array indices
        x_raw = numpy.array(
            sources[TableColumn.X_CENTROID], dtype=numpy.float64) - 1.0
        y_raw = numpy.array(
            sources[TableColumn.Y_CENTROID], dtype=numpy.float64) - 1.0

        # Mask boundary stars based on 0-indexed bounds
        mask = (
            numpy.isfinite(x_raw) & numpy.isfinite(y_raw) &
            (x_raw > h_size) & (x_raw < (clean_data.shape[1] - 1 - h_size)) &
            (y_raw > h_size) & (y_raw < (clean_data.shape[0] - 1 - h_size))
        )

        stars_tbl = Table()
        stars_tbl['x'] = x_raw[mask]
        stars_tbl['y'] = y_raw[mask]

        return extract_stars(
            nd_data, stars_tbl, size=config.custom_psf_size_pixels)

    @staticmethod
    def execute_background_extraction(base: StarbugBase) -> Tuple[
                ExitStates, numpy.ndarray | None]:
        result_state: ExitStates = base.bgd_estimate()
        if result_state != ExitStates.EXIT_SUCCESS:
            p_error("Failed to execute bgd_estimate")
            return result_state, None

        result_state = base.bgd_subtraction()
        if result_state != ExitStates.EXIT_SUCCESS:
            p_error("Failed to execute bgd_subtraction")
            return result_state, None

        return ExitStates.EXIT_SUCCESS, base.residues

    @staticmethod
    def execute_custom_e_psf(config: StarBugMainConfig) -> ExitStates:
        """
        generates an epsf from photutils.
        follows the sample code from:
            https://photutils.readthedocs.io/en/latest/user_guide/
            epsf_building.html
        :param config: the config object
        :type config: StarBugMainConfig
        :return: success if complete.
        :rtype: ExitStates
        """

        # read in image from fits file
        base: StarbugBase = StarbugBase(
            config.fits_images[0], config, ap_file=None, bkg_file=None)
        data: numpy.ndarray = base.main_image().data

        # locate stars.
        sources: Table | None = CustomPSF.get_psf_sources(
            data, config, base.full_width_half_max)

        # remove background from the data
        (exit_states, data_bkg_removed) = (
            CustomPSF.execute_background_extraction(base))
        if exit_states != ExitStates.EXIT_SUCCESS:
            return exit_states

        # generate psf.
        assert sources is not None
        assert data_bkg_removed is not None
        result: EPSFBuildResult = CustomPSF.generate_epsf(
            sources, data, config)

        # extract e-PSF from the builder/
        epsf: ImagePSF = result.epsf
        fitted_stars: EPSFStars = result.fitted_stars

        output_dir: str | None = config.output_file
        _, b_name, _ = split_file_name(base.f_name)
        assert output_dir is not None
        return CustomPSF.write_files_to_disk(
            output_dir, epsf, fitted_stars, b_name)

    @staticmethod
    def get_psf_sources(
            data: numpy.ndarray, config: StarBugMainConfig,
            full_width_half_max: float) -> Table | None:
        """
        extracts sources based off DAOStarFinder.
        :param data: the image data
        :type data: numpy.ndarray
        :param config: the main config
        :type config: StarBugMainConfig
        :param full_width_half_max: the full width 1/2 max value.
        :type full_width_half_max: float
        :return: the sources as a table format. contains columns of:
            ['id', 'x_centroid', 'y_centroid', 'sharpness', 'roundness1',
            'roundness2', 'n_pixels', 'peak', 'flux', 'mag']
        :rtype: Table
        """
        # determine threshold.
        median_stat: float
        std: float
        _, median_stat, std = sigma_clipped_stats(data, sigma=config.sigma_sky)

        # locate stars in fits image
        detector: DetectionRoutine = DetectionRoutine(
            sig_src=config.sigma_source,
            sig_sky=config.sigma_sky,
            full_width_half_max=full_width_half_max,
            sharp_lo=config.sharp_cutoff_low,
            sharp_hi=config.sharp_cutoff_high,
            round_1_hi=config.round1_cutoff_high,
            round_2_hi=config.round2_cutoff_high,
            smooth_lo=config.smooth_low,
            smooth_hi=config.smooth_high,
            ricker_r=config.ricker_wavelet_radius,
            do_bgd_2d=config.do_bgd_2d,
            do_con_vl=config.do_convolution,
            box_size=config.background_box_size,
            clean_src=config.clean_sources,
            verbose=config.verbose_logs)

        return detector(data.copy())

    @staticmethod
    def write_files_to_disk(
            output_dir: str, epsf: ImagePSF,
            fitted_stars: EPSFStars, fits_file_name: str) -> ExitStates:
        """
        writes the new psf and the detected stars into files.

        :param output_dir: the output dir
        :type output_dir: str
        :param epsf: the psf object
        :type epsf: ImagePSF
        :param fitted_stars: the stars being fitted.
        :type fitted_stars: EPSFStars
        :param fits_file_name: the fits file name.
        :type fits_file_name: str
        :return: success if done
        :rtype: ExitStates
        """
        # write new psf into a .fits file for further use.
        new_psf_header: Header = Header()
        assert output_dir is not None
        file_name: str = os.path.join(
            output_dir, f"{fits_file_name}_custom{FileExtensions.CUSTOM_PSF}")
        ImageHDU(data=epsf.data, header=new_psf_header).writeto(
            file_name, overwrite=True)

        # write detected stars as a .ap file
        star_data = []
        for star in fitted_stars:
            star_data.append(
                {
                    TableColumn.ID:
                        star.id if hasattr(star, TableColumn.ID) else None,
                    TableColumn.X_FIT: star.center[0],
                    TableColumn.Y_FIT: star.center[1],
                    TableColumn.FLUX_FIT:
                        getattr(star, TableColumn.FLUX, None),
                    "cutout_x": star.origin[0],
                    "cutout_y": star.origin[1],
                }
            )

        # Create an Astropy Table
        stars_table = Table(star_data)
        stars_table.remove_column(TableColumn.ID)
        custom_file_name: str = os.path.join(
            output_dir,
            f"{fits_file_name}_custom_fit_stars{FileExtensions.AP}")
        export_table(stars_table, custom_file_name, header=new_psf_header)

        return ExitStates.EXIT_SUCCESS
