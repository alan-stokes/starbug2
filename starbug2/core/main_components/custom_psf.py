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
import numpy
import os
from astropy.io.fits import Header, ImageHDU
from astropy.nddata import NDData
from photutils.psf import (
    EPSFBuilder, extract_stars, EPSFStars, EPSFBuildResult, ImagePSF)
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, Column

from routines.detection_routines import DetectionRoutine
from starbug2.constants import TableColumn, FileExtensions, ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.utilities.utils import export_table, split_file_name
from utilities.filters import STAR_BUG_FILTERS


class CustomPSF:
    """
    utilises photutils to build a e-PSF using the image
    """

    @staticmethod
    def generate_epsf(
            sources: Table, data: numpy.ndarray,
            config: StarBugMainConfig) -> EPSFBuildResult:
        """
        takes sources and image data as well as hsize and config and generates
        an epsf result.

        :param sources: the sources Table.
        :type sources: Table
        :param data: the image data.
        :type data: numpy.ndarray
        :param config: the system config.
        :type config: StarBugMainConfig
        :return: the epsf result object.
        :rtype: EPSFBuildResult
        """
        h_size: float = float((config.custom_psf_size_pixels - 1) / 2)

        # locate stars not on the edge of the image (within capture area).
        assert sources is not None
        stars: EPSFStars = CustomPSF.extract_stars(
            sources, h_size, data, config)

        # build e-PSF.
        epsf_builder: EPSFBuilder = EPSFBuilder(
            oversampling=4, maxiters=3, progress_bar=config.verbose_logs)
        result: EPSFBuildResult = epsf_builder(stars)
        return result

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

        # generate psf.
        assert sources is not None
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
            'roundness2', 'n_pixels', 'peak', 'flux', 'mag', 'daofind_mag']
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
    def extract_stars(
            sources: Table, h_size: float, data: numpy.ndarray,
            config: StarBugMainConfig) -> EPSFStars:
        """ extracts stars from the image.

        :param sources: the star locations.
        :type sources: Table
        :param h_size: the size of the window for detecting stars.
        :type h_size: float
        :param data: the image data.
        :type data: numpy.ndarray
        :param config: the config
        :type config: StarBugMainConfig
        :return: the extracted stars
        :rtype: EPSFStars
        """

        # determine states
        mean_val: float
        median_val: float
        std_val: float
        mean_val, median_val, std_val = sigma_clipped_stats(
            data, sigma=config.sigma_sky)
        data -= median_val

        # extract stars from the modified data.
        nd_data: NDData = NDData(data=data)

        assert sources is not None
        x: Column = sources[TableColumn.X_CENTROID]
        y: Column = sources[TableColumn.Y_CENTROID]
        mask: numpy.ndarray = (
            (x > h_size) & (x < (data.shape[1] - 1 - h_size)) &
            (y > h_size) & (y < (data.shape[0] - 1 - h_size)))

        stars_tbl = Table()
        stars_tbl[TableColumn.X] = x[mask]
        stars_tbl[TableColumn.Y] = y[mask]

        return extract_stars(
            nd_data, stars_tbl, size=config.custom_psf_size_pixels)

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
