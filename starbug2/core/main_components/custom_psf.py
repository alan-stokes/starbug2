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
from astropy.io.fits import ImageHDU, PrimaryHDU, Header, BinTableHDU
from astropy.nddata import NDData
from photutils.detection import DAOStarFinder
from photutils.psf import (
    EPSFBuilder, extract_stars, EPSFStars, EPSFBuildResult, ImagePSF)
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, Column

from starbug2.constants import TableColumn, FileExtensions, ExitStates
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.utilities.utils import export_table


class CustomPSF:
    """
    utilises photutils to build a e-PSF using the image
    """

    @staticmethod
    def execute_custom_e_psf(config: StarBugMainConfig) -> ExitStates:
        """
        generates a epsf from photutils.
        follows the sample code from:
            https://photutils.readthedocs.io/en/latest/user_guide/
            epsf_building.html
        :param config:
        :return:
        """

        # read in image from fits file
        base: StarbugBase = StarbugBase(
            config.fits_images[0], config, ap_file=None, bkg_file=None)
        data: ImageHDU | PrimaryHDU = base.main_image()

        # determine threshold.
        median_stat: float
        std: float
        _, median_stat, std = sigma_clipped_stats(data, sigma=config.sigma_sky)
        threshold = std * config.sigma_source

        # locate stars in fits image
        finder: DAOStarFinder = DAOStarFinder(
            threshold=threshold, fwhm=config.full_width_half_max)
        sources: Table | None = finder(data)
        hsize: float = float((config.custom_psf_size_pixels - 1) / 2)

        # locate stars not on the edge of the image (within capture area).
        stars: EPSFStars = CustomPSF.extract_stars()

        # build e-PSF.
        epsf_builder: EPSFBuilder = EPSFBuilder(
            oversampling=4, maxiters=3, progress_bar=config.verbose_logs)
        result: EPSFBuildResult = epsf_builder(stars)

        # extract e-PSF from the builder/
        epsf: ImagePSF = result.epsf
        fitted_stars: EPSFStars = result.fitted_stars

        output_dir: str | None = config.output_file
        assert output_dir is not None
        return CustomPSF.write_files_to_disk(output_dir, epsf, fitted_stars)

    @staticmethod
    def extract_stars(sources, h_size, data, config) -> EPSFStars:
        assert sources is not None
        x: Column = sources[TableColumn.X_CENTROID]
        y: Column = sources[TableColumn.Y_CENTROID]
        mask: numpy.ndarray = (
            (x > h_size) & (x < (data.shape[1] - 1 - h_size)) &
            (y > h_size) & (y < (data.shape[0] - 1 - h_size)))

        stars_tbl = Table()
        stars_tbl[TableColumn.X] = x[mask]
        stars_tbl[TableColumn.Y] = y[mask]

        # determine states
        mean_val: float
        median_val: float
        std_val: float
        mean_val, median_val, std_val = sigma_clipped_stats(
            data, sigma=config.sigma_sky)
        data -= median_val

        # extract stars from the modified data.
        nd_data: NDData = NDData(data=data)
        return extract_stars(
            nd_data, stars_tbl, size=config.custom_psf_size_pixels)

    @staticmethod
    def write_files_to_disk(
            output_dir: str, epsf: ImagePSF,
        fitted_stars: EPSFStars) -> ExitStates:
        # write new psf into a .fits file for further use.
        new_psf_header: Header = Header()
        assert output_dir is not None
        file_name: str = os.path.join(
            output_dir, f"custom{FileExtensions.CUSTOM_PSF}")
        BinTableHDU(data=epsf, header=new_psf_header).writeto(
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
        custom_file_name: str = os.path.join(
            output_dir, f"custom_fit_stars{FileExtensions.AP}")
        export_table(stars_table, custom_file_name, header=new_psf_header)

        return ExitStates.EXIT_SUCCESS