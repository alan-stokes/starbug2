import numpy
from astropy.io.fits import ImageHDU, PrimaryHDU
from astropy.nddata import NDData
from photutils.detection import DAOStarFinder
from photutils.psf import (
    EPSFBuilder, extract_stars, EPSFStars, EPSFBuildResult, ImagePSF)
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, Column

from constants import TableColumn
from core.star_bug_config import StarBugMainConfig
from core.starbug_main import StarbugBase


class CustomPSF:
    """
    utilises photutils to build a e-PSF using the image
    """

    @staticmethod
    def execute_custom_e_psf(config: StarBugMainConfig):
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
        assert sources is not None
        x: Column = sources[TableColumn.X_CENTROID]
        y: Column = sources[TableColumn.Y_CENTROID]
        mask: numpy.ndarray = (
            (x > hsize) & (x < (data.shape[1] - 1 - hsize)) &
            (y > hsize) & (y < (data.shape[0] - 1 - hsize)))

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
        stars: EPSFStars = extract_stars(
            nd_data, stars_tbl, size=config.custom_psf_size_pixels)

        # build e-PSF.
        epsf_builder: EPSFBuilder = EPSFBuilder(
            oversampling=4, maxiters=3, progress_bar=config.verbose_logs)
        result: EPSFBuildResult = epsf_builder(stars)

        # extract e-PSF from the builder/
        epsf: ImagePSF = result.epsf
        fitted_stars: EPSFStars = result.fitted_stars




