import os
from typing import List, Callable, Tuple, Any, cast

import numpy as np
from astropy.io.fits import (
    HDUList, open, Header, ImageHDU, BinTableHDU, PrimaryHDU)
from astropy.stats import sigma_clipped_stats
from astropy.table import Table, Column, vstack
from astropy.wcs import WCS
from photutils.datasets import make_model_image

from photutils.psf import ImagePSF
from starbug2.core.constants import (
    TableColumn, ExitStates, HeaderTags, SourceFlags, Units, ImageHeaderTags,
    MIRI_STRING, MIRI_IMAGE, DetectorLengths, NIRCAM, BGD_FILE, AP_FILE)
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.interfaces.star_bug_interface import StarBugInterface
from starbug2.routines.psf_phot_routine import PSFPhotRoutine
from starbug2.utilities.filters import STAR_BUG_FILTERS, FilterStruct
from starbug2.utilities.utils import (
    flux2mag, p_error, parse_unit, warn, get_data_path,
    get_mj_ysr2jy_scale_factor, reindex)


class Photometry:

    @staticmethod
    def determine_f_name(filter_string: str | None, f_name: str | None,
                         info: dict[str, str]) -> str | None:
        """
        Finds the file name for the PSF if not provided
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :param f_name: Filename of a PSF fits image
        :type f_name: str
        :param info: Metadata dictionary detailing execution properties.
        :type info: dict[str, str]
        :return: the correct psf file name
        :rtype: str
        """
        assert filter_string is not None
        if not f_name:
            filter_struct: FilterStruct | None = (
                STAR_BUG_FILTERS.get(filter_string))
            if filter_struct:
                dt_name: str = info[ImageHeaderTags.DETECTOR]
                if dt_name == "NRCALONG":
                    dt_name = "NRCA5"
                if dt_name == "NRCBLONG":
                    dt_name = "NRCB5"
                if dt_name == "MULTIPLE":
                    if (filter_struct.instr == NIRCAM
                            and filter_struct.length == DetectorLengths.SHORT):
                        dt_name = "NRCA1"
                    elif (filter_struct.instr == NIRCAM and
                          filter_struct.length == DetectorLengths.LONG):
                        dt_name = "NRCA5"
                    elif filter_struct.instr == MIRI_STRING:
                        dt_name = ""
                if dt_name == MIRI_IMAGE:
                    dt_name = ""
                f_name = "%s/%s%s.fits" % (
                    get_data_path(), filter_string, dt_name)
            else:
                f_name = None
        return f_name

    # noinspection SpellCheckingInspection
    @staticmethod
    def load_psf(
            filter_string: str | None, info: dict[str, str],
            log: Callable[[str], None],
            f_name: str | None = None) -> Tuple[ExitStates, np.ndarray | None]:
        """
        Load a PSF_FILE to be used during photometry.
        :param filter_string: Name of the photometric filter band used.
        :type filter_string: str | None
        :param info: Metadata dictionary detailing execution properties.
        :type info: dict[str, str]
        :param log: Callable logging function for status tracking.
        :type log: Callable[[str], None]
        :param f_name: Filename of a PSF fits image
        :type f_name: str
        :return: The exit status state and the psf
        :rtype: Tuple[ExitStates, np.ndarray | None]
        """
        status: ExitStates = ExitStates.EXIT_SUCCESS
        f_name = Photometry.determine_f_name(filter_string, f_name, info)

        psf: np.ndarray | None = None
        if f_name is not None and os.path.exists(f_name):
            fp: HDUList = open(f_name)

            if fp[0].data is None:
                p_error(
                    "There is a version mismatch between starbug and "
                    "webbpsf. Please reinitialise with: starbug2 --init.\n")
                quit("Fatal error, quitting\n")

            psf = fp[0].data
            fp.close()
            log("loaded PSF_FILE='%s'\n" % f_name)
        else:
            p_error("PSF_FILE='%s' does not exist\n" % f_name)
            status = ExitStates.EXIT_FAIL
        return status, psf

    @staticmethod
    def photometry_routine(
            filter_string: str | None, wcs: WCS | None,
            config: StarBugMainConfig, main_image: ImageHDU | PrimaryHDU,
            log: Callable[[str], None], image: HDUList | None,
            info: dict[str, str], background: ImageHDU | PrimaryHDU | None,
            header: Header, detections: Table | None, psf: np.ndarray | None,
            full_width_half_max: float, ap_file: str | None,
            background_file: str | None, out_dir: str | None,
            b_name: str | None) -> Tuple[int, Table | None, np.ndarray | None]:
        """
        Full photometry routine
        Saves the result as a table self._psf_catalogue,
        Additionally it appends a residual Image onto the
        self._residuals HDUList

        :return: success state, the psf_catalogue, and residuals
        :rtype Tuple[int, Table | None, np.ndarray | None]
        """
        if filter_string is None or wcs is None:
            return ExitStates.EXIT_FAIL, None, None

        if main_image:
            log("\nRunning PSF Photometry\n")

            # lock the types.
            image_data: np.ndarray
            error: np.ndarray
            bgd: np.ndarray | None
            mask: np.ndarray
            image_data, error, bgd, mask = (
                StarBugInterface.prepare_image_arrays(
                    image, log, background, header, main_image))

            if bgd is None:
                clipped_median: float
                _, clipped_median, _ = (
                    sigma_clipped_stats(image_data, sigma=config.sigma_sky))
                bgd = np.ones(main_image.shape) * clipped_median
                log(
                    "-> no background file loaded, measuring sigma "
                    "clipped median\n")
            assert bgd is not None

            # Collect relevant files and data
            if detections is None:
                p_error("unable to run photometry: no source list loaded\n")
                return ExitStates.EXIT_FAIL, None, None

            if psf is None:
                result: ExitStates
                (result, psf) = Photometry.load_psf(
                    filter_string, info, log,
                    os.path.expandvars(config.psf_file_override))
                if result != ExitStates.EXIT_SUCCESS:
                    p_error("unable to run photometry: no PSF loaded\n")
                    return ExitStates.EXIT_FAIL, None, None

            psf_mask: np.ndarray = ~np.isfinite(psf)
            if psf_mask.sum():
                assert psf is not None
                psf[psf_mask] = 0
                log("-> masking INF pixels in PSF_FILE\n")

            psf_model: ImagePSF = ImagePSF(data=psf)
            size: int
            if config.psf_fit_size > 0:
                size = config.psf_fit_size
            else:
                size = psf_model.shape[0]
            if not size % 2:
                size -= 1
            log("-> psf size: %d\n" % size)

            # Sort out Init guesses
            app_hot_r: float = float(config.aperture_phot_radius)
            if not app_hot_r or app_hot_r <= 0:
                app_hot_r = 3.0

            init_guesses: Table = detections.copy()

            # noinspection DuplicatedCode
            if TableColumn.X_CENTROID in init_guesses.colnames:
                init_guesses.rename_column(
                    TableColumn.X_CENTROID, TableColumn.X_INIT)
            if TableColumn.Y_CENTROID in init_guesses.colnames:
                init_guesses.rename_column(
                    TableColumn.Y_CENTROID, TableColumn.Y_INIT)
            if TableColumn.X_DET in init_guesses.colnames:
                init_guesses.rename_column(
                    TableColumn.X_DET, TableColumn.X_INIT)
            if TableColumn.Y_DET in init_guesses.colnames:
                init_guesses.rename_column(
                    TableColumn.Y_DET, TableColumn.Y_INIT)

            init_guesses = init_guesses[init_guesses[TableColumn.X_INIT] >= 0]
            init_guesses = init_guesses[init_guesses[TableColumn.Y_INIT] >= 0]
            init_guesses = init_guesses[
                init_guesses[TableColumn.X_INIT]
                < main_image.header[HeaderTags.NAXIS1]]
            init_guesses = init_guesses[
                init_guesses[TableColumn.Y_INIT]
                < main_image.header[HeaderTags.NAXIS2]]

            # Allow tables that don't have the correct columns through
            required: List[str] = [
                TableColumn.X_INIT, TableColumn.Y_INIT, TableColumn.FLUX,
                filter_string, TableColumn.FLAG]
            for notfound in set(required) - set(init_guesses.colnames):
                dtype = np.uint16 if notfound == TableColumn.FLAG else float
                init_guesses.add_column(
                    Column(np.zeros(len(init_guesses)),
                           name=notfound, dtype=dtype))

            init_guesses = init_guesses[required]
            init_guesses.remove_column(TableColumn.FLUX)
            init_guesses.rename_column(filter_string, "ap_%s" % filter_string)

            # Run Fit
            min_separation: float = config.critical_separation
            if not min_separation:
                min_separation = min(5.0, 2.5 * full_width_half_max)

            if config.force_centroid_position:
                phot: PSFPhotRoutine = PSFPhotRoutine(
                    psf_model, size, min_separation=min_separation,
                    app_hot_r=app_hot_r, background=bgd, force_fit=1,
                    verbose=config.verbose_logs)
                psf_cat: Table = phot(
                    image_data, init_params=init_guesses, error=error,
                    mask=mask)
                psf_cat[TableColumn.FLAG] |= SourceFlags.SRC_FIX

            else:
                phot = PSFPhotRoutine(
                    psf_model, size, min_separation=min_separation,
                    app_hot_r=app_hot_r, background=bgd, force_fit=0,
                    verbose=config.verbose_logs)
                psf_cat = phot(
                    image_data, init_params=init_guesses, error=error,
                    mask=mask)

                if not psf_cat:
                    return ExitStates.EXIT_FAIL, None, None

                # Setting position max variation
                max_y_dev: float
                unit: int
                max_y_dev, unit = parse_unit(config.max_xy_deviation)
                if unit is not None:
                    if unit == Units.DEG:
                        max_y_dev *= 60
                        unit = Units.ARCMIN
                    if unit == Units.ARCMIN:
                        max_y_dev *= 60
                        unit = Units.ARCSEC
                    if unit == Units.ARCSEC:
                        if not header.get(ImageHeaderTags.PIXAR_A2):
                            warn(
                                "MAX_XYDEV is units arcseconds, but starbug "
                                "cannot locate a pixel scale in the header."
                                " Please use syntax MAX_XYDEV=%sp to set "
                                "change to pixels\n" % max_y_dev)
                        else:
                            max_y_dev /= np.sqrt(
                                header.get(ImageHeaderTags.PIXAR_A2))

                if max_y_dev > 0:
                    log(
                        "-> position fit threshold: %.2gpix\n" % max_y_dev)
                    phot = PSFPhotRoutine(
                        psf_model, size, min_separation=min_separation,
                        app_hot_r=app_hot_r, background=bgd, force_fit=1,
                        verbose=config.verbose_logs)
                    ii: np.ndarray = psf_cat[TableColumn.XY_DEV] > max_y_dev
                    fixed_centres: Table = psf_cat[ii][
                        [TableColumn.X_INIT, TableColumn.Y_INIT,
                         "ap_%s" % filter_string, TableColumn.FLAG]]
                    if len(fixed_centres):
                        log("-> forcing positions for deviant sources\n")
                        fixed_cat: Table = phot(
                            image_data, init_params=fixed_centres,
                            error=error, mask=mask)
                        # ABS. why are we using such an aggressive type check
                        # here?
                        fixed_cat[TableColumn.FLAG] |= (
                            np.uint16(SourceFlags.SRC_FIX))
                        psf_cat.remove_rows(ii.tolist())
                        psf_cat = vstack((psf_cat, fixed_cat))
                    else:
                        log("-> no deviant sources\n")
            ra: np.ndarray
            dec: np.ndarray
            ra, dec = wcs.all_pix2world(
                psf_cat[TableColumn.X_FIT], psf_cat[TableColumn.Y_FIT], 0)
            psf_cat.add_column(Column(ra, name=TableColumn.RA), index=2)
            psf_cat.add_column(Column(dec, name=TableColumn.DEC), index=3)

            mag: float
            mag_err: float
            mag, mag_err = flux2mag(
                psf_cat[TableColumn.FLUX], psf_cat[TableColumn.E_FLUX])

            filter_string: str = (
                filter_string if filter_string else TableColumn.MAG)
            psf_cat.add_column(
                mag + config.zero_point_magnitude, name=filter_string)
            psf_cat.add_column(mag_err, name="e%s" % filter_string)
            psf_catalogue = psf_cat

            # verify catalogue isn't none
            assert psf_catalogue is not None

            psf_catalogue.meta = dict(header.items())
            psf_catalogue.meta[AP_FILE] = ap_file
            psf_catalogue.meta[BGD_FILE] = background_file

            reindex(psf_catalogue)

            file_name: str = (
                "%s/%s-psf.fits" % (out_dir, b_name))
            log("--> %s\n" % file_name)
            BinTableHDU(
                data=psf_catalogue,
                header=header).writeto(file_name, overwrite=True)

            # Residual Image
            residuals: np.ndarray | None = None
            if config.generate_residual_image:
                log("-> generating residual\n")
                _tmp: Table = psf_cat[
                    TableColumn.X_FIT, TableColumn.Y_FIT,
                    TableColumn.FLUX].copy()
                _tmp.rename_columns(
                    (TableColumn.X_FIT, TableColumn.Y_FIT),
                    (TableColumn.X_0, TableColumn.Y_0))
                stars: np.ndarray = make_model_image(
                    image_data.shape, psf_model, _tmp,
                    model_shape=(size, size))
                residual: np.ndarray = image_data - (bgd + stars)
                residuals = (
                    residual / get_mj_ysr2jy_scale_factor(main_image))
                header: Header = header
                header.update(wcs.to_header())
                ImageHDU(
                    data=cast(Any, residuals),
                    name="RES", header=header).writeto(
                    "%s/%s-res.fits" % (out_dir, b_name),
                    overwrite=True)
            return ExitStates.EXIT_SUCCESS, psf_catalogue, residuals
        return ExitStates.EXIT_FAIL, None, None
