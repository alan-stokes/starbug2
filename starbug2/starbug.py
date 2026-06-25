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

import os
import sys
from os import getenv
from typing import Final, Optional, Tuple, Dict, List, cast, Any

from astropy.wcs import (
    WCS, NoConvergence, SingularMatrixError, InconsistentAxisTypesError,
    InvalidTransformError)
import numpy as np
from astropy.io.fits import (
    PrimaryHDU, ImageHDU, HDUList, Header, open, BinTableHDU)
from astropy.table import hstack, Column, vstack, Table
from astropy.stats import sigma_clipped_stats
from photutils.datasets import make_model_image
from photutils.psf import ImagePSF
from starbug2.constants import (
    FILTER, STAR_BUG, CALIBRATION_LV, DETECTOR, TELESCOPE, INSTRUMENT, BUN_IT,
    PIXAR_A2, PIXAR_SR, SCI, BGD, RES, VERBOSE_TAG, AP_FILE,
    BGD_FILE, FITS_EXTENSION, JWST, DQ, AREA, WHT, RA,
    DEC, X_CENTROID, Y_CENTROID, SHORT, LONG, NIRCAM, STAR_BUG_MIRI, SRC_FIX,
    DEG, ARCMIN, ARCSEC, DQ_DO_NOT_USE,
    DQ_SATURATED, NAXIS1, NAXIS2, ERR, EXIT_SUCCESS, EXIT_FAIL,
    NIRCAM_STRING,
    STARBUG_DATA_DIR, FLUX, E_FLUX, X_INIT, Y_INIT,
    X_DET, Y_DET, FLAG, XY_DEV, X_FIT, Y_FIT, MAG, X_0, Y_0,
    DEFAULT_FULL_WIDTH_HALF_MAX, FLUX_DET)
from starbug2.filters import STAR_BUG_FILTERS, FilterStruct
from starbug2.routines.app_hot_routine import APPhotRoutine
from starbug2.routines.background_estimate_routine import (
    BackGroundEstimateRoutine)
from starbug2.routines.detection_routines import DetectionRoutine
from starbug2.routines.psf_phot_routine import PSFPhotRoutine
from starbug2.routines.source_properties import SourceProperties
from starbug2.star_bug_config import StarBugMainConfig
from starbug2.star_bug_interface import StarBugInterface
from starbug2.utils import (
    collapse_header, parse_unit, get_version, ext_names, printf,
    split_file_name, p_error, warn, import_table, get_mj_ysr2jy_scale_factor,
    flux2mag, reindex, export_table)


class StarbugBase(StarBugInterface):
    """
    StarbugBase is the overall container for the photometry package. It holds
    the active image, the parameter file and the output images/tables.
    It is self-contained enough to simply run "photometry" and everything
    should just take care of itself from there on.
    """

    MIN_MAG: Final[int] = 27
    MAX_MAG: Final[int] = 18

    @staticmethod
    def get_data_path() -> str:
        """
        returns the data path.
        :return: the data path
        """
        env_path: str | None = getenv(STARBUG_DATA_DIR)
        return (env_path if env_path else
                "%s/.local/share/starbug" % (getenv("HOME")))

    @staticmethod
    def sort_output_names(
            f_name: str, param_output: Optional[str]=None) -> (
                Tuple[str, str, str]):
        """
        This is a useful function that looks at both an input file and a set
        output and figures out how to name output files. If param_output looks
        like a directory then the output will be set to that directory with
        the basename of f_name. If param_output looks like a file, then the
        output basename will take that form

        :param f_name: Filename to use as the core of the output
        :type f_name: str
        :param param_output: This is the OUTPUT parameter in the parameter
        file. It can be an output directory or output filename. If None
        (default) then it will be ignored
        :type param_output: str
        :return: a tuple of (
             The output directory,
             The output file basename (e.g. path/to/output.txt -> "output"),
             The file extension split from the inputs)
        :rtype tuple of (str, str, str)
        """

        out_dir: str = ""
        b_name: str = ""
        extension: str = ""
        if f_name:
            out_dir, b_name, extension = split_file_name(f_name)
            if (tmp_out_name := param_output) and tmp_out_name != '.':
                inner_out_dir, inner_b_name, _= split_file_name(tmp_out_name)
                if os.path.exists(out_dir) and os.path.isdir(out_dir):
                    out_dir = inner_out_dir
                else:
                    p_error("unable to locate output directory \"%s\"\n" %
                            inner_out_dir)
                if inner_b_name:
                    b_name = inner_b_name

        return out_dir, b_name, extension

    def __init__(
            self, f_name: str, config: StarBugMainConfig,
            ap_file, bkg_file, verbose) -> None:
        """
         Star bug init.

        :param f_name: FITS image file name
        :type f_name: str
        :param config: the starbug config
        :type config: StarBugMainConfig
        """
        # defaults.
        self._config = config
        self._f_name: Optional[str] = None
        self._out_dir: Optional[str] = None
        self._b_name: Optional[str] = None
        self._image:  Optional[HDUList] = None
        self._filter: str | None = None
        self._header: Header | None = None
        self._wcs: Optional[WCS] = None
        self._stage: float = 0.0
        self._detections: Optional[Table] = None
        self._n_hdu: int = -1
        self._unit: Optional[str] = None
        self._background: Optional[ImageHDU | PrimaryHDU] = None
        self._residuals: np.ndarray | None = None
        self._psf_catalogue: Optional[Table] = None
        self._source_stats: Optional[np.ndarray] = None
        self._psf: Optional[np.ndarray] = None

        # overridden configs
        self._ap_file = ap_file
        self._background_file = bkg_file
        self._verbose = verbose
        self._full_width_half_max = config.full_width_half_max

        # process options.

        ## Load the fits image
        self.load_image(f_name)

        if ap_file is not None:
            ## Load the source list if given
            self.load_ap_file(ap_file)
        if bkg_file is not None:
            self.load_bgd_file(bkg_file)


    def log(self, msg: str) -> None:
        """
        Print message if in verbose mode

        :param msg: Message to print out
        :type msg: str
        :return: None
        """
        if self._config.verbose_logs:
            printf(msg)
            sys.stdout.flush()


    def load_image(self, f_name: str | None) -> None:
        """
        Given f_name, load the image into starbug to be worked on.

        :param f_name: Filename of fits image (with any number of extensions).
            If using a non-standard HDU index, set the name or index of the
            extension with "HDU_NAME=XXX" in the parameter file.
        :type f_name: str | None
        :return: None
        """
        self._f_name = f_name
        if f_name:
            #########################################
            # Sorting out the file names and what not
            #########################################
            extension: str
            self._out_dir, self._b_name, extension = self.sort_output_names(
                f_name, self._config.output_file)

            if extension == FITS_EXTENSION:
                if os.path.exists(f_name):
                    self.log("loaded: \"%s\"\n" % f_name)
                    self._image = open(f_name)

                    ## Force assigning _nHDU
                    main_image: ImageHDU | PrimaryHDU = self.main_image
                    self._header = main_image.header

                    self.log(
                        "-> using image HDU: %d (%s)\n" % (
                            self._n_hdu, main_image.name))

                    if main_image.data is None:
                        warn("Image seems to be empty.\n")

                    if ((val := main_image.header.get(TELESCOPE)) is None
                            or (val.find(JWST)<0)):
                        warn("Telescope not JWST, "
                             "there may be undefined behaviour.\n")

                    self._filter = self._config.custom_filter
                    assert self._filter is not None
                    if ((FILTER in main_image.header) and
                        (main_image.header[FILTER] in
                         STAR_BUG_FILTERS.keys())):
                        self._filter = main_image.header[FILTER]
                        assert self._filter is not None
                        if self._full_width_half_max < 0:
                            self._full_width_half_max = (STAR_BUG_FILTERS[
                                self._filter].full_width_half_max)
                    if self._filter:
                        self.log("-> photometric band: %s\n" % self._filter)
                    else:
                        warn("Unable to determine image filter\n")

                    if DETECTOR in self.info.keys():
                        self.log(
                            "-> detector module: %s\n" % self.info[DETECTOR])
                    else:
                        warn("Unable to determine Telescope DETECTOR.\n")

                    if BUN_IT in main_image.header:
                        self._unit = main_image.header[BUN_IT]
                    else:
                        warn("Unable to determine image BUNIT.\n")

                    self._wcs = WCS(self.main_image.header)

                    ## I NEED TO DETERMINE BETTER WHAT STAGE IT IS IN
                    extension_names: List[str] = ext_names(self._image)
                    if DQ in extension_names:
                        if AREA in extension_names:
                            self._stage = 2.0
                        else:
                            self._stage = 2.5
                    elif WHT in extension_names:
                        self._stage = 3.0
                    elif CALIBRATION_LV in self.main_image.header:
                        self._stage = self.main_image.header[CALIBRATION_LV]
                    else:
                        warn("Unable to determine calibration level, "
                             "assuming stage 3\n")
                        self._stage = 3.0
                    self.log("-> pipeline stage: %d\n" % self._stage)

                else:
                    warn("fits file \"%s\" does not exist\n" % f_name)
            else:
                warn("included file must be FITS format\n")

    def load_ap_file(self, f_name: Optional[str]=None) -> None:
        """
        Load an AP_FILE to be used during photometry
        
        :param f_name: Filename for fits table containing source coordinates.
         These coordinates can be x-centroid / y-centroid, x_init / y_init,
          x_0, y_0 or RA/DEC. The latter is used if starbug gets "USE_WCS=1" 
          in the parameter file.
        :type f_name: str
        :return: None
        """
        if not f_name: 
            f_name = self._ap_file
        if os.path.exists(f_name):
            self._detections = import_table(f_name)
            if self._detections is None or self._wcs is None:
                raise Exception("could not read the ap file")

            column_names: set[str] = set(self._detections.colnames)

            self.log("loaded AP_FILE='%s'\n" % f_name)

            if self._config.use_wcs_values:
                if len(column_names & {RA, DEC}) == 2:
                    self.log("-> using RA-DEC coordinates\n")
                    try:
                        xy: Any = self._wcs.all_world2pix(
                            self._detections[RA], self._detections[DEC], 0)
                    except (NoConvergence, MemoryError, SingularMatrixError,
                            InconsistentAxisTypesError, ValueError,
                            InvalidTransformError) as e:
                        warn(f"Something went wrong converting WCS to pixels "
                             f"({e}), trying wcs_world2pix next.\n")
                        xy: Any = self._wcs.wcs_world2pix(
                            self._detections[RA], self._detections[DEC], 0)
                    if X_CENTROID in column_names: 
                        self._detections.remove_column(X_CENTROID)
                    if Y_CENTROID in column_names: 
                        self._detections.remove_column(Y_CENTROID)
                    self._detections.add_columns(
                        xy, names=[X_CENTROID, Y_CENTROID], indexes=[0, 0])
                else:
                    warn("No 'RA' or 'DEC' found in AP_FILE\n")

            elif len({X_0, Y_0} & column_names) == 2:
                self._detections.rename_columns(
                    (X_0, Y_0), (X_CENTROID, Y_CENTROID))
            elif len({X_INIT, Y_INIT} & column_names) == 2:
                self._detections.rename_columns(
                    (X_INIT, Y_INIT), (X_CENTROID, Y_CENTROID))

            if len({X_CENTROID, Y_CENTROID} & 
                   set(self._detections.colnames)) == 2:
                mask: np.ndarray = (
                    (self._detections[X_CENTROID] >= 0)
                    & (self._detections[X_CENTROID] < self.main_image.shape[1])
                    & (self._detections[Y_CENTROID] >= 0)
                    & (self._detections[Y_CENTROID] < self.main_image.shape[0])
                )

                # cant figure how to resolve this typing
                self._detections.remove_rows(~mask) # noqa
                self.log(
                    "-> loaded %d sources from AP_FILE\n" % 
                    len(self._detections))
            else:
                warn("Unable to determine physical coordinates"
                     " from detections table\n")
        else:
            p_error("AP_FILE='%s' does not exists\n" % f_name)

    def load_bgd_file(self, f_name: Optional[str]=None) -> None:
        """
        Load a BGD_FILE to be used during photometry

        :param f_name: Filename of fits image the same dimensions as the
                       main image
        :type f_name: str
        :return: None
        """
        if f_name is None:
            f_name = self._config.background_file
        if f_name is None:
            return
        if os.path.exists(f_name):
            self._background = open(f_name)[1]
            self.log("loaded BGD_FILE='%s'\n" % f_name)
        else:
            p_error("BGD_FILE='%s' does not exist\n" % f_name)

    # noinspection SpellCheckingInspection
    def load_psf(self, f_name: Optional[str]=None) -> int:
        """
        Load a PSF_FILE to be used during photometry

        :param f_name: Filename of a PSF fits image
        :type f_name: str
        :return: the status
        :rtype int
        """
        status: int = EXIT_SUCCESS
        assert self._filter is not None
        if not f_name:
            filter_struct: FilterStruct | None = (
                STAR_BUG_FILTERS.get(self._filter))
            if filter_struct:
                dt_name: str = self.info[DETECTOR]
                if dt_name == "NRCALONG":
                    dt_name = "NRCA5"
                if dt_name == "NRCBLONG":
                    dt_name = "NRCB5"
                if dt_name == "MULTIPLE":
                    if (filter_struct.instr == NIRCAM
                            and filter_struct.length == SHORT):
                        dt_name = "NRCA1"
                    elif (filter_struct.instr == NIRCAM and
                          filter_struct.length == LONG):
                        dt_name = "NRCA5"
                    elif filter_struct.instr == STAR_BUG_MIRI:
                        dt_name = ""
                if dt_name == "MIRIMAGE":
                    dt_name = ""
                f_name = "%s/%s%s.fits" % (
                    StarbugBase.get_data_path(), self._filter, dt_name)
            else:
                status = EXIT_FAIL

        if f_name is not None and os.path.exists(f_name):
            fp: HDUList = open(f_name)

            if fp[0].data is None: 
                p_error(
                    "There is a version mismatch between starbug and "
                    "webbpsf. Please reinitialise with: starbug2 --init.\n")
                quit("Fatal error, quitting\n")

            ####hmm
            self._psf = fp[0].data
            fp.close()
            self.log("loaded PSF_FILE='%s'\n" % f_name)
        else:
            p_error("PSF_FILE='%s' does not exist\n" % f_name)
            status = EXIT_FAIL
        return status

    def prepare_image_arrays(self) -> (
            Tuple[np.ndarray, np.ndarray, np.ndarray | None, np.ndarray]):
        """
        Make a copy of the original image, and prepare the other image arrays

        :return: tuple of image, error, bgd, mask
        :rtype: tuple of int /float, np.ndarray, np.ndarray or None, int
        """

        # Collect scale factor
        scale_factor: int | float
        if self.header.get(BUN_IT) == "MJy/sr":
            scale_factor = get_mj_ysr2jy_scale_factor(self.main_image)
            self.log(
                "-> converting unit from MJy/sr to Jr with factor: %e\n"
                % scale_factor)
        else:
            scale_factor = 1

        image: np.ndarray = self.main_image.data.copy() * scale_factor

        # scale by area
        extension_names: list[str] = ext_names(self._image)
        assert self._image is not None
        if AREA in extension_names:
            ## AREA distortion correction
            image *= self._image[AREA].data

        # collect and scale error
        error: np.ndarray
        if ERR in extension_names and np.shape(self._image[ERR]):
            error = self._image[ERR].data.copy() * scale_factor
        else:
            error = np.sqrt(np.abs(image))
        
        # create mask
        mask: np.ndarray
        if DQ in extension_names:
            mask = self._image[DQ].data & (DQ_DO_NOT_USE | DQ_SATURATED)
            mask = mask.astype(bool)
        else:
            mask = (np.isnan(image) | np.isnan(error))

        # collect and scale background array
        bgd: np.ndarray | None
        if self._background is not None:
            bgd = self._background.data.copy() * scale_factor
        else:
            bgd = None

        return image, error, bgd, mask

    def detect(self) -> int:
        """
        Full source detection routine. Saves the result as a table
        self._detections
        :return: status
        :rtype: int
        """
        self.log("Detecting Sources\n")
        status: int = EXIT_SUCCESS
        assert self._filter is not None
        if self.main_image:
            filter_struct: FilterStruct | None = (
                STAR_BUG_FILTERS.get(self._filter))

            full_width_half_max: float
            if self._full_width_half_max > 0:
                full_width_half_max = self._full_width_half_max
            elif filter_struct:
                full_width_half_max = filter_struct.full_width_half_max
            else:
                full_width_half_max = DEFAULT_FULL_WIDTH_HALF_MAX

            # noinspection SpellCheckingInspection
            detector: DetectionRoutine = DetectionRoutine(
                sig_src=self._config.sigma_source,
                sig_sky=self._config.sigma_sky,
                full_width_half_max=full_width_half_max,
                sharp_lo=self._config.sharp_cutoff_low,
                sharp_hi=self._config.sharp_cutoff_high,
                round_1_hi=self._config.round1_cutoff_high,
                round_2_hi=self._config.round2_cutoff_high,
                smooth_lo=self._config.smooth_low,
                smooth_hi=self._config.smooth_high,
                ricker_r=self._config.ricker_wavelet_radius,
                do_bgd_2d=self._config.do_bgd_2d,
                do_con_vl=self._config.do_convolution,
                box_size=self._config.background_box_size,
                clean_src=self._config.clean_sources,
                verbose=self._verbose)

            self._detections = detector(self.main_image.data.copy())[
                 X_CENTROID, Y_CENTROID, "sharpness", "roundness1",
                 "roundness2"]

            # check for insane states
            if self._detections is None or self._wcs is None:
                return EXIT_FAIL

            ra: np.ndarray
            dec: np.ndarray
            ra, dec = self._wcs.all_pix2world(
                self._detections[X_CENTROID], self._detections[Y_CENTROID], 0)
            self._detections.add_column( Column(ra, name=RA), index=2)
            self._detections.add_column( Column(dec, name=DEC), index=3)
            self._detections.meta=dict(self.header.items())

            # noinspection SpellCheckingInspection
            self._detections.meta.update({"ROUNTINE": "DETECT"})
            self.aperture_photometry()

        else:
            p_error("Something went wrong.\n")
            status = EXIT_FAIL
        return status


    # noinspection SpellCheckingInspection
    def aperture_photometry(self) -> int:
        """
        executes aperture photometry
        :return: 0 for success 1 for failure
        :rtype int
        """
        if self._detections is None:
            p_error("No detection source file loaded (-d file-ap.fits)\n")
            return EXIT_FAIL
        if self._image is None:
            p_error("No image provided")
            return EXIT_FAIL
        if len({"x_0", "y_0", "x_init", "y_init", X_CENTROID, Y_CENTROID} &
               set(self._detections.colnames)) < 2:
            p_error("No pixel coordinates in source file\n")
            return EXIT_FAIL
        if self._filter is None:
            p_error("no filter name")
            return EXIT_FAIL

        new_columns: tuple[str, str, str, str, str, str | None, str] = (
            "smoothness","flux","eflux","sky", "flag",
            self._filter, "e%s" % self._filter)
        self._detections.remove_columns(
            set(new_columns) & set(self._detections.colnames))


        #######################
        # APERTURE PHOTOMETRY #
        #######################
        self.log("\nRunning Aperture Photometry\n")

        image: np.ndarray
        error: np.ndarray
        mask: np.ndarray
        image, error, _, mask = self.prepare_image_arrays()

        #######################
        # Aperture Correction #
        #######################
        ap_corr_f_name: Optional[str] = None
        if _ap_corr_f_name := self._config.ap_corr_file_override:
            ap_corr_f_name = _ap_corr_f_name
        elif   self.info.get(INSTRUMENT) == NIRCAM_STRING:
            ap_corr_f_name = (
                "%s/apcorr_nircam.fits" % StarbugBase.get_data_path())
        elif self.info.get(INSTRUMENT) == "MIRI":
            ap_corr_f_name = (
                "%s/apcorr_miri.fits" % StarbugBase.get_data_path())

        if ap_corr_f_name:
            self.log("-> apcorr file: %s\n" % ap_corr_f_name)
        else:
            warn("No apcorr file available for instrument\n")

        radius: float = float(self._config.aperture_phot_radius)
        ee_frac: float = float(self._config.encircled_energy_fraction)
        sky_in: float = float(self._config.sky_annulus_inner_radius)
        sky_out: float = float(self._config.sky_annulus_outer_radius)


        if ee_frac >= 0:
            radius: float = APPhotRoutine.radius_from_enc_energy(
                self._filter, ee_frac, ap_corr_f_name)
            if radius > 0:
                self.log(
                    "-> calculating aperture radius from encircled energy\n")

        if radius <= 0:
            if (radius := self._full_width_half_max) > 0:
                self.log("-> using FWHM as aperture radius\n")
            else:
                self.log(
                    "No valid aperture radius was detected. defaulting "
                    "to default value of 2")
                radius = 2.0

        ap_corr: float = APPhotRoutine.calc_ap_corr(
            self._filter, radius, table_f_name=ap_corr_f_name,
            verbose=self._verbose)

        ##################
        # Run Photometry #
        ##################
        app_hot: APPhotRoutine = APPhotRoutine(
            radius, sky_in, sky_out, verbose=bool(self._verbose))

        dq_flags: np.ndarray | None
        if DQ in ext_names(self._image):
            dq_flags = self._image[DQ].data.copy()
        else:
            dq_flags = None
        ap_cat: Table = app_hot(
            image, self._detections, error=error, dq_flags=dq_flags,
            ap_corr=ap_corr, sig_sky=self._config.sigma_sky)


        filter_string: str = self._filter if self._filter else "mag"

        # extract magitudes
        mag: float
        mag_err: float
        mag, mag_err = flux2mag(ap_cat[FLUX], ap_cat[E_FLUX])

        # add columsn to the catalogue
        ap_cat.add_column(Column(
            mag + self._config.zero_point_magnitude, filter_string))
        ap_cat.add_column(Column(
            mag_err, "e%s" % filter_string))

        # update detections
        self._detections = hstack((self._detections, ap_cat))

        # check for insanitiy
        if self._detections is None:
            return EXIT_FAIL

        if self._config.clean_sources:
            detections_length = len(self._detections)
            if (smooth_lo := self._config.smooth_low) != "":
                self._detections.remove_rows(
                    self._detections["smoothness"] < smooth_lo)
            if (smooth_hi := self._config.smooth_high) != "":
                self._detections.remove_rows(
                    self._detections["smoothness"] > smooth_hi)
            if len(self._detections) != detections_length:
                self.log("-> removing %d sources outside SMOOTH range\n"
                         % (detections_length - len(self._detections)))

        reindex(self._detections)
        self._detections.meta[FILTER] = self._filter

        f_name = "%s/%s-ap.fits" % (self._out_dir, self._b_name)
        self.log("--> %s\n" % f_name)
        export_table(self._detections, f_name, header=self.header)

        return EXIT_SUCCESS


    def bgd_estimate(self) -> int:
        """
        Estimate the background of the active image
        Saves the result as an ImageHDU self._background
        :return: the status.
        :rtype: int
        """
        self.log("\nEstimating Diffuse Background\n")
        status: int = EXIT_SUCCESS
        assert self._filter is not None
        if self._detections:
            source_list: Table = self._detections.copy()

            filter_struct: FilterStruct | None = (
                STAR_BUG_FILTERS.get(self._filter))
            full_width_half_max: float
            if self._full_width_half_max > 0:
                full_width_half_max = self._config.full_width_half_max
            elif filter_struct:
                full_width_half_max = filter_struct.full_width_half_max
            else:
                full_width_half_max = 2.0

            # noinspection DuplicatedCode
            if X_INIT in source_list.colnames:
                source_list.rename_column(X_INIT, X_CENTROID)
            if Y_INIT in source_list.colnames:
                source_list.rename_column(Y_INIT, Y_CENTROID)
            if X_DET in source_list.colnames:
                source_list.rename_column(X_DET, X_CENTROID)
            if Y_DET in source_list.colnames:
                source_list.rename_column(Y_DET, Y_CENTROID)
            if FLUX_DET in source_list.colnames:
                source_list.rename_column(FLUX_DET, FLUX)
            mask: np.ndarray = ~(np.isnan(source_list[X_CENTROID])
                     | np.isnan(source_list[Y_CENTROID]))


            bgd: BackGroundEstimateRoutine = BackGroundEstimateRoutine(
                source_list[mask],
                box_size=self._config.background_box_size,
                full_width_half_max=full_width_half_max,
                sig_sky=self._config.sigma_sky,
                bgd_r=self._config.bgd_radius,
                profile_scale=self._config.profile_scaling_factor,
                profile_slope=self._config.profile_slope,
                verbose=self._verbose)
            header: Header = self.header

            # check for insanity
            if self._wcs is None:
                return EXIT_FAIL

            header.update(self._wcs.to_header())

            # get image data
            image_data = bgd(
                self.main_image.data.copy(),
                output=self._config.bgd_check_file)
            assert image_data is not None

            self._background = ImageHDU(
                data=image_data.background,
                header=header)

            # check for insanity
            if self._background is None:
                return EXIT_FAIL

            f_name = "%s/%s-bgd.fits"%(self._out_dir, self._b_name)
            self.log("--> %s\n" % f_name)
            self._background.writeto(f_name, overwrite=True)

        else:
            p_error("unable to estimate background, no source list loaded\n")
            status = EXIT_FAIL
        return status



    def bgd_subtraction(self) -> int:
        """
        Internally subtract a background array from an image array
        :return: 0 for success, 1 otherwise
        :rtype int.
        """
        self.log("Subtracting Background\n")

        if self._background is None or self._wcs is None:
            p_error("No background array loaded (-b file-bgd.fits)\n")
            return EXIT_FAIL
        array: np.ndarray = self.main_image.data - self._background.data
        self._residuals = array

        assert self._image is not None
        self._image[self._n_hdu].data = array
        header: Header = self.header
        header.update(self._wcs.to_header())

        # having to cast to any as the ImageHUD expects an 'array.pyi' and
        # not a ndarray. Note this is being used as a glorified writer
        ImageHDU(
            data=cast(Any, self._residuals), name="RES",
            header=header).writeto(
                "%s/%s-res.fits" % (self._out_dir, self._b_name),
            overwrite=True)
        return EXIT_SUCCESS

    # noinspection SpellCheckingInspection
    def photometry_routine(self) -> int:
        """
        Full photometry routine
        Saves the result as a table self._psf_catalogue,
        Additionally it appends a residual Image onto the
        self._residuals HDUList

        :return: 0 for success, 1 otherwise
        :rtype int
        """
        if self._filter is None or self._wcs is None:
            return EXIT_FAIL

        if self.main_image:
            self.log("\nRunning PSF Photometry\n")

            # lock the types.
            image: np.ndarray
            error: np.ndarray
            bgd: np.ndarray | None
            mask: np.ndarray
            image, error, bgd, mask = self.prepare_image_arrays()

            if bgd is None:
                clipped_median: float
                _, clipped_median, _ = (
                    sigma_clipped_stats(image, sigma=self._config.sigma_sky))
                bgd = np.ones(self.main_image.shape) * clipped_median
                self.log(
                    "-> no background file loaded, measuring sigma "
                    "clipped median\n")
            assert bgd is not None

            ###################################
            # Collect relevant files and data #
            ###################################

            if self._detections is None:
                p_error("unable to run photometry: no source list loaded\n")
                return EXIT_FAIL

            if (self._psf is None
                and self.load_psf(
                    os.path.expandvars(self._config.psf_file_override))):
                p_error("unable to run photometry: no PSF loaded\n")
                return EXIT_FAIL

            psf_mask: np.ndarray = ~np.isfinite(self._psf)
            if psf_mask.sum():
                assert self._psf is not None
                self._psf[psf_mask] = 0
                self.log("-> masking INF pixels in PSF_FILE\n")

            psf_model: ImagePSF = ImagePSF(data=self._psf)
            size: int
            if self._config.psf_fit_size > 0:
                size = self._config.psf_fit_size
            else:
                size = psf_model.shape[0]
            if not size % 2:
                size -= 1
            self.log("-> psf size: %d\n" % size)

            #########################
            # Sort out Init guesses #
            #########################
            app_hot_r: float = float(self._config.aperture_phot_radius)
            if not app_hot_r or app_hot_r <= 0:
                app_hot_r = 3.0

            init_guesses: Table = self._detections.copy()

            # noinspection DuplicatedCode
            if X_CENTROID in init_guesses.colnames:
                init_guesses.rename_column(X_CENTROID, X_INIT)
            if Y_CENTROID in init_guesses.colnames:
                init_guesses.rename_column(Y_CENTROID, Y_INIT)
            if X_DET in init_guesses.colnames:
                init_guesses.rename_column(X_DET, X_INIT)
            if Y_DET in init_guesses.colnames:
                init_guesses.rename_column(Y_DET, Y_INIT)

            init_guesses = init_guesses[ init_guesses[X_INIT] >=0 ]
            init_guesses = init_guesses[ init_guesses[Y_INIT] >=0 ]
            init_guesses = init_guesses[
                init_guesses[X_INIT] < self.main_image.header[NAXIS1]]
            init_guesses=init_guesses[
                init_guesses[Y_INIT] < self.main_image.header[NAXIS2]]

            ######
            # Allow tables that don't have the correct columns through
            ######
            required: List[str] = [X_INIT, Y_INIT, FLUX, self._filter, FLAG]
            for notfound in  set(required) - set(init_guesses.colnames):
                dtype = np.uint16 if notfound == FLAG else float
                init_guesses.add_column(
                    Column(np.zeros(len(init_guesses)),
                           name=notfound, dtype=dtype))

            init_guesses = init_guesses[required]
            init_guesses.remove_column(FLUX)
            init_guesses.rename_column(self._filter, "ap_%s" % self._filter)
            
            ###########
            # Run Fit #
            ###########

            min_separation: float = self._config.critical_separation
            if not min_separation:
                min_separation = min(5.0, 2.5 * self._full_width_half_max)

            if self._config.force_centroid_position:
                phot: PSFPhotRoutine = PSFPhotRoutine(
                    psf_model, size, min_separation=min_separation,
                    app_hot_r=app_hot_r, background=bgd, force_fit=1,
                    verbose=self._verbose)
                psf_cat: Table = phot(
                    image, init_params=init_guesses, error=error,
                    mask=mask)
                psf_cat[FLAG] |= SRC_FIX

            else:
                phot: PSFPhotRoutine = PSFPhotRoutine(
                    psf_model, size, min_separation=min_separation,
                    app_hot_r=app_hot_r, background=bgd, force_fit=0,
                    verbose=self._verbose)
                psf_cat: Table = phot(
                    image, init_params=init_guesses, error=error,
                    mask=mask)

                if not psf_cat:
                    return EXIT_FAIL

                ##################################
                # Setting position max variation #
                ##################################
                max_y_dev: float
                unit: int
                max_y_dev, unit = parse_unit(self._config.max_xy_deviation)
                if unit is not None:
                    if unit == DEG:
                        max_y_dev *= 60
                        unit = ARCMIN
                    if unit == ARCMIN:
                        max_y_dev *= 60
                        unit = ARCSEC
                    if unit == ARCSEC:
                        if not self.header.get(PIXAR_A2):
                            warn(
                                "MAX_XYDEV is units arcseconds, but starbug "
                                "cannot locate a pixel scale in the header."
                                " Please use syntax MAX_XYDEV=%sp to set "
                                "change to pixels\n" % max_y_dev)
                        else:
                            max_y_dev /= np.sqrt(self.header.get(PIXAR_A2))

                if max_y_dev > 0:
                    self.log(
                        "-> position fit threshold: %.2gpix\n" % max_y_dev)
                    phot: PSFPhotRoutine = PSFPhotRoutine(
                        psf_model, size, min_separation=min_separation,
                        app_hot_r=app_hot_r, background=bgd, force_fit=1,
                        verbose=self._verbose)
                    ii: bool = psf_cat[XY_DEV] > max_y_dev
                    fixed_centres: Table = psf_cat[ii][
                        [X_INIT, Y_INIT, "ap_%s" % self._filter, FLAG]]
                    if len(fixed_centres):
                        self.log("-> forcing positions for deviant sources\n")
                        fixed_cat: Table = phot(
                            image, init_params=fixed_centres,
                            error=error, mask=mask)
                        fixed_cat[FLAG] |= SRC_FIX
                        psf_cat.remove_rows(ii)
                        psf_cat = vstack((psf_cat, fixed_cat))
                    else:
                        self.log("-> no deviant sources\n")
            ra: np.ndarray
            dec: np.ndarray
            ra, dec = self._wcs.all_pix2world(
                psf_cat[X_FIT], psf_cat[Y_FIT], 0)
            psf_cat.add_column(Column(ra, name=RA), index=2)
            psf_cat.add_column(Column(dec, name=DEC), index=3)

            mag: float
            mag_error: float
            mag, mag_err = flux2mag(psf_cat[FLUX],psf_cat[E_FLUX])

            filter_string: str = self._filter if self._filter else MAG
            psf_cat.add_column(
                mag + self._config.zero_point_magnitude, name=filter_string)
            psf_cat.add_column(mag_err, name="e%s" % filter_string)
            self._psf_catalogue = psf_cat

            # verify catalogue isnt none
            assert self._psf_catalogue is not None

            self._psf_catalogue.meta = dict(self.header.items())
            self._psf_catalogue.meta[AP_FILE]=self._ap_file
            self._psf_catalogue.meta[BGD_FILE]=self._background_file

            reindex(self._psf_catalogue)

            file_name: str = (
                "%s/%s-psf.fits" % (self._out_dir, self._b_name))
            self.log("--> %s\n" % file_name)
            BinTableHDU(
                data=self._psf_catalogue,
                header=self.header).writeto(file_name, overwrite=True)

            ##################
            # Residual Image #
            ##################

            if self._config.generate_residual_image:
                self.log("-> generating residual\n")
                _tmp: Table = psf_cat[X_FIT, Y_FIT, FLUX].copy()
                _tmp.rename_columns( (X_FIT, Y_FIT), (X_0, Y_0))
                stars: np.ndarray = make_model_image(
                    image.shape, psf_model, _tmp, model_shape=(size,size))
                residual: np.ndarray = image - (bgd + stars)
                self._residuals = (
                    residual / get_mj_ysr2jy_scale_factor(self.main_image))
                header: Header = self.header
                header.update(self._wcs.to_header())
                ImageHDU(
                    data=cast(Any, self._residuals),
                    name="RES", header=header).writeto(
                        "%s/%s-res.fits" % (self._out_dir, self._b_name),
                        overwrite=True)
            return EXIT_SUCCESS
        return EXIT_FAIL

    def source_geometry(self) -> None:
        """
        Calculate source geometry stats for a given image and source list
        :return: None
        """
        if self._detections is None:
            p_error("No source file loaded\n")
            return

        if self._filter is None:
            p_error("no filter string provided\n")
            return

        self.log("Running Source Geometry\n")
        slist: Table = self._filter_detections()

        sp: SourceProperties = SourceProperties(
            self.main_image.data, slist,
            verbose=self._verbose)
        stat: Table = sp(
            full_width_half_max=STAR_BUG_FILTERS[
                self._filter].full_width_half_max,
            do_crowd=self._config.calculate_crowding_metric)

        self._source_stats = hstack((slist, stat))
        f_name: str = "%s/%s-stat.fits" % (self._out_dir, self._b_name)
        self.log("--> %s\n" % f_name)
        reindex(Table(self._source_stats))
        BinTableHDU(
            data=self._source_stats, header=self.header).writeto(
                f_name, overwrite=True)

    # noinspection SpellCheckingInspection
    def verify(self) -> int:
        """
        This simple function verifies that everything necessary has been
        loaded properly

        :return: int where 0 on success, 1 on fail
        :rtype int
        """

        status: int = EXIT_SUCCESS
        
        self.log("Checking internal systems..\n")

        if not self._filter:
            warn("No FILTER set, please set in parameter file or "
                 "use \"-s FILTER=XXX\"\n")
            status = EXIT_FAIL

        d_name: str = os.path.expandvars(StarbugBase.get_data_path())
        if not os.path.exists(d_name):
            warn("Unable to locate STARBUG_DATDIR='%s'\n" % d_name)

        if self._out_dir is not None and not os.path.exists(self._out_dir):
            warn("Unable to locate OUTPUT='%s'\n" % self._out_dir)
            status = EXIT_FAIL
        
        if self._image is None or self.main_image.data is None:
            warn("Image did not load correctly\n")
            status = EXIT_FAIL

        if self._ap_file and self._detections is not None:
            test = self._filter_detections()
            if not len(test):
                warn("Detection file empty or no sources overlap the image.\n")
                status = EXIT_FAIL

        return status


    def _filter_detections(self) -> Table:
        """
        filters the detections based on some fixed constraints.
        :return: the filtered detections
        """
        assert self._detections is not None
        detections: Table = self._detections[[X_CENTROID,Y_CENTROID]].copy()
        detections = detections[ detections[X_CENTROID]>=0 ]
        detections = detections[ detections[Y_CENTROID]>=0 ]
        detections = detections[
            detections[X_CENTROID] < self.main_image.header[NAXIS1]]
        return detections[detections[Y_CENTROID] <
                self.main_image.header[NAXIS2]]

    def __getstate__(self) -> dict[str, Any]:
        """
        extracts the inner state of this class. deleting image or/and
         background if it's there.
        :return: the internal state with those bits filtered away
        """
        assert self._image is not None
        self._image.close()
        state:  dict[str, Any] = self.__dict__.copy()
        if "_image" in state:
            ##Sorry but we cant have that
            del state["_image"]
            ## This currently doesnt get reloaded
        if "_background" in state:
            del state["_background"]
            
        return state

    def __setstate__(self, state) -> None:
        self.__dict__.update(state)
        v:  int = int(self._verbose)
        self._verbose = 0
        self.load_image(self._f_name)
        self._verbose = v

    @property
    def header(self) -> Header:
        """
        Construct relevant base header information for routine products

        :return:  Header file containing a series of relevant information
        :rtype: Header
        """
        head: Dict[str, str | float] = {
            STAR_BUG: get_version(),
            CALIBRATION_LV: self._stage
        }

        if self._filter:
            head[FILTER] = self._filter

        # add the basic params
        for fits_key, (property_name, _) in (
                StarBugMainConfig.MAIN_PARAM_FILE_MAP.items()):
            value = getattr(self._config, property_name)
            if value is None:
                head[fits_key] = ""
            else:
                head[fits_key] = value

        # add the changed ones
        head[AP_FILE] = self._ap_file
        head[BGD_FILE] = self._background_file
        head[VERBOSE_TAG] = self._verbose

        # add info
        head.update(self.info)
        return collapse_header(head)


    @property
    def info(self) -> dict[str, str]:
        """
        Get some useful information from the image header file.

        :return: extracted keys and elements from the image header.
        :rtype: dict of str, to str.
        """
        out: dict[str, str] = {}
        keys: list[str] = [
            FILTER, DETECTOR, TELESCOPE, INSTRUMENT, BUN_IT, PIXAR_A2,
            PIXAR_SR]
        if self._image:
            for hdu in self._image:
                out.update(
                    { (key,hdu.header[key]) for key in keys
                      if key in hdu.header})
        return out

    @property
    def main_image(self) -> ImageHDU | PrimaryHDU:
        # noinspection SpellCheckingInspection
        """
        automagically find the main image array to use
        Order of importance is:
        > self._nHDU (if set)
        > param[ HDUNAME ]
        > SCI, BGD, RES
        > first ImageHDU
        > first ImageHDU
        > image[0]

        :return: the main image array.
        :rtype: HDUList
        """
        assert self._image is not None
        if self._n_hdu >= 0:
            return self._image[self._n_hdu]
        e_names: list[str] = ext_names(self._image)

        ## HDU_NAME in param file
        n: str = str(self._config.hdu_name)
        if n and n in e_names:
            self._n_hdu = e_names.index(n)
            return self._image[n]

        ##index?
        if isinstance(n, (int, float, np.number)):
            self._n_hdu = int(n)
            return self._image[self._n_hdu]

        ## SCI, BGD, RES (common names)
        for name in (SCI, BGD, RES):
            name: str
            if name in e_names:
                self._n_hdu = e_names.index(name)
                return self._image[name]

        ## First ImageHDU
        #ABS ARE WE SURE WE WANT TO LOOK FOR A INDEX WITH A ENUMERATE INDEX?
        assert self._image is not None
        for index, hdu in enumerate(self._image):
            index: int
            hdu: ImageHDU | PrimaryHDU | BinTableHDU
            if isinstance(hdu, ImageHDU):
                self._n_hdu = index
                return hdu

        self._n_hdu = 0
        return self._image[0]

    @property
    def filter(self) -> str | None:
        return self._filter

    @property
    def n_hdu(self) -> int:
        return self._n_hdu

    @property
    def image(self) -> HDUList | None:
        return self._image

    @image.setter
    def image(self, new_image: HDUList) -> None:
        self._image = new_image

    @property
    def psf_catalogue(self) -> Table | None:
        return self._psf_catalogue

    @property
    def psf(self) -> np.ndarray | None:
        return self._psf

    @property
    def f_name(self) -> str | None:
        return self._f_name

    @property
    def detections(self) -> Table | None:
        return self._detections

    @detections.setter
    def detections(self, new_detections: Table) -> None:
        self._detections = new_detections

    @property
    def out_dir(self) -> str | None:
        return self._out_dir
