import os
import sys
from os import getenv

from astropy.wcs import WCS, NoConvergence
import numpy as np
from astropy.io import fits
from astropy.table import hstack, Column, vstack
from astropy.stats import sigma_clipped_stats
from photutils.datasets import make_model_image
from photutils.psf import FittableImageModel

from starbug2.constants import (
    VERBOSE, FILTER, STAR_BUG, CALIBRATION_LV, DETECTOR, TELESCOPE,
    INSTRUMENT, BUN_IT, PIXAR_A2, PIXAR_SR, HDU_NAME, SCI, BGD, RES,
    VERBOSE_TAG, AP_FILE, BGD_FILE, OUTPUT, FITS_EXTENSION, JWST, FWHM, DQ,
    AREA, WHT, USE_WCS, RA, DEC, X_CENTROID, Y_CENTROID, SHORT, LONG, NIRCAM,
    MIRI, SRC_FIX, CRIT_SEP, FORCE_POS, DEG, ARCMIN, ARCSEC, MAX_XY_DEV,
    DQ_DO_NOT_USE, DQ_SATURATED, NAXIS1, NAXIS2, CALC_CROWD, ERR,
    EXIT_SUCCESS, EXIT_FAIL, APCORR_FILE, APPHOT_R, ENCENERGY, SKY_RIN,
    SKY_ROUT, SIGSKY, ZP_MAG, CLEANSRC, QUIETMODE, BOX_SIZE, BGD_R,
    PROF_SCALE, PROF_SLOPE, BGD_CHECKFILE, PSF_FILE, PSF_SIZE, GEN_RESIDUAL,
    NIRCAM_STRING, STARBUG_DATA_DIR)
from starbug2.filters import filters
from starbug2.param import load_params, load_default_params
from starbug2.routines.app_hot_routine import APPhotRoutine
from starbug2.routines.background_estimate_routine import (
    BackGroundEstimateRoutine)
from starbug2.routines.detection_routines import DetectionRoutine
from starbug2.routines.psf_phot_routine import PSFPhotRoutine
from starbug2.routines.source_properties import SourceProperties
from starbug2.utils import (
    collapse_header, parse_unit, get_version, ext_names, printf,
    split_file_name, p_error, warn, import_table, get_mj_ysr2jy_scale_factor,
    flux2mag, reindex, export_table)


class StarbugBase(object):
    """
    StarbugBase is the overall container for the photometry package. It holds
    the active image, the parameter file and the output images/tables.
    It is self-contained enough to simply run "photometry" and everything
    should just take care of itself from there on.
    """

    @staticmethod
    def get_data_path():
        """
        returns the data path.
        :return: the data path
        """
        env_path = getenv(STARBUG_DATA_DIR)
        return (env_path if env_path else
                "%s/.local/share/starbug" % (getenv("HOME")))

    @staticmethod
    def sort_output_names(f_name, param_output=None):
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

        out_dir = ""
        b_name = ""
        extension = ""
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

    def __init__(self, f_name, p_file=None, options=None):
        """
         Star bug init.

        :param f_name: FITS image file name
        :param p_file: parameter file name
        :param options: extra options to load into starbug
        """
        # defaults.
        self._f_name = None
        self._out_dir = None
        self._b_name = None
        self._image = None
        self._filter = None
        self._header = None
        self._info = None
        self._wcs = None
        self._stage = 0
        self._detections = None
        self._n_hdu = -1
        self._unit = None
        self._background = None
        self._residuals = None
        self._psf_catalogue = None
        self._source_stats = None
        self._psf = None

        # process options.
        if options is None:
            options = {}

        if not p_file:
            if os.path.exists("starbug.param"):
                p_file = "starbug.param"
            else:
                p_file = None
        self._options = load_params(p_file)
        self._options.update(options)

        ## Load the fits image
        self.load_image(f_name)

        if self._options[AP_FILE]:
            ## Load the source list if given
            self.load_ap_file()
        if self._options[BGD_FILE]:
            self.load_bgd_file()

    @property
    def header(self):
        """
        Construct relevant base header information for routine products

        :return:  Header file containing a series of relevant information
        :rtype: fits.Header
        """
        head = {
            STAR_BUG: get_version(),
            CALIBRATION_LV: self._stage
        }

        if self._filter:
            head[FILTER] = self._filter
        head.update(self._options)
        head.update(self._info)
        return collapse_header(head)


    @property
    def info(self):
        """
        Get some useful information from the image header file.

        :return: extracted keys and elements from the image header.
        :rtype: dict of str, to str.
        """
        out = {}
        keys = (FILTER, DETECTOR, TELESCOPE, INSTRUMENT,
                BUN_IT, PIXAR_A2, PIXAR_SR)
        if self._image:
            for hdu in self._image:
                out.update(
                    { (key,hdu.header[key]) for key in keys
                      if key in hdu.header})
        return out

    @property
    def image(self):
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

        if self._n_hdu >= 0:
            return self._image[self._n_hdu]
        e_names = ext_names(self._image)

        ## HDU_NAME in param file
        n = self._options[HDU_NAME]
        if n and n in e_names:
            self._n_hdu = e_names.index(n)
            return self._image[n]

        ##index?
        if isinstance(n, (int, float, np.number)):
            self._n_hdu = int(n)
            return self._image[self._n_hdu]

        ## SCI, BGD, RES (common names)
        for name in (SCI, BGD, RES):
            if name in e_names:
                self._n_hdu = e_names.index(name)
                return self._image[name]

        ## First ImageHDU
        #ABS ARE WE SURE WE WANT TO LOOK FOR A INDEX WITH A ENUMERATE INDEX?
        for index, hdu in enumerate(self._image):
            if isinstance(hdu, fits.ImageHDU):
                self._n_hdu = e_names.index(index)
                return hdu

        self._n_hdu = 0
        return self._image[0]

    def log(self, msg):
        """
        Print message if in verbose mode

        :param msg: Message to print out
        :type msg: str
        :return: None
        """
        if self._options[VERBOSE_TAG]:
            printf(msg)
            sys.stdout.flush()


    def load_image(self, f_name):
        """
        Given f_name, load the image into starbug to be worked on.

        :param f_name: Filename of fits image (with any number of extensions).
            If using a non-standard HDU index, set the name or index of the
            extension with "HDU_NAME=XXX" in the parameter file.
        :type f_name: str
        :return: None
        """
        self._f_name = f_name
        if f_name:
            #########################################
            # Sorting out the file names and what not
            #########################################
            self._out_dir, self._b_name, extension = self.sort_output_names(
                f_name, self._options.get(OUTPUT))

            if extension == FITS_EXTENSION:
                if os.path.exists(f_name):
                    self.log("loaded: \"%s\"\n" % f_name)
                    self._image = fits.open(f_name)

                    # ABS WTF
                    _ = self.image ## Force assigning _nHDU

                    self.log(
                        "-> using image HDU: %d (%s)\n" % (
                            self._n_hdu, self._image.name))

                    if self._image.data is None:
                        warn("Image seems to be empty.\n")

                    if ((val := self._header.get(TELESCOPE)) is None
                            or (val.find(JWST)<0)):
                        warn("Telescope not JWST, "
                             "there may be undefined behaviour.\n")

                    self._filter = self._options.get(FILTER)
                    if ((FILTER in self._header) and
                            (self._header[FILTER] in filters.keys())):
                        self._filter = self._header[FILTER]
                        if self._options[FWHM] < 0:
                            self._options[FWHM ] = filters[self._filter].pFWHM
                    if self._filter:
                        self.log("-> photometric band: %s\n" % self._filter)
                    else:
                        warn("Unable to determine image filter\n")

                    if DETECTOR in self._info.keys():
                        self.log(
                            "-> detector module: %s\n" % self._info[DETECTOR])
                    else:
                        warn("Unable to determine Telescope DETECTOR.\n")

                    if BUN_IT in self._image.header:
                        self._unit = self._image.header[BUN_IT]
                    else:
                        warn("Unable to determine image BUNIT.\n")

                    self._wcs = WCS(self.image.header)

                    ## I NEED TO DETERMINE BETTER WHAT STAGE IT IS IN
                    extension_names = ext_names(self._image)
                    if DQ in extension_names:
                        if AREA in extension_names:
                            self._stage = 2
                        else:
                            self._stage = 2.5
                    elif WHT in extension_names:
                        self._stage = 3
                    elif CALIBRATION_LV in self.image.header:
                        self._stage = self.image.header[CALIBRATION_LV]
                    else:
                        warn("Unable to determine calibration level, "
                             "assuming stage 3\n")
                        self._stage = 3
                    self.log("-> pipeline stage: %d\n" % self._stage)

                else:
                    warn("fits file \"%s\" does not exist\n" % f_name)
            else:
                warn("included file must be FITS format\n")

    def load_ap_file(self, f_name=None):
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
            f_name = self._options[AP_FILE]
        if os.path.exists(f_name):
            self._detections = import_table(f_name)
            column_names = set(self._detections.col_names)

            self.log("loaded AP_FILE='%s'\n" % f_name)

            if self._options.get(USE_WCS):
                if len(column_names & {RA, DEC}) == 2:
                    self.log("-> using RA-DEC coordinates\n")
                    try:
                        xy = self._wcs.all_world2pix(
                            self._detections[RA], self._detections[DEC], 0)
                    except (NoConvergence, Exception) as e:
                        warn(f"Something went wrong converting WCS to pixels "
                             f"({e}), trying wcs_world2pix next.\n")
                        xy = self._wcs.wcs_world2pix(
                            self._detections[RA], self._detections[DEC], 0)
                    if X_CENTROID in column_names: 
                        self._detections.remove_column(X_CENTROID)
                    if Y_CENTROID in column_names: 
                        self._detections.remove_column(Y_CENTROID)
                    self._detections.add_columns(
                        xy, names=(X_CENTROID, Y_CENTROID), indexes=[0, 0])
                else:
                    warn("No 'RA' or 'DEC' found in AP_FILE\n")

            elif len({"x_0", "y_0"} & column_names) == 2:
                self._detections.rename_columns(
                    ("x_0", "y_0"), (X_CENTROID, Y_CENTROID))
            elif len({"x_init", "y_init"} & column_names) == 2:
                self._detections.rename_columns(
                    ("x_init", "y_init"), (X_CENTROID, Y_CENTROID))

            if len({X_CENTROID, Y_CENTROID} & 
                   set(self._detections.col_names)) == 2:
                mask = (
                    (self._detections[X_CENTROID] >= 0)
                    & (self._detections[X_CENTROID] < self._image.shape[1])
                    & (self._detections[Y_CENTROID] >= 0)
                    & (self._detections[Y_CENTROID] < self._image.shape[0])
                )
                self._detections.remove_rows(~mask)
                self.log(
                    "-> loaded %d sources from AP_FILE\n" % 
                    len(self._detections))
            else:
                warn("Unable to determine physical coordinates"
                     " from detections table\n")
        else: p_error("AP_FILE='%s' does not exists\n" % f_name)

    def load_bgd_file(self, f_name=None):
        """
        Load a BGD_FILE to be used during photometry

        :param f_name: Filename of fits image the same dimensions as the
                       main image
        :type f_name: str
        :return:
        """
        if not f_name: 
            f_name = self._options[BGD_FILE]
        if os.path.exists(f_name):
            self._background = fits.open(f_name)[1]
            self.log("loaded BGD_FILE='%s'\n" % f_name)
        else: p_error("BGD_FILE='%s' does not exist\n" % f_name)

    # noinspection SpellCheckingInspection
    def load_psf(self, f_name=None):
        """
        Load a PSF_FILE to be used during photometry

        :param f_name: Filename of a PSF fits image
        :type f_name: str
        :return: the status
        :rtype int
        """
        status = 0
        if not f_name:
            filter_string = filters.get(self._filter)
            if filter_string:
                dt_name = self._info[DETECTOR]
                if dt_name == "NRCALONG":
                    dt_name = "NRCA5"
                if dt_name == "NRCBLONG":
                    dt_name = "NRCB5"
                if dt_name == "MULTIPLE":
                    if (filter_string.instr == NIRCAM
                            and filter_string.length == SHORT):
                        dt_name = "NRCA1"
                    elif (filter_string.instr == NIRCAM and
                          filter_string.length == LONG):
                        dt_name = "NRCA5"
                    elif filter_string.instr == MIRI:
                        dt_name = ""
                if dt_name == "MIRIMAGE":
                    dt_name = ""
                f_name = "%s/%s%s.fits" % (
                    StarbugBase.get_data_path(), self._filter, dt_name)
            else:
                status = 1
        if os.path.exists(f_name):
            fp = fits.open(f_name)

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
            status = 1
        return status

    def prepare_image_arrays(self):
        """
        Make a copy of the original image, and prepare the other image arrays

        :return: tuple of image, error, bgd, mask
        :rtype: tuple of int /float, np.array, np.array or None, int
        """

        # Collect scale factor
        if self.header.get(BUN_IT) == "MJy/sr":
            scale_factor = get_mj_ysr2jy_scale_factor(self.image)
            self.log(
                "-> converting unit from MJy/sr to Jr with factor: %e\n"
                % scale_factor)
        else:
            scale_factor = 1

        image = self.image.data.copy() * scale_factor

        # scale by area
        extension_names = ext_names(self._image)
        if AREA in extension_names:
            ## AREA distortion correction
            image *= self._image[AREA].data

        # collect and scale error
        if ERR in extension_names and np.shape(self._image[ERR]):
            error = self._image[ERR].data.copy() * scale_factor
        else:
            error = np.sqrt(np.abs(image))
        
        # create mask
        if DQ in extension_names:
            mask = self._image[DQ].data & (DQ_DO_NOT_USE | DQ_SATURATED)
            mask = mask.astype(bool)
        else:
            mask = (np.isnan(image) | np.isnan(error))

        # collect and scale background array
        if self._background is not None:
            bgd = self._background.data.copy() * scale_factor
        else:
            bgd = None

        return image, error, bgd, mask

    def detect(self):
        """
        Full source detection routine. Saves the result as a table
        self._detections
        :return: status
        :rtype: int
        """
        self.log("Detecting Sources\n")
        status = 0
        if self.image:
            filter_map = filters.get(self._filter)
            if self._options[FWHM] > 0:
                full_width_half_max = self._options[FWHM]
            elif filter_map:
                full_width_half_max = filter_map.pFWHM
            else:
                full_width_half_max = 2

            # noinspection SpellCheckingInspection
            detector=DetectionRoutine(
                sig_src=self._options["SIGSRC"],
                sig_sky=self._options["SIGSKY"],
                full_width_half_max=full_width_half_max,
                sharp_lo=self._options["SHARP_LO"],
                sharp_hi=self._options["SHARP_HI"],
                round_1_hi=self._options["ROUND1_HI"],
                round_2_hi=self._options["ROUND2_HI"],
                smooth_lo=self._options["SMOOTH_LO"],
                smooth_hi=self._options["SMOOTH_HI"],
                ricker_r=self._options["RICKER_R"],
                do_bgd_2d=self._options["DOBGD2D"],
                do_con_vl=self._options["DOCONVL"],
                box_size=int(self._options["BOX_SIZE"]),
                clean_src=self._options["CLEANSRC"],
                verbose=self._options["VERBOSE"])

            self._detections = detector(self.image.data.copy())[
                 X_CENTROID, Y_CENTROID, "sharpness", "roundness1",
                 "roundness2"]

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
            status=1
        return status


    # noinspection SpellCheckingInspection
    def aperture_photometry(self):
        """
        executes aperture photometry
        :return: 0 for success 1 for failure
        :rtype int
        """
        if self._detections is None:
            p_error("No detection source file loaded (-d file-ap.fits)\n")
            return EXIT_FAIL
        if len({"x_0", "y_0", "x_init", "y_init", X_CENTROID, Y_CENTROID} &
               set(self._detections.col_names)) < 2:
            p_error("No pixel coordinates in source file\n")
            return EXIT_FAIL

        new_columns = (
            "smoothness","flux","eflux","sky", "flag",
            self._filter, "e%s" % self._filter)
        self._detections.remove_columns(
            set(new_columns) & set(self._detections.col_names))


        #######################
        # APERTURE PHOTOMETRY #
        #######################
        self.log("\nRunning Aperture Photometry\n")

        image, error, _, mask = self.prepare_image_arrays()

        #######################
        # Aperture Correction #
        #######################
        ap_corr_f_name=None
        if _ap_corr_f_name := self._options.get(APCORR_FILE):
            ap_corr_f_name = _ap_corr_f_name
        elif   self._info.get(INSTRUMENT) == NIRCAM_STRING:
            ap_corr_f_name = (
                "%s/apcorr_nircam.fits" % StarbugBase.get_data_path())
        elif self._info.get(INSTRUMENT) == "MIRI":
            ap_corr_f_name = (
                "%s/apcorr_miri.fits" % StarbugBase.get_data_path())

        if ap_corr_f_name:
            self.log("-> apcorr file: %s\n" % ap_corr_f_name)
        else:
            warn("No apcorr file available for instrument\n")

        radius = self._options[APPHOT_R]
        ee_frac = self._options[ENCENERGY]
        sky_in = self._options[SKY_RIN]
        sky_out = self._options[SKY_ROUT]

        if ee_frac >= 0:
            radius = APPhotRoutine.radius_from_enc_energy(
                self._filter, ee_frac, ap_corr_f_name)
            if radius > 0:
                self.log(
                    "-> calculating aperture radius from encircled energy\n")

        if radius <= 0:
            if (radius := self._options[FWHM]) > 0:
                self.log("-> using FWHM as aperture radius\n")
            else:
                radius = 2

        ap_corr = APPhotRoutine.calc_ap_corr(
            self._filter, radius, table_f_name=ap_corr_f_name,
            verbose=self._options[VERBOSE])

        ##################
        # Run Photometry #
        ##################
        app_hot = APPhotRoutine(
            radius, sky_in, sky_out, verbose=self._options[VERBOSE])

        if DQ in ext_names(self._image):
            dq_flags = self._image[DQ].data.copy()
        else:
            dq_flags = None
        ap_cat = app_hot(
            image, self._detections, error=error, dqflags=dq_flags,
            apcorr=ap_corr, sig_sky=self._options[SIGSKY])


        filter_string = self._filter if self._filter else "mag"
        mag, mag_err = flux2mag(ap_cat["flux"], ap_cat["eflux"])
        ap_cat.add_column(Column(
            mag + self._options.get(ZP_MAG), filter_string))
        ap_cat.add_column(Column(
            mag_err, "e%s" % filter_string))
        self._detections = hstack((self._detections,ap_cat))

        if self._options.get(CLEANSRC):
            detections_length = len(self._detections)
            if (smooth_lo := self._options.get("SMOOTH_LO")) != "":
                self._detections.remove_rows(
                    self._detections["smoothness"] < smooth_lo)
            if (smooth_hi := self._options.get("SMOOTH_HI")) != "":
                self._detections.remove_rows(
                    self._detections["smoothness"] > smooth_hi)
            if len(self._detections) != detections_length:
                self.log("-> removing %d sources outside SMOOTH range\n"
                         % (detections_length - len(self._detections)))

        reindex(self._detections)
        self._detections.meta[FILTER] = self._filter

        if not self._options.get(QUIETMODE):
            f_name = "%s/%s-ap.fits" % (self._out_dir, self._b_name)
            self.log("--> %s\n" % f_name)
            export_table(self._detections, f_name, header=self.header)

        return EXIT_SUCCESS


    def bgd_estimate(self):
        """
        Estimate the background of the active image
        Saves the result as an ImageHDU self._background

        ABS: do we need to return this if its always 1? should it not be 0?
        :return: the status. which seems to always be 1.
        """
        self.log("\nEstimating Diffuse Background\n")
        status = 1
        if self._detections:
            source_list = self._detections.copy()

            _f = filters.get(self._filter)
            if self._options[FWHM] > 0:
                full_width_half_max = self._options[FWHM]
            elif _f:
                full_width_half_max = _f.pFWHM
            else:
                full_width_half_max = 2

            if "x_init" in source_list.col_names:
                source_list.rename_column("x_init", X_CENTROID)
            if "y_init" in source_list.col_names:
                source_list.rename_column("y_init", Y_CENTROID)
            if "x_det" in source_list.col_names:
                source_list.rename_column("x_det", X_CENTROID)
            if "y_det" in source_list.col_names:
                source_list.rename_column("y_det", Y_CENTROID)
            if "flux_det" in source_list.col_names:
                source_list.rename_column("flux_det", "flux")
            mask = ~(np.isnan(source_list[X_CENTROID])
                     | np.isnan(source_list[Y_CENTROID]))


            bgd = BackGroundEstimateRoutine(
                source_list[mask], box_size=int(self._options[BOX_SIZE]),
                full_width_half_max=full_width_half_max,
                sig_sky=self._options[SIGSKY],
                bgd_r=self._options[BGD_R],
                profile_scale=self._options[PROF_SCALE],
                profile_slope=self._options[PROF_SLOPE],
                verbose=self._options[VERBOSE])
            header = self.header
            header.update(self._wcs.to_header())
            self._background = fits.ImageHDU(
                data=bgd(
                    self.image.data.copy(),
                    output=self._options.get(BGD_CHECKFILE)),
                header=header)
            if not self._options.get(QUIETMODE):
                f_name = "%s/%s-bgd.fits"%(self._out_dir, self._b_name)
                self.log("--> %s\n" % f_name)
                self._background.writeto(f_name, overwrite=True)

        else:
            p_error("unable to estimate background, no source list loaded\n")
            status = 1
        return status



    def bgd_subtraction(self):
        """
        Internally subtract a background array from an image array
        :return: 0 for success, 1 otherwise
        :rtype int.
        """
        self.log("Subtracting Background\n")

        if self._background is None:
            p_error("No background array loaded (-b file-bgd.fits)\n")
            return EXIT_FAIL
        array = self.image.data - self._background.data
        self._residuals = array
        self._image[self._n_hdu].data = array
        header = self.header
        header.update(self._wcs.to_header())
        fits.ImageHDU(
            data=self._residuals, name="RES", header=header).writeto(
                "%s/%s-res.fits" % (self._out_dir, self._b_name),
            overwrite=True)
        return EXIT_SUCCESS

    # noinspection SpellCheckingInspection
    def photometry(self):
        """
        Full photometry routine
        Saves the result as a table self._psf_catalogue,
        Additionally it appends a residual Image onto the
        self._residuals HDUList

        :return: 0 for success, 1 otherwise
        :rtype int
        """
        if self.image:
            self.log("\nRunning PSF Photometry\n")

            image, error, bgd, mask = self.prepare_image_arrays()

            if bgd is None:
                _, median, _ = (
                    sigma_clipped_stats(image, sigma=self._options[SIGSKY]))
                bgd = np.ones(self.image.shape) * median
                self.log(
                    "-> no background file loaded, measuring sigma "
                    "clipped median\n")

            ###################################
            # Collect relevant files and data #
            ###################################

            if self._detections is None:
                p_error("unable to run photometry: no source list loaded\n")
                return EXIT_FAIL

            if (self._psf is None
                and self.load_psf(
                    os.path.expandvars(self._options[PSF_FILE]))):
                p_error("unable to run photometry: no PSF loaded\n")
                return EXIT_FAIL

            psf_mask = ~np.isfinite(self._psf)
            if psf_mask.sum():
                self._psf[psf_mask] = 0
                self.log("-> masking INF pixels in PSF_FILE\n")

            psf_model = FittableImageModel(self._psf)
            if self._options[PSF_SIZE] > 0:
                size = int(self._options[PSF_SIZE])
            else:
                size = psf_model.shape[0]
            if not size % 2:
                size -= 1
            self.log("-> psf size: %d\n" % size)

            #########################
            # Sort out Init guesses #
            #########################
            app_hot_r = self._options.get(APPHOT_R)
            if not app_hot_r or app_hot_r <= 0:
                app_hot_r = 3

            init_guesses = self._detections.copy()
            if X_CENTROID in init_guesses.col_names:
                init_guesses.rename_column(X_CENTROID, "x_init")
            if Y_CENTROID in init_guesses.col_names:
                init_guesses.rename_column(Y_CENTROID, "y_init")
            if "x_det" in init_guesses.col_names:
                init_guesses.rename_column("x_det", "x_init")
            if "y_det" in init_guesses.col_names:
                init_guesses.rename_column("y_det", "y_init")

            init_guesses = init_guesses[ init_guesses["x_init"] >=0 ]
            init_guesses = init_guesses[ init_guesses["y_init"] >=0 ]
            init_guesses = init_guesses[
                init_guesses["x_init"] < self.image.header[NAXIS1]]
            init_guesses=init_guesses[
                init_guesses["y_init"] < self.image.header[NAXIS2]]

            ######
            # Allow tables that don't have the correct columns through
            ######
            required = ["x_init","y_init","flux",self._filter, "flag"]
            for notfound in  set(required) - set(init_guesses.col_names):
                dtype = np.uint16 if notfound == "flag" else float
                init_guesses.add_column(
                    Column(np.zeros(len(init_guesses)),
                           name=notfound, dtype=dtype))

            init_guesses = init_guesses[required]
            init_guesses.remove_column("flux")
            init_guesses.rename_column(self._filter, "ap_%s" % self._filter)
            
            ###########
            # Run Fit #
            ###########

            min_separation = self._options.get(CRIT_SEP)
            if not min_separation:
                min_separation = min(5, 2.5 * self._options.get(FWHM))

            if self._options[FORCE_POS]:
                phot = PSFPhotRoutine(
                    psf_model, size, min_separation=min_separation,
                    app_hot_r=app_hot_r, background=bgd, force_fit=1,
                    verbose=self._options[VERBOSE])
                psf_cat = phot(
                    image,init_params=init_guesses, error=error, mask=mask)
                psf_cat["flag"] |= SRC_FIX

            else:
                phot = PSFPhotRoutine(
                    psf_model, size, min_separation=min_separation,
                    app_hot_r=app_hot_r, background=bgd, force_fit=0,
                    verbose=self._options[VERBOSE])
                psf_cat = phot(
                    image,init_params=init_guesses, error=error, mask=mask)

                if not psf_cat:
                    return EXIT_FAIL

                ##################################
                # Setting position max variation #
                ##################################
                max_y_dev, unit = parse_unit(self._options[MAX_XY_DEV])
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
                    phot = PSFPhotRoutine(
                        psf_model, size, min_separation=min_separation,
                        app_hot_r=app_hot_r, background=bgd, force_fit=1,
                        verbose=self._options[VERBOSE])
                    ii = np.where(psf_cat["xydev"] > max_y_dev)
                    fixed_centres = psf_cat[ii][
                        ["x_init", "y_init", "ap_%s" % self._filter, "flag"]]
                    if len(fixed_centres):
                        self.log("-> forcing positions for deviant sources\n")
                        fixed_cat = phot(
                            image, init_params=fixed_centres, error=error,
                            mask=mask)
                        fixed_cat["flag"] |= SRC_FIX
                        psf_cat.remove_rows(ii)
                        psf_cat = vstack((psf_cat, fixed_cat))
                    else:
                        self.log("-> no deviant sources\n")

            ra, dec = self._wcs.all_pix2world(
                psf_cat["x_fit"], psf_cat["y_fit"], 0)
            psf_cat.add_column(Column(ra, name=RA), index=2)
            psf_cat.add_column(Column(dec, name=DEC), index=3)

            mag, mag_err = flux2mag(psf_cat["flux"],psf_cat["eflux"])

            filter_string = self._filter if self._filter else "mag"
            psf_cat.add_column(
                mag + self._options.get(ZP_MAG), name=filter_string)
            psf_cat.add_column(mag_err, name="e%s" % filter_string)
            self._psf_catalogue = psf_cat
            self._psf_catalogue.meta = dict(self.header.items())
            self._psf_catalogue.meta[AP_FILE]=self._options[AP_FILE]
            self._psf_catalogue.meta[BGD_FILE]=self._options[BGD_FILE]

            reindex(self._psf_catalogue)
            if not self._options.get(QUIETMODE):
                file_name = "%s/%s-psf.fits" % (self._out_dir, self._b_name)
                self.log("--> %s\n" % file_name)
                fits.BinTableHDU(
                    data=self._psf_catalogue,
                    header=self.header).writeto(file_name, overwrite=True)

            ##################
            # Residual Image #
            ##################

            if self._options[GEN_RESIDUAL]:
                self.log("-> generating residual\n")
                _tmp = psf_cat["x_fit", "y_fit", "flux"].copy()
                _tmp.rename_columns( ("x_fit", "y_fit"), ("x_0","y_0"))
                stars = make_model_image(
                    image.shape, psf_model, _tmp, model_shape=(size,size))
                residual = image - (bgd + stars)
                self._residuals = (
                    residual / get_mj_ysr2jy_scale_factor(self.image))
                header = self.header
                header.update(self._wcs.to_header())
                fits.ImageHDU(
                    data=self._residuals, name="RES", header=header).writeto(
                        "%s/%s-res.fits" % (self._out_dir, self._b_name),
                        overwrite=True)
            return EXIT_SUCCESS

    def source_geometry(self):
        """
        Calculate source geometry stats for a given image and source list
        :return: None
        """
        if self._detections is None:
            p_error("No source file loaded\n")
        else:
            self.log("Running Source Geometry\n")
            slist = self._filter_detections()

            sp = SourceProperties(
                self.image.data, slist, verbose=self._options[VERBOSE])
            stat = sp(
                fwhm=filters[self._filter].pFWHM,
                do_crowd=self._options[CALC_CROWD])
            
            self._source_stats = hstack((slist, stat))
            f_name = "%s/%s-stat.fits"%(self._out_dir, self._b_name)
            self.log("--> %s\n" % f_name)
            reindex(self._source_stats)
            fits.BinTableHDU(
                data=self._source_stats, header=self.header).writeto(
                    f_name, overwrite=True)

    # noinspection SpellCheckingInspection
    def verify(self):
        """
        This simple function verifies that everything necessary has been
        loaded properly

        :return: int where 0 on success, 1 on fail
        :rtype int
        """

        status = 0
        
        self.log("Checking internal systems..\n")

        if not self._filter:
            warn("No FILTER set, please set in parameter file or "
                 "use \"-s FILTER=XXX\"\n")
            status = 1

        d_name = os.path.expandvars(StarbugBase.get_data_path())
        if not os.path.exists(d_name):
            warn("Unable to locate STARBUG_DATDIR='%s'\n" % d_name)

        if not os.path.exists(self._out_dir):
            warn("Unable to locate OUTPUT='%s'\n" % self._out_dir)
            status = 1

        tmp=load_default_params()
        if set(tmp.keys()) - set(self._options.keys()):
            warn("Parameter file version mismatch. "
                 "Run starbug2 --update-param to update\n")
            status = 1
        
        if self._image is None or self.image.data is None:
            warn("Image did not load correctly\n")
            status = 1

        if self._options[AP_FILE] and self._detections is not None:
            test = self._filter_detections()
            if not len(test):
                warn("Detection file empty or no sources overlap the image.\n")
                status = 1

        return status

    def _filter_detections(self):
        """
        filters the detections based on some fixed constraints.
        :return: the filtered detections
        """
        detections = self._detections[[X_CENTROID,Y_CENTROID]].copy()
        detections = detections[ detections[X_CENTROID]>=0 ]
        detections = detections[ detections[Y_CENTROID]>=0 ]
        detections = detections[
            detections[X_CENTROID] < self.image.header[NAXIS1]]
        return detections[ detections[Y_CENTROID] < self.image.header[NAXIS2]]

    def __getstate__(self):
        """
        extracts the inner state of this class. deleting image or/and
         background if it's there.
        :return: the internal state with those bits filtered away
        """
        self._image.close()
        state = self.__dict__.copy()
        if "_image" in state:
            ##Sorry but we cant have that
            del state["_image"]
            ## This currently doesnt get reloaded
        if "_background" in state:
            del state["_background"]
            
        return state

    def __setstate__(self, state):
        self.__dict__.update(state)
        v = self._options[VERBOSE]
        self._options[VERBOSE] = 0
        self.load_image(self._f_name)
        self._options[VERBOSE] = v
