import os
import sys
import numpy as np

from astropy.stats import sigma_clip
from astropy.table import Column, Table

from photutils.aperture import (
    CircularAperture, CircularAnnulus, aperture_photometry)

from starbug2.constants import (
    SRC_GOOD, DQ_DO_NOT_USE, DQ_SATURATED, SRC_BAD, DQ_JUMP_DET, SRC_JMP,
    X_CENTROID, Y_CENTROID, E_FLUX, FLUX, FILTER_LOWER, EXIT_FAIL, EE_FRACTION,
    RADIUS, AP_CORR)
from starbug2.utils import printf, p_error, warn


class APPhotRoutine:

    @staticmethod
    def calc_ap_corr(filter_string, radius, table_f_name=None, verbose=0):
        """
        Using CRDS ap_corr table, fit a curve to the radius vs ap_corr
        columns and then return ap_corr to respective input radius

        :param filter_string: the filter string
        :type filter_string: str
        :param radius: the radius
        :type radius: float
        :param table_f_name: the table file name
        :type table_f_name: str
        :param verbose: int for verbose.
        :type verbose: int
        :return: the ap_corr table.
        :rtype: np.array
        """

        if not table_f_name or not os.path.exists(table_f_name):
            return EXIT_FAIL
        tmp = Table.read(table_f_name, format="fits")

        if FILTER_LOWER in tmp.col_names:
            t_ap_corr = tmp[(tmp[FILTER_LOWER] == filter_string)]
        else:
            t_ap_corr = tmp


        if "pupil" in t_ap_corr.col_names:
            t_ap_corr = t_ap_corr[ t_ap_corr["pupil"] == "CLEAR"]

        ap_corr = np.interp(radius, t_ap_corr[RADIUS], t_ap_corr[AP_CORR])
        if verbose:
            printf("-> estimating aperture correction: %.3g\n" % ap_corr)
        return ap_corr


    @staticmethod
    def ap_corr_from_enc_energy(
        filter_string, encircled_energy, table_f_name=None, verbose=0):
        """
        Rather than fitting radius to the AP_CORR CRDS, use the closes
        Encircled energy value

        :param filter_string: the filter string
        :type filter_string: str
        :param encircled_energy: the encircled energy
        :type encircled_energy: np.array
        :param table_f_name: the table file name
        :type table_f_name: str
        :param verbose: int for verbose.
        :type verbose: int
        :return: tuple of ap_corr and radius or ExitFail
        :rtype: tuple of np.array and np.array or int.
        """
        if not table_f_name or not os.path.exists(table_f_name):
            return EXIT_FAIL

        tmp = Table.read(table_f_name, format="fits")

        if FILTER_LOWER in tmp.col_names:
            t_ap_corr = tmp[(tmp[FILTER_LOWER] == filter_string)]
        else: t_ap_corr = tmp

        line = t_ap_corr[
            (np.abs(t_ap_corr[EE_FRACTION] - encircled_energy)).argmin()]
        if verbose:
            printf(
                "-> best matching encircled energy %.1f,"
                " with radius %g pixels\n" % (
                    line[EE_FRACTION], line[RADIUS]))
            printf("-> using aperture correction: %f\n" % line[AP_CORR])

        return line[AP_CORR], line[RADIUS]


    @staticmethod
    def radius_from_enc_energy(filter_string, ee_frac, table_f_name):
        """
        """
        if not table_f_name or not os.path.exists(table_f_name):
            return -1
        t_ap_corr = Table.read(table_f_name, format="fits")

        if len({EE_FRACTION, RADIUS} & set(t_ap_corr.col_names)) != 2:
            return -1

        # Crop down table
        if FILTER_LOWER in t_ap_corr.col_names:
            t_ap_corr=t_ap_corr[(t_ap_corr[FILTER_LOWER] == filter_string)]

        if "pupil" in t_ap_corr.col_names: # Crop down table
            t_ap_corr=t_ap_corr[ t_ap_corr["pupil"] == "CLEAR"]

        return np.interp(ee_frac, t_ap_corr[EE_FRACTION], t_ap_corr[RADIUS])

    def __init__(self, radius, sky_in, sky_out, verbose=0):
        """
        Aperture photometry called by starbug

        :param radius: Pixel radius of photometric aperture
        :type radius: float
        :param sky_in: Pixel radius of inner sky annuli
        :type sky_in: float
        :param sky_out: Pixel radius of outer sky annuli
        :type sky_out: float
        :param verbose: Set whether to print verbose output information
        :type verbose: bool
        """
        if sky_in < radius:
            warn("Sky annulus radii must be larger than aperture radius.\n")
            sky_in = radius + 1

        if sky_in >= sky_out:
            warn("Sky annulus outer radii must be larger than the inner.\n")
            sky_out = sky_in + 1

        self.radius = radius
        self.sky_in = sky_in
        self.sky_out = sky_out
        self.catalogue = Table(None)
        self.verbose = verbose

    def __call__(self, image, detections, **kwargs):
        return self.run(image, detections, **kwargs)

    def run(self, image, detections, error=None, dq_flags=None, ap_corr=1.0,
            sig_sky=3):
        """
        Forced aperture photometry on a list of detections are an
         `astropy.table.Table` with columns x_centroid y_centroid or x_0 y_0
        This will add extra columns into this table ap_flux ap_sky

        :param detections: Table with source positions containing
                           x_centroid, y_centroid columns
        :type detections: astropy.table.Table
        :param image: 2D Image data to run photometry on
        :type image: numpy.ndarray
        :param error: 2D Image array containing photometric error per pixel
        :type error: numpy.ndarray
        :param dq_flags: 2D Image array containing JWST data quality flags per
                         pixel
        :type dq_flags: numpy.ndarray
        :param ap_corr: Aperture correction to be applied to the flux
        :type ap_corr: float
        :param sig_sky: Sigma threshold above the median to clip from sky
                        apertures
        :type sig_sky: float
        :return: Photometry catalogue
        :rtype: astropy.table.Table
        """
        if len({X_CENTROID, Y_CENTROID} & set(detections.colnames)) == 2:
            pos = [(line[X_CENTROID],line[Y_CENTROID]) for line in detections]
        elif len({"x_0", "y_0"} & set(detections.colnames))==2:
            pos = [(line["x_0"], line["y_0"]) for line in detections]
        elif len({"x_init", "y_init"} & set(detections.colnames)) == 2:
            pos=[(line["x_init"], line["y_init"]) for line in detections]
        else:
            p_error(
                "Cannot identify position in detection catalogue ("
                "x_0/x_centroid)\n")
            return None

        mask = np.isnan(image)
        if error is None:
            error = np.sqrt(image)

        apertures = CircularAperture(pos, self.radius)
        smooth_apertures = CircularAperture(
            pos, min(1.5 * self.radius, self.sky_in))
        annulus_aperture = CircularAnnulus(
            pos, r_in=self.sky_in, r_out=self.sky_out)

        self.log(
            "-> apertures: %.2g (%.2g - %.2g)\n" % (
                self.radius, self.sky_in, self.sky_out))
        phot = aperture_photometry(
            image, (apertures, smooth_apertures), error=error, mask=mask)
        self.catalogue=(
            Table(np.full((len(pos), 4), np.nan),
                  names=("smoothness", "flux", "eflux", "sky")))

        self.log("-> calculating sky values\n")
        masks = annulus_aperture.to_mask(method="center")
        dat = list(map(lambda a : a.multiply(image), masks))

        try:
            dat = np.array(dat).astype(float)
        except (ValueError, TypeError) as e:
            ## Cases where the array is inhomogeneous
            ## If annulus reaches the edge of the image, it will create a
            # mask the wrong shape. If for whatever reason the point lies
            # outside the image, it will have None in the list, this needs
            # to be caught too
            warn(
                f"Ran into issues with the sky annuli. {e},"
                f" trying to fix them..\n")
            size = np.max([np.shape(d) for d in dat if d is not None])
            for i, d in enumerate(dat):
                if d is None:
                    dat[i] = np.zeros((size,size))
                elif (shape := np.shape(d)) != (size, size):
                    dat[i] = np.zeros((size,size))
                    dat[i][:shape[0],:shape[1]] += d
            dat = np.array(dat)

        mask = (dat > 0 & np.isfinite(dat))
        dat[~mask] = np.nan
        dat = sigma_clip(dat.reshape(dat.shape[0],-1), sigma=sig_sky, axis=1)
        self.catalogue["sky"] = np.ma.median(dat, axis=1).filled(fill_value=0)
        std = np.ma.std(dat, axis=1)

        e_poisson = phot["aperture_sum_err_0"]
        esky_scatter = apertures.area * std ** 2
        esky_mean = (std ** 2 * apertures.area ** 2) / annulus_aperture.area

        self.catalogue[E_FLUX] = np.sqrt(
            e_poisson ** 2 + esky_scatter ** 2 + esky_mean ** 2)
        self.catalogue[FLUX] = (
            ap_corr * (phot["aperture_sum_0"] - (
            self.catalogue["sky"] * apertures.area)))

        self.catalogue[FLUX][self.catalogue[FLUX] == 0] = np.nan

        ######################
        # Source "smoothness", the gradient of median pixel values within the
        # two test apertures
        ######################

        self.catalogue["smoothness"] = (
            (phot["aperture_sum_1"] / smooth_apertures.area)
            / (phot["aperture_sum_0"] / apertures.area))

        col = Column(
            np.full(len(apertures), SRC_GOOD), dtype=np.uint16, name="flag")
        if dq_flags is not None:
            self.log("-> flagging unlikely sources\n")
            for i, mask in enumerate(apertures.to_mask(method="center")):
                tmp = mask.multiply(dq_flags)
                if tmp is not None:
                    dat = np.array(tmp,dtype=np.uint32)
                    if np.sum( dat & (DQ_DO_NOT_USE | DQ_SATURATED)):
                        col[i] |= SRC_BAD
                    if np.sum( dat & DQ_JUMP_DET):
                        col[i] |= SRC_JMP
        self.catalogue.add_column(col)
        return self.catalogue


    def log(self, msg):
        """
        log message if in verbose mode

        :param msg: message to log
        :return: None
        """
        if self.verbose:
            printf(msg)
            sys.stdout.flush()