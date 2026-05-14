"""
Core routines for StarbugII.
"""
import sys
from typing import overload

import numpy as np
from astropy.table import Column, hstack
from photutils.aperture import (
    CircularAperture, aperture_photometry)
from photutils.psf import PSFPhotometry, SourceGrouper
from starbug2.utils import printf, p_error, warn

class _Grouper(SourceGrouper):
    """
    Overloaded SourceGrouper. This class gives a starbug warning into stderr
     if there are more than CRITICAL_VAL source in any given group. Then
     returns the original __call__ function results.

    Parameters:
    -----------
    min_separation : float > 0
        The minimum distance (in pixels) such that any two sources
        separated by less than this distance will be placed in the same
        group if the ``min_size`` criteria is also met.
    """
    CRITICAL_VAL = 25
    min_separation = 0

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __call__(self, *args, **kwargs):
        res = super().__call__(*args, **kwargs)
        if n := sum(np.bincount(res) > self.CRITICAL_VAL):
            warn("Source grouper has %d groups larger than %d. Consider"
                 " reducing \"CRIT_SEP=%g\" or fitting might take a long"
                 " time.\n" % (n,self.CRITICAL_VAL, self.min_separation))
        return res

class PSFPhotRoutine(PSFPhotometry):
    def __init__(self, psf_model, fit_shape, app_hot_r=3, min_separation=8,
                 force_fit=False, background=None, verbose=1):
        """
        PSF Photometry routine called by starbug

        :param psf_model: Model PSF to be used in the fitting
        :type psf_model: photutils.psf.FittableImageModel
        :param fit_shape: Size of PSF to use in pixels. Must be less than or
                          equal to the size of psf_model
        :type fit_shape: int or tuple
        :param min_separation: Minimum source separation for source grouper
                              (pixels)
        :type min_separation: float
        :param app_hot_r: Aperture radius to be used in initial guess
                          photometry
        :type app_hot_r: float
        :param force_fit: Conduct forced centroid PSF fitting
        :type force_fit: bool
        :param background: 2D array with the same dimensions as the data used
                           in fitting
        :type background: numpy.ndarray
        :param verbose: Show verbose outputs
        :type verbose: bool or int
        """
        self._verbose = verbose
        self._force_fit = force_fit
        self._background = background

        grouper = _Grouper(min_separation)

        if force_fit:
            psf_model.x_0.fixed = True
            psf_model.y_0.fixed = True

        super().__init__(
            psf_model=psf_model, fit_shape=fit_shape, finder=None,
            progress_bar=verbose, aperture_radius=app_hot_r, grouper=grouper)

        if self._verbose:
            printf("-> source group separation: %g\n" % min_separation)

    def __call__(self, *args, **kwargs):
        """
        runs the psf phot routine.
        :param args: args
        :type args: dict
        :param kwargs: extra args
        :type kwargs: dict
        :return: processed table
        :rtype: astropy.table.Table
        """
        return self.do_photometry(*args, **kwargs)


    def do_photometry(
            self, image, init_params=None, error=None, mask=None,
            progress_bar=False):
        """
        does the photometry
        :param image: the image to process.
        :type image: astropy.table.Table
        :param init_params: the init params.
        :type init_params: dict of str, str
        :param error: the error.
        :type error: ????
        :param mask: the mask.
        :type mask: np.array
        :param progress_bar: the progress bar.
        :type progress_bar: ??????
        :return: the processed table.
        :rtype:  astropy.table.Table
        """

        if init_params is None or len(init_params) == 0:
            p_error("Must include source list\n")
            return None

        ### Removing completely masked sources
        apertures = CircularAperture(
            [(l["x_init"],l["y_init"]) for l in init_params],
            self.aperture_radius)
        ap_masks = aperture_photometry(~mask, apertures)
        init_params.remove_rows(ap_masks["aperture_sum"] == 0)

        ## bad errors should be big not small
        error[error == 0] = sys.maxsize

        if self._background is not None:
            image = image - self._background
        if self._verbose:
            printf("-> fitting %d sources\n"%len(init_params))
        cat = super().__call__(
            image, mask=mask, init_params=init_params, error=error)

        d = np.sqrt((
            (cat["x_init"] - cat["x_fit"]) ** 2.0 +
            (cat["y_init"]-cat["y_fit"]) ** 2.0))
        cat.add_column(Column(d, name="xydev"))

        if "flux_err" not in cat.col_names:
            cat.add_column(Column(np.full(len(cat), np.nan), name="eflux"))
            warn("Something went wrong with PSF error fitting\n")
        else:
            cat.rename_column("flux_err","eflux")

        cat.rename_column("flux_fit", "flux")

        keep=["x_fit", "y_fit", "flux", "eflux", "xydev", "qfit"]
        return hstack((init_params, cat[keep]))