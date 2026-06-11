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

import numpy as np
from typing import cast, Any
from photutils.datasets import make_model_image, make_random_models_table
from photutils.psf import ImagePSF
from astropy.table import Table, hstack
from astropy.io import fits
from scipy.optimize import curve_fit

from starbug2.constants import (
    X_0, Y_0, MAG, FLUX, X_DET, Y_DET, FLUX_DET, STATUS, ZP_MAG, ID,
    X_CENTROID, Y_CENTROID, REC, EXIT_SUCCESS)
from starbug2.matching.generic_match import GenericMatch

try:
    import matplotlib.pyplot as plt
except ImportError:
    import matplotlib; matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

from starbug2.utils import (
    printf, p_error, get_mj_ysr2jy_scale_factor, warn)


class ArtificialStarsIII:
    # not found
    NULL = 0

    # found
    DETECT = 1


    """
    ast
    """
    def __init__(self, starbug, index=-1):
        ## Initials the starbug instance
        self.starbug = starbug
        _ = self.starbug.main_image
        psf_success = self.starbug.load_psf()

        if psf_success != EXIT_SUCCESS:
            warn("the psf file was not loaded. Expected failure.")
            raise Exception("the psf file failed to load.")

        self.psf = ImagePSF(self.starbug.psf)
        self.index = index

    def __call__(self, *args, **kwargs):
        return self.auto_run(*args, **kwargs)

    def auto_run(
            self, n_tests, stars_per_test=1, sub_image_size=-1,
            mag_range=(18, 27), loading_buffer=None, autosave=-1, skip_phot=0,
            skip_background=0):
        """
        The main entry point into the artificial star test.
        This handles everything except the results compilation at the end.

        :param n_tests: Number of tests to run.
        :type n_tests: int
        :param stars_per_test: Number of stars to inject per test.
        :type stars_per_test: int
        :param sub_image_size: In prep.
        :type sub_image_size: int
        :param mag_range: Length two list or tuple containing the magnitude
                          range of injected stars. These will be uniformly
                          sampled from within this range.
        :type mag_range: tuple[float, float] or list[float]
        :param loading_buffer: Length 3 array of shared memory to increment a
                              loading bar between multiple subprocesses.
        :type loading_buffer: numpy.ndarray
        :param autosave: Auto quick saving output frequency.
        :type autosave: int
        :param skip_phot: If true then ignore the PSF phot step and use
                          aperture fluxes instead.
        :type skip_phot: bool or int
        :param skip_background: If true then ignore the background
                                subtraction step.
        :type skip_background: bool or int
        :return: Full raw test results. Injected initial properties with
                 measured values.
        :rtype: astropy.table.Table
        """

        test_result = Table(
            np.full((n_tests * stars_per_test, 8), np.nan),
            names=[X_0, Y_0, MAG, FLUX, X_DET, Y_DET, FLUX_DET, STATUS])
        scale_factor = get_mj_ysr2jy_scale_factor(self.starbug.main_image)
        base_image = self.starbug.image.copy()
        base_shape = np.copy(self.starbug.main_image.shape)
        stars_per_test = int(stars_per_test)
        passed = 0

        z_p = (
            self.starbug.options.get(ZP_MAG) if
            self.starbug.options.get(ZP_MAG) else 0)
        buffer = 0

        if mag_range[0] - mag_range[1] >= 0:
            warn("Detected magnitude range in wrong order,"
                 " put bright limit first\n")
            return None

        if any(base_shape < sub_image_size):
            sub_image_size = min(base_shape)
            p_error("sub image size greater than image size, setting to "
                    "'safe' value %d.\n" % sub_image_size)
            
        for test in range(1, int(n_tests) + 1):
            image = base_image.__deepcopy__()

            shape = image[self.starbug.n_hdu].shape

            source_list = make_random_models_table(
                stars_per_test,
                { X_0 : [buffer, shape[0] - buffer],
                  Y_0 : [buffer, shape[1] - buffer],
                  MAG :  mag_range
                })
            source_list.add_column(
                10.0 ** ( (z_p - source_list[MAG]) / 2.5 ) , name=FLUX)
            source_list.remove_column(ID)

            star_overlay = (
                make_model_image(
                    shape, self.psf, source_list,
                    model_shape=self.psf.data.shape)
                / scale_factor)
            image[self.starbug.n_hdu].data += star_overlay
            self.starbug.image = image

            result = self.single_test(
                source_list, skip_phot=skip_phot,
                skip_background=skip_background)

            passed += sum(result[STATUS])
            test_result[
                (test - 1) * stars_per_test:
                test * stars_per_test] = result

            if loading_buffer is not None:
                loading_buffer[0] += 1
                loading_buffer[2] = int(100 * passed / (test * stars_per_test))

            if autosave > 0 and not test % autosave:
                # noinspection SpellCheckingInspection
                test_result.write(
                    "sbast-autosave%d.tmp" % self.index, overwrite=True,
                    format="fits")

            # is this necessary?
            del image
        return test_result

    def single_test(self, contains, skip_phot=0, skip_background=0):
        """
        Conduct a single test on an image with a set of initial source
        properties.

        :param contains: Table of initial source properties to be injected
                         into the image. This table must contain the columns
                         ("x_0", "y_0", "flux").
        :type contains: astropy.table.Table
        :param skip_phot: Skip the PSF phot routine.
        :type skip_phot: bool or int
        :param skip_background: Skip the background estimation and
                                subtraction step.
        :type skip_background: bool or int
        :return: Table horizontally stacked with the initial inputs and the
                 detection and  photometric results. Plus column named
                 "status", an integer flag whether the source was detected or
                 not.
        :rtype: astropy.table.Table
        """
        test_result = Table(
            np.full((len(contains), 4), np.nan),
            names=[X_DET, Y_DET, FLUX_DET, STATUS])

        threshold = 2

        #Run detection on the image
        if not self.starbug.detect():
            det = self.starbug.detections

            #Check for detection in output
            for i, src in enumerate(contains): # type: ignore
                separations = np.sqrt(
                    (src[X_0] - det[X_CENTROID]) ** 2
                    + (src[Y_0] - det[Y_CENTROID]) ** 2)
                best_match = np.argmin(separations)
                if separations[best_match] < threshold:
                    test_result[X_DET][i] = det[X_CENTROID][best_match]
                    test_result[Y_DET][i] = det[Y_CENTROID][best_match]
                    test_result[FLUX_DET][i] = det[FLUX][best_match]
                    test_result[STATUS][i] = self.DETECT
                else:
                    test_result[STATUS][i] = self.NULL

            # Run background
            if (sum(test_result[STATUS])
                and (skip_background or not self.starbug.bgd_estimate())):
                # estimate if there were detections
                self.starbug.detections = test_result

                # Run PSF photometry on detected sources
                if not skip_phot and not self.starbug.photometry_routine():
                    # noinspection SpellCheckingInspection
                    self.starbug.psf_catalogue.rename_columns(
                        ("x_init", "y_init", "xydev"),
                        ("_x_init", "_y_init", "_xydev"))
                    matched = GenericMatch(threshold=threshold)(
                        [contains, self.starbug.psf_catalogue],
                        cartesian=True)
                    test_result[FLUX_DET] = (
                        matched[:len(test_result)]["flux_2"])
        return hstack((contains, test_result))


def get_completeness(test_result):
    """
    Compile the results into magnitude binned values of recovery fraction
    and flux error.

    :param test_result: The output from auto_run.
    :type test_result: astropy.table.Table
    :return: A table containing percent completeness as a function of
             magnitude.
    :rtype: astropy.table.Table
    """

    bins = np.arange(
        np.floor(min(test_result[MAG])),
        np.ceil(max(test_result[MAG])),
        0.1)
    percents = np.zeros(len(bins))
    errors = np.zeros(len(bins))
    offsets = np.zeros(len(bins))
    means = np.zeros(len(bins))
    
    i_bins = np.digitize( test_result[MAG], bins=bins)
    for i in range(max(i_bins)):
        binned = test_result[ (i_bins==i) ]
        if binned:
            percents[i] = float(sum(binned[STATUS])) / len(binned)

        mag_inj = -2.5 * np.log10( binned[FLUX])
        mag_det = -2.5 * np.log10( binned[FLUX_DET])
        errors[i] = np.nanstd( mag_inj - mag_det )
        means[i] = np.nanmean( mag_inj - mag_det )
        offsets[i] = np.nanmedian(binned[FLUX] / binned[FLUX_DET])


    out = Table(
        [bins, percents, errors, offsets],
        names=("mag", "rec", "err", "off"),
        dtype=(float, float, float, float))
    return out

def get_spatial_completeness(test_result, image, res=10):
    """
    Produce an image array showing the spatially dependent recovery fraction.

    :param test_result: The output from auto_run.
    :type test_result: astropy.table.Table
    :param image: 2D image array to take the shape from.
    :type image: numpy.ndarray
    :param res: The resolution of the spatial bins.
    :type res: int
    :return: A 2D array the same shape as the image input, where pixel values
        show the fraction of injected sources recovered in this bin.
    :rtype: numpy.ndarray
    """
    if not image: return None
    x_bins = np.arange(
        min(test_result[X_0]), max(test_result[X_0]), int(res))
    y_bins = np.arange(min(test_result[Y_0]),max(test_result[Y_0]), int(res))
    percents = np.zeros(image.shape)
    for xi in x_bins[:-1]:
        xo = xi + res
        for yi in y_bins[:-1]:
            yo = yi + res
            mask = (
                (test_result[X_0] >= xi) &
                (test_result[X_0] < xo) &
                (test_result[Y_0] >= yi) &
                (test_result[Y_0] < yo))
            binned = test_result[mask]
            if len(binned):
                percents[int(xi): int(xo), int(yi): int(yo)] = (
                    float(sum(binned[STATUS]) / len(binned)))
    return percents

def estimate_completeness_mag(ast):
    """
    Estimate the completeness level of the artificial star test.

    :param ast: Output of Artificial_Stars.get_completeness, table must
               contain columns (mag, rec).
    :type ast: astropy.table.Table
    :return: A tuple containing:
        - **fit** (*list*): The fitting parameters to the logistic curve
                            $f(x) = \frac{l}{1 + \exp(-k(x - x_0))}$ formatted
                            as ``[l, x_0, k]``.
        - **complete** (*list*): Magnitude of 70% and 50% completeness.
    :rtype: tuple[list, list]
    """
    fit = [None, None, None]
    completeness = [None, None, None]
    fn_i = lambda y, l, k, xo: xo - (np.log((l / y) - 1) / k)

    if len(set(ast.colnames) & {MAG, REC}) == 2:
        try:
            # need the *_ as the return tuple can be multiple sizes. The *_
            # allows the IDE to not freak out, especially as we don't care
            # about the rest of the return values.
            fit, *_ = curve_fit(
                scurve, ast[MAG], ast[REC], [1, -1, np.median(ast[MAG])])
            completeness = (fn_i(0.9, *fit), fn_i(0.7, *fit), fn_i(0.5, *fit))
        except (RuntimeError, ValueError) as e:
            warn(f"Unable to fit completeness fractions: {e}\n")
    else:
        p_error("Input table must have columns 'mag' and 'rec'\n")
    return fit, completeness

def scurve(x, l, k, xo):
    """
    S-curve function to fit completeness results to.

    math::  f(x) = \\frac{l}{1 + \\exp(-k(x - x_0))}

    :param x: Magnitude range or array to input into the function.
    :type x: list or numpy.ndarray
    :param l: Maximum value asymptote (typically representing maximum
              completeness, near 1.0).
    :type l: float
    :param xo: The inflection point of the curve (the magnitude where
               completeness is 50%).
    :type xo: float
    :param k: The logistic growth rate or steepness of the curve.
    :type k: float
    :return: Calculated function value(s) matching the shape of the
             input ``x``.
    :rtype: float or numpy.ndarray
    """
    return l / (1 + np.exp(-k * (x - xo)))

def compile_results(raw, image=None, plot_ast=None, filter_string="m"):
    """
    Compile all the raw data into usable results

    :param raw:
    :param image:
    :param plot_ast:
    :param filter_string:
    :return:
    """

    completeness = get_completeness(raw)
    _cfit, _completeness = estimate_completeness_mag(completeness)
    spatial_completeness = get_spatial_completeness(raw, image, res=10)

    head = {
        "COMPLETE_FN":"F(x)=l/(1+exp(-k(x-xo)))", "l":_cfit[0], "k":_cfit[1],
        "xo":_cfit[2] }
    for i, frac in enumerate((90, 70, 50)):
        if _completeness[i] and not np.isnan(_completeness[i]):
            printf(
                "-> complete to %d%%: %s=%.2f\n" % (
                    frac, filter_string, _completeness[i]))
            head["COMPLETE %d%%" % frac] = _completeness[i]

    # needed for spatial_completeness as it expects an 'array.pyi',
    # got 'ndarray' instead.
    results= fits.HDUList(
        [fits.PrimaryHDU(header=fits.Header(head)),
         fits.BinTableHDU(data=completeness, name="AST"),
         fits.BinTableHDU(data=raw, name="RAW"),
         fits.ImageHDU(data=cast(Any, spatial_completeness), name="CMP")])

    if plot_ast:
        fig, ax = plt.subplots(1, figsize=(3.5,3), dpi=300)
        ax.scatter(completeness[MAG], completeness[REC], c='k', lw=0, s=8)
        ax.plot(completeness[MAG], scurve(completeness[MAG],*_cfit),
                c='g',
                label=r"$f(x)=\frac{%.2f}{1+e^{%.2f("r"x-%.2f)}}$" % (
                    _cfit[0],-_cfit[1],_cfit[2]))
        ax.axvline(
            _completeness[0], c="seagreen", ls='--',
            label=("90%%:%.2f" % _completeness[0]), lw=0.75)
        ax.axvline(
            _completeness[1], c="seagreen", ls='-.',
            label=("70%%:%.2f" % _completeness[1]), lw=0.75)
        ax.axvline(
            _completeness[2], c="seagreen", ls=':',
            label=("50%%:%.2f" % _completeness[2]), lw=0.75)
        ax.scatter(_completeness, (0.9, 0.7, 0.5), marker='*', c='teal', s=10)
        ax.tick_params(direction="in",top=True, right=True)
        ax.set_title("Artificial Star Test")
        ax.set_xlabel(filter_string)
        ax.set_ylabel("Fraction Recovered")
        ax.set_yticks([0,.25,.5,.75,1])
        ax.legend(loc="lower left", frameon=False, fontsize=8)
        plt.tight_layout()
        fig.savefig(plot_ast, dpi=300)
        printf("--> %s\n" % plot_ast)
    return results
