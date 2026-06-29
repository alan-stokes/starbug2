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
import numpy as np
from typing import cast, Any, Final, List, Tuple, Callable, Dict

from astropy.io.fits import HDUList
from photutils.datasets import make_model_image, make_random_models_table
from photutils.psf import ImagePSF
from astropy.table import Table, hstack, QTable
from astropy.io import fits
from astropy import units
from astropy.units import Quantity
from scipy.optimize import curve_fit
from matplotlib.figure import Figure
from matplotlib.axes import Axes

from starbug2.core.constants import ExitStates, TableColumn
from starbug2.matching.generic_match import GenericMatch
from starbug2.interfaces.star_bug_interface import StarBugInterface


try:
    import matplotlib.pyplot as plt
except ImportError:
    import matplotlib
    matplotlib.use("TkAgg")
    import matplotlib.pyplot as plt

from starbug2.utilities.utils import (
    printf, p_error, get_mj_ysr2jy_scale_factor, warn)


class ArtificialStars:
    # not found
    NULL: Final[int] = 0

    # found
    DETECT: Final[int] = 1

    # the number of columns in the test table.
    N_COLUMNS: Final[int] = 8

    # the column names of the table
    TEST_TABLE_COLUMN_NAMES: Final[List[str]] = [
        TableColumn.X_0, TableColumn.Y_0, TableColumn.MAG, TableColumn.FLUX,
        TableColumn.X_DET, TableColumn.Y_DET, TableColumn.FLUX_DET,
        TableColumn.STATUS]

    # mag ranges
    MAG_RANGE_LOW: Final[int] = 18
    MAG_RANGE_HIGH: Final[int] = 27

    """
    ast
    """
    def __init__(self,
                 starbug: StarBugInterface,
                 index: int = -1) -> None:
        # Initials the starbug instance
        self._starbug: StarBugInterface = starbug
        _ = self._starbug.main_image
        psf_success: int = self._starbug.load_psf()

        if psf_success != ExitStates.EXIT_SUCCESS:
            warn("the psf file was not loaded. Expected failure.")
            raise Exception("the psf file failed to load.")

        self._psf: ImagePSF = ImagePSF(self._starbug.psf)
        self._index: int = index

    def execute_ast(
            self, n_tests: int, stars_per_test, sub_image_size: int,
            mag_range: Tuple[int, int], loading_buffer: np.ndarray | None,
            autosave: int, skip_phot: bool | int, skip_background: bool | int,
            zp_mag: float, save_image: bool,
            save_image_path: str, ast_seed: int | None) -> Table | None:
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
        :param zp_mag: the zero point magnitude
        :type zp_mag: float
        :param save_image: bool saying if we should save the added star image.
        :type save_image: bool
        :param save_image_path: the path for the save image.
        :type save_image_path: str
        :param ast_seed: the seed to set the artificial stars with. or None.
        :type ast_seed: int | None
        :return: Full raw test results. Injected initial properties with
                 measured values.
        :rtype: astropy.table.Table
        """
        return self._auto_run(
            n_tests, stars_per_test, sub_image_size, mag_range, loading_buffer,
            autosave, skip_phot, skip_background, zp_mag, save_image,
            save_image_path, ast_seed)

    def _add_stars(
            self, base_image: fits.HDUList, stars_per_test: int, buffer: int,
            mag_range: Tuple[int, int], zp_mag: float,
            scale_factor: float | int, save_image: bool, test: int,
            save_image_path: str, ast_seed: int | None) -> QTable:
        """
        adds new stars to the image.
        :param base_image: copy of the current image
        :type base_image: fits.HDUList
        :param stars_per_test: the number of stars to put in per test.
        :type stars_per_test: int
        :param buffer: the buffer
        :type buffer: int
        :param mag_range: the magnitude ranges
        :type mag_range: Tuple[int, int]
        :param zp_mag: the zero point magnitude.
        :type zp_mag: float
        :param scale_factor: the scale factor.
        :type scale_factor: float | int
        :param save_image: bool saying if we should save the added star image.
        :type save_image: bool
        :param test: the test id
        :type test: int
        :param save_image_path: the path to save the image to
        :type save_image_path: str
        :param ast_seed: the seed to set the artificial stars with. or None.
        :type ast_seed: int | None
        :return: the location of the new stars.
        :rtype: astrophy.QTable
        """

        image: fits.HDUList = base_image.__deepcopy__()

        shape: list[int, int] = image[self._starbug.n_hdu].shape # noqa

        source_list: QTable = make_random_models_table(
            stars_per_test, {
                TableColumn.X_0: [buffer, shape[0] - buffer],
                TableColumn.Y_0: [buffer, shape[1] - buffer],
                TableColumn.MAG:  mag_range
            }, ast_seed
        )
        source_list.add_column(
            10.0 ** ((zp_mag - source_list[TableColumn.MAG]) / 2.5),
            name=TableColumn.FLUX)
        source_list.remove_column(TableColumn.ID)

        star_overlay: np.ndarray = (
            make_model_image(
                shape, self._psf, source_list,
                model_shape=self._psf.data.shape)
            / scale_factor)
        image[self._starbug.n_hdu].data += star_overlay
        self._starbug.image = image

        if save_image:
            image.writeto(os.path.join(
                save_image_path, f"inserted_image_for_test_{test}.fits"))

        return source_list

    def _execute_test(
            self, base_image: fits.HDUList, stars_per_test: int, buffer: int,
            mag_range:  Tuple[int, int], zp_mag: float,
            scale_factor: float | int, skip_phot: bool | int,
            skip_background: bool | int, passed: int, test_result: Table,
            test: int, loading_buffer: np.ndarray | None,
            autosave: int, save_image: bool, save_image_path: str,
            ast_seed: int | None) -> None:
        """
        executes a test.
        :param base_image: copy of the current image
        :type base_image: fits.HDUList
        :param stars_per_test: the number of stars to put in per test.
        :type stars_per_test: int
        :param buffer: the buffer
        :type buffer: int
        :param mag_range: the magnitude ranges
        :type mag_range: Tuple[int, int]
        :param zp_mag: the zero point magnitude.
        :type zp_mag: float
        :param scale_factor: the scale factor.
        :type scale_factor: float | int
        :param skip_phot: check for if we should skip phot
        :type skip_phot: bool
        :param skip_background: check for if we should skip background
                                reduction.
        :type skip_background: bool
        :param passed: int to accum if the test works
        :type passed: int
        :param test_result:  the result location for the tests output
        :type test_result: Table
        :param test: the test id
        :type test: int
        :param loading_buffer: the loading buffer
        :type loading_buffer: np.ndarray | None
        :param autosave: the auto save value.
        :type autosave: int
        :param save_image: bool saying if we should save the added star image.
        :type save_image: bool
        :param save_image_path: the path to save the new image to.
        :type save_image_path: str
        :param ast_seed: the seed to set the artificial stars with. or None.
        :type ast_seed: int | None
        :return: None
        """
        source_list = self._add_stars(
            base_image, stars_per_test, buffer, mag_range, zp_mag,
            scale_factor, save_image, test, save_image_path, ast_seed)

        result: Table = self.single_test(
            source_list, skip_phot=skip_phot,
            skip_background=skip_background)

        # save system memory.
        image: HDUList | None = self._starbug.image
        assert image is not None
        image.close()

        # process result.
        passed += sum(result[TableColumn.STATUS])
        test_result[
            (test - 1) * stars_per_test:
            test * stars_per_test] = result

        if loading_buffer is not None:
            loading_buffer[0] += 1
            loading_buffer[2] = int(
                100 * passed / (test * stars_per_test))

        if autosave > 0 and not test % autosave:
            # noinspection SpellCheckingInspection
            test_result.write(
                "sbast-autosave%d.tmp" % self._index, overwrite=True,
                format="fits")

    def _auto_run(
            self,  n_tests: int, stars_per_test: int, sub_image_size: int,
            mag_range: Tuple[int, int], loading_buffer: np.ndarray | None,
            autosave: int, skip_phot: bool | int, skip_background: bool | int,
            zp_mag: float, save_image: bool,
            save_image_path: str, ast_seed: int | None) -> Table | None:
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
        :param zp_mag: the zero point magnitude
        :type zp_mag: float
        :param save_image: bool saying if we should save the added star image.
        :type save_image: bool
        :param save_image_path: path for the saving of the image
        :type save_image_path: str
        :param ast_seed: the seed to set the artificial stars with. or None.
        :type ast_seed: int | None
        :return: Full raw test results. Injected initial properties with
                 measured values.
        :rtype: astropy.table.Table
        """

        test_result: Table = Table(
            np.full((n_tests * stars_per_test, self.N_COLUMNS), np.nan),
            names=self.TEST_TABLE_COLUMN_NAMES)
        scale_factor: float | int = (
            get_mj_ysr2jy_scale_factor(self._starbug.main_image()))

        current_image = self._starbug.image
        assert current_image is not None
        base_image: fits.HDUList = current_image.copy()
        base_shape: np.ndarray = np.copy(self._starbug.main_image().shape)
        stars_per_test: int = int(stars_per_test)
        passed: int = 0
        buffer: int = 0

        if mag_range[0] - mag_range[1] >= 0:
            warn("Detected magnitude range in wrong order,"
                 " put bright limit first\n")
            return None

        if any(base_shape < sub_image_size):
            sub_image_size = min(base_shape)
            p_error("sub image size greater than image size, setting to "
                    "'safe' value %d.\n" % sub_image_size)

        for test in range(1, int(n_tests) + 1):
            self._execute_test(
                base_image, stars_per_test, buffer, mag_range, zp_mag,
                scale_factor, skip_phot, skip_background, passed, test_result,
                test, loading_buffer, autosave, save_image, save_image_path,
                ast_seed)
        return test_result

    def single_test(
            self, contains: Table, skip_phot: bool | int = 0,
            skip_background: bool | int = 0) -> Table:
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
        test_result: Table = Table(
            np.full((len(contains), 4), np.nan),
            names=[TableColumn.X_DET, TableColumn.Y_DET, TableColumn.FLUX_DET,
                   TableColumn.STATUS])

        threshold: Quantity = 2 * units.arcsec

        # Run detection on the image
        end_state: ExitStates
        det: Table | None
        end_state = self._starbug.detect()
        det = self._starbug.detections
        if end_state == ExitStates.EXIT_SUCCESS:
            assert self._starbug is not None
            assert det is not None

            # Check for detection in output
            for i, src in enumerate(contains):  # type: ignore
                separations: np.ndarray = np.sqrt(
                    (src[TableColumn.X_0] - det[TableColumn.X_CENTROID]) ** 2
                    + (src[TableColumn.Y_0] -
                       det[TableColumn.Y_CENTROID]) ** 2) * threshold.unit
                best_match: int = np.argmin(separations)  # noqa
                if separations[best_match] < threshold:
                    test_result[TableColumn.X_DET][i] = (
                        det[TableColumn.X_CENTROID][best_match])
                    test_result[TableColumn.Y_DET][i] = (
                        det[TableColumn.Y_CENTROID][best_match])
                    test_result[TableColumn.FLUX_DET][i] = (
                        det[TableColumn.FLUX][best_match])
                    test_result[TableColumn.STATUS][i] = self.DETECT
                else:
                    test_result[TableColumn.STATUS][i] = self.NULL

            # Run background
            if (sum(test_result[TableColumn.STATUS])
                and (skip_background
                     or not self._starbug.bgd_estimate())):

                # estimate if there were detections
                self._starbug.detections = test_result

                if skip_phot:
                    return hstack((contains, test_result))

                # Run PSF photometry on detected sources
                self._starbug.photometry_routine()
                # noinspection SpellCheckingInspection
                psf_catalogue = self._starbug.psf_catalogue
                assert psf_catalogue is not None
                psf_catalogue.rename_columns(
                    (TableColumn.X_INIT, TableColumn.Y_INIT,
                     TableColumn.XY_DEV),
                    (TableColumn.X_INIT, TableColumn.Y_INIT,
                     TableColumn.XY_DEV_))
                matched: Table = GenericMatch(threshold=threshold)(
                    [contains, psf_catalogue],
                    cartesian=True)
                test_result[TableColumn.FLUX_DET] = (
                    matched[:len(test_result)][TableColumn.FLUX_2])
        return hstack((contains, test_result))


def get_completeness(test_result: Table) -> Table:
    """
    Compile the results into magnitude binned values of recovery fraction
    and flux error.

    :param test_result: The output from auto_run.
    :type test_result: astropy.table.Table
    :return: A table containing per cent completeness as a function of
             magnitude.
    :rtype: astropy.table.Table
    """

    bins: np.ndarray = np.arange(
        np.floor(np.nanmin(test_result[TableColumn.MAG])),
        np.ceil(np.nanmax(test_result[TableColumn.MAG])),
        0.1)
    percents: np.ndarray = np.zeros(len(bins))
    errors: np.ndarray = np.zeros(len(bins))
    offsets: np.ndarray = np.zeros(len(bins))
    means: np.ndarray = np.zeros(len(bins))

    i_bins: np.ndarray = np.asarray(np.digitize(
        test_result[TableColumn.MAG], bins=bins))
    for i in range(max(i_bins)):
        binned: Table = test_result[(i_bins == i)]
        if binned:
            percents[i] = float(sum(binned[TableColumn.STATUS])) / len(binned)

        mag_inj: np.ndarray = -2.5 * np.log10(binned[TableColumn.FLUX])
        mag_det: np.ndarray = -2.5 * np.log10(binned[TableColumn.FLUX_DET])
        errors[i] = np.nanstd(mag_inj - mag_det)
        means[i] = np.nanmean(mag_inj - mag_det)
        offsets[i] = np.nanmedian(
            binned[TableColumn.FLUX] / binned[TableColumn.FLUX_DET])

    out: Table = Table(
        [bins, percents, errors, offsets],
        names=(TableColumn.MAG, TableColumn.REC, TableColumn.ERR_LOWER,
               TableColumn.OFF),
        dtype=(float, float, float, float))
    return out


def get_spatial_completeness(
        test_result: Table, image: np.ndarray | None,
        res: int = 10) -> np.ndarray | None:
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
    if image is None:
        return None

    x_bins: np.ndarray = np.arange(
        min(test_result[TableColumn.X_0]),
        max(test_result[TableColumn.X_0]), int(res))
    y_bins: np.ndarray = np.arange(
        min(test_result[TableColumn.Y_0]),
        max(test_result[TableColumn.Y_0]), int(res))
    percents: np.ndarray = np.zeros(image.shape)

    xi: int
    for xi in x_bins[:-1]:
        xo: int = xi + res
        yi: int
        for yi in y_bins[:-1]:
            yo: int = yi + res
            mask: np.ndarray = (
                (test_result[TableColumn.X_0] >= xi) &
                (test_result[TableColumn.X_0] < xo) &
                (test_result[TableColumn.Y_0] >= yi) &
                (test_result[TableColumn.Y_0] < yo))
            binned: Table = test_result[mask]
            if len(binned):
                percents[int(xi): int(xo), int(yi): int(yo)] = (
                    float(np.sum(binned[TableColumn.STATUS])) / len(binned))
    return percents


def estimate_completeness_mag(ast: Table) -> (
        Tuple[Tuple[float, float, float] | None,
              Tuple[float, float, float] | None]):
    """
    Estimate the completeness level of the artificial star test.

    :param ast: Output of Artificial_Stars.get_completeness, table must
               contain columns (mag, rec).
    :type ast: astropy.table.Table
    :return: A tuple containing:
        - **fit** (*list*): The fitting parameters to the logistic curve
                            $f(x) = \frac{l}{1 + exp(-k(x - x_0))}$ formatted
                            as ``[l, x_0, k]``.
        - **complete** (*list*): Magnitude of 70% and 50% completeness.
    :rtype: tuple[list, list]
    """
    fit: Tuple[float, float, float] | None = None
    completeness: Tuple[float, float, float] | None = None

    # Syntax: Callable[[Param1Type, Param2Type, ...], ReturnType]
    fn_i: Callable[[float, float, float, float], float] = (
        lambda y, limit, k, xo: xo - (np.log((limit / y) - 1) / k)
    )

    if len(set(ast.colnames) & {TableColumn.MAG, TableColumn.REC}) == 2:
        try:
            # need the *_ as the return tuple can be multiple sizes. The *_
            # allows the IDE to not freak out, especially as we don't care
            # about the rest of the return values.
            fit, *_ = curve_fit(
                scurve, ast[TableColumn.MAG], ast[TableColumn.REC],
                [1, -1, np.median(ast[TableColumn.MAG])])
            assert fit is not None
            completeness = (fn_i(0.9, *fit), fn_i(0.7, *fit), fn_i(0.5, *fit))
        except (RuntimeError, ValueError) as e:
            warn(f"Unable to fit completeness fractions: {e}\n")
    else:
        p_error("Input table must have columns 'mag' and 'rec'\n")
    return fit, completeness


def scurve(
        x: np.ndarray, limit: float, k: float,
        xo: float) -> float | np.ndarray:
    """
    S-curve function to fit completeness results to.

    math::  f(x) = \\frac{l}{1 + \\exp(-k(x - x_0))

    :param x: Magnitude range or array to input into the function.
    :type x: list or numpy.ndarray
    :param limit: Maximum value asymptote (typically representing maximum
              completeness, near 1.0).
    :type limit: float
    :param xo: The inflection point of the curve (the magnitude where
               completeness is 50%).
    :type xo: float
    :param k: The logistic growth rate or steepness of the curve.
    :type k: float
    :return: Calculated function value(s) matching the shape of the
             input ``x``.
    :rtype: float or numpy.ndarray
    """
    return limit / (1 + np.exp(-k * (x - xo)))


def compile_results(
        raw: Table,
        image: np.ndarray | None = None,
        plot_ast: str | None = None,
        filter_string: str = "m") -> fits.HDUList:
    """
    Compile all the raw data into usable results

    :param raw:raw data
    :type raw: astro.table.table
    :param image: the image data
    :type image: np.ndarray
    :param plot_ast: the save plot file name
    :type plot_ast: str or None
    :param filter_string: the filter string
    :type filter_string: str
    :return: the results
    :rtype: fits.HDUList
    """

    completeness_raw: Table = get_completeness(raw)
    cfit: Tuple[float, float, float]
    completeness: Tuple[float, float, float]
    cfit, completeness = estimate_completeness_mag(completeness_raw)
    spatial_completeness: np.ndarray | None = (
        get_spatial_completeness(raw, image, res=10))

    head: Dict[str, str | float] = {
        "COMPLETE_FN": "F(x)=l/(1+exp(-k(x-xo)))", "l": cfit[0],
        "k": cfit[1], "xo": cfit[2]}
    for i, frac in enumerate((90, 70, 50)):
        if completeness[i] and not np.isnan(completeness[i]):
            printf(
                "-> complete to %d%%: %s=%.2f\n" % (
                    frac, filter_string, completeness[i]))
            head["COMPLETE %d%%" % frac] = str(completeness[i])

    # needed for spatial_completeness as it expects an 'array.pyi',
    # got 'ndarray' instead.
    results: fits.HDUList = fits.HDUList(
        [fits.PrimaryHDU(header=fits.Header(head)),
         fits.BinTableHDU(data=completeness_raw, name="AST"),
         fits.BinTableHDU(data=raw, name="RAW"),
         fits.ImageHDU(data=cast(Any, spatial_completeness), name="CMP")])

    if plot_ast:
        fig: Figure
        ax: Axes
        fig, ax = plt.subplots(1, figsize=(3.5, 3), dpi=300)
        ax.scatter(
            completeness_raw[TableColumn.MAG],
            completeness_raw[TableColumn.REC], c='k', lw=0, s=8)
        ax.plot(completeness_raw[TableColumn.MAG],
                scurve(completeness_raw[TableColumn.MAG], *cfit),
                c='g',
                label=r"$f(x)=\frac{%.2f}{1+e^{%.2f("r"x-%.2f)}}$" % (
                    cfit[0], cfit[1], cfit[2]))
        ax.axvline(
            completeness[0], c="seagreen", ls='--',
            label=("90%%:%.2f" % completeness[0]), lw=0.75)
        ax.axvline(
            completeness[1], c="seagreen", ls='-.',
            label=("70%%:%.2f" % completeness[1]), lw=0.75)
        ax.axvline(
            completeness[2], c="seagreen", ls=':',
            label=("50%%:%.2f" % completeness[2]), lw=0.75)
        ax.scatter(completeness, (0.9, 0.7, 0.5), marker='*', c='teal', s=10)
        ax.tick_params(direction="in", top=True, right=True)
        ax.set_title("Artificial Star Test")
        ax.set_xlabel(filter_string)
        ax.set_ylabel("Fraction Recovered")
        ax.set_yticks([0, .25, .5, .75, 1])
        ax.legend(loc="lower left", frameon=False, fontsize=8)
        plt.tight_layout()
        fig.savefig(plot_ast, dpi=300)
        printf("--> %s\n" % plot_ast)
    return results
