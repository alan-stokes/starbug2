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
from typing import Final, Tuple, Dict, List, cast, Any

from astropy.wcs import (
    WCS, NoConvergence, SingularMatrixError, InconsistentAxisTypesError,
    InvalidTransformError)
import numpy as np
from astropy.io.fits import (
    PrimaryHDU, ImageHDU, HDUList, Header, open, BinTableHDU)
from astropy.table import hstack, Table
from starbug2.core.constants import (
    HeaderTags, ImageHeaderTags, SCI, BGD, RES, VERBOSE_TAG, AP_FILE, BGD_FILE,
    FITS_EXTENSION, DQ, AREA, WHT, ExitStates, TableColumn)
from starbug2.core.main_components.aperture_photometry import (
    AperturePhotometry)
from starbug2.core.main_components.detect import Detect

from starbug2.core.main_components.photometry import Photometry
from starbug2.utilities.filters import STAR_BUG_FILTERS
from starbug2.routines.background_estimate_routine import (
    BackGroundEstimateRoutine)
from starbug2.routines.source_properties import SourceProperties
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.interfaces.star_bug_interface import StarBugInterface
from starbug2.utilities.utils import (
    collapse_header, get_version, ext_names, printf, split_file_name, p_error,
    warn, import_table, reindex, get_data_path)


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
    def sort_output_names(
        f_name: str, param_output: str | None = None
    ) -> Tuple[str, str, str]:
        """
        This is a useful function that looks at both an input file and a set
        output and figures out how to name output files. If param_output looks
        like a directory then the output will be set to that directory with
        the basename of f_name. If param_output looks like a file, then the
        output basename will take that form.

        :param f_name: Filename to use as the core of the output
        :type f_name: str
        :param param_output: This is the OUTPUT parameter in the parameter
            file. It can be an output directory or output filename. If None
            (default) then it will be ignored
        :type param_output: str
        :return: A tuple of (The output directory, The output file basename,
            The file extension split from the inputs)
        :rtype: tuple of (str, str, str)
        """
        out_dir: str = ""
        b_name: str = ""
        extension: str = ""
        if f_name:
            out_dir, b_name, extension = split_file_name(f_name)
            if (tmp_out_name := param_output) and tmp_out_name != '.':
                inner_out_dir, inner_b_name, _ = split_file_name(tmp_out_name)
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
        ap_file: str | None, bkg_file: str | None
    ) -> None:
        """
        Star bug initialisation.

        :param f_name: FITS image file name
        :type f_name: str
        :param config: The starbug configuration object
        :type config: StarBugMainConfig
        :param ap_file: Optional aperture coordinates file path
        :type ap_file: str or None
        :param bkg_file: Optional background reference file path
        :type bkg_file: str or None
        """
        # Defaults
        self._config = config
        self._f_name: str | None = None
        self._out_dir: str | None = None
        self._b_name: str | None = None
        self._image: HDUList | None = None
        self._filter: str | None = None
        self._header: Header | None = None
        self._wcs: WCS | None = None
        self._stage: float = 0.0
        self._detections: Table | None = None
        self._n_hdu: int = -1
        self._unit: str | None = None
        self._background: ImageHDU | PrimaryHDU | None = None
        self._residuals: np.ndarray | None = None
        self._psf_catalogue: Table | None = None
        self._source_stats: np.ndarray | None = None
        self._psf: np.ndarray | None = None

        # Overridden configs
        self._ap_file: str | None = ap_file
        self._background_file: str | None = bkg_file
        self._full_width_half_max: float = config.full_width_half_max

        # Process options
        # Load the fits image
        self.load_image(f_name)

        if ap_file is not None:
            # Load the source list if given
            self.load_ap_file(ap_file)
        if bkg_file is not None:
            self.load_bgd_file(bkg_file)

    def log(self, msg: str) -> None:
        """
        Print message if in verbose mode.

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
            # Sorting out the file names and what not
            extension: str
            self._out_dir, self._b_name, extension = self.sort_output_names(
                f_name, self._config.output_file)

            if extension == FITS_EXTENSION:
                if os.path.exists(f_name):
                    self.log("loaded: \"%s\"\n" % f_name)
                    self._image = open(f_name)

                    # Force assigning _nHDU
                    main_image: ImageHDU | PrimaryHDU = self.main_image()
                    self._header = main_image.header

                    self.log(
                        "-> using image HDU: %d (%s)\n" % (
                            self._n_hdu, main_image.name))

                    if main_image.data is None:
                        warn("Image seems to be empty.\n")

                    if ((val := main_image.header.get(
                            ImageHeaderTags.TELESCOPE)) is None
                            or (val.find(ImageHeaderTags.JWST) < 0)):
                        warn("Telescope not JWST, "
                             "there may be undefined behaviour.\n")

                    self._filter = self._config.custom_filter
                    assert self._filter is not None
                    if ((HeaderTags.FILTER in main_image.header) and
                        (main_image.header[HeaderTags.FILTER] in
                         STAR_BUG_FILTERS.keys())):
                        self._filter = main_image.header[HeaderTags.FILTER]
                        assert self._filter is not None
                        if self._full_width_half_max < 0:
                            self._full_width_half_max = (
                                STAR_BUG_FILTERS[
                                    self._filter].full_width_half_max)
                    if self._filter:
                        self.log("-> photometric band: %s\n" % self._filter)
                    else:
                        warn("Unable to determine image filter\n")

                    if ImageHeaderTags.DETECTOR in self.info.keys():
                        self.log(
                            "-> detector module: %s\n" %
                            self.info[ImageHeaderTags.DETECTOR])
                    else:
                        warn("Unable to determine Telescope DETECTOR.\n")

                    if ImageHeaderTags.BUN_IT in main_image.header:
                        self._unit = main_image.header[ImageHeaderTags.BUN_IT]
                    else:
                        warn("Unable to determine image BUNIT.\n")

                    self._wcs = WCS(self.main_image().header)

                    # Determine calculation stage level
                    extension_names: List[str] = ext_names(self._image)
                    if DQ in extension_names:
                        if AREA in extension_names:
                            self._stage = 2.0
                        else:
                            self._stage = 2.5
                    elif WHT in extension_names:
                        self._stage = 3.0
                    elif HeaderTags.CALIBRATION_LV in self.main_image().header:
                        self._stage = (self.main_image().header[
                            HeaderTags.CALIBRATION_LV])
                    else:
                        warn("Unable to determine calibration level, "
                             "assuming stage 3\n")
                        self._stage = 3.0
                    self.log("-> pipeline stage: %d\n" % self._stage)

                else:
                    warn("fits file \"%s\" does not exist\n" % f_name)
            else:
                warn("included file must be FITS format\n")

    def load_ap_file(self, f_name: str | None = None) -> None:
        """
        Load an AP_FILE to be used during photometry.

        :param f_name: Filename for fits table containing source coordinates.
            These coordinates can be x-centroid / y-centroid, x_init / y_init,
            x_0, y_0 or RA/DEC. The latter is used if starbug gets "USE_WCS=1"
            in the parameter file.
        :type f_name: str
        :return: None
        """
        if not f_name:
            f_name = self._ap_file
        if f_name and os.path.exists(f_name):
            self._detections = import_table(f_name)
            if self._detections is None or self._wcs is None:
                raise Exception("could not read the ap file")

            column_names: set[str] = set(self._detections.colnames)

            self.log("loaded AP_FILE='%s'\n" % f_name)

            if self._config.use_wcs_values:
                if len(column_names & {TableColumn.RA, TableColumn.DEC}) == 2:
                    self.log("-> using RA-DEC coordinates\n")
                    try:
                        xy: Any = self._wcs.all_world2pix(
                            self._detections[TableColumn.RA],
                            self._detections[TableColumn.DEC], 0)
                    except (NoConvergence, MemoryError, SingularMatrixError,
                            InconsistentAxisTypesError, ValueError,
                            InvalidTransformError) as e:
                        warn(f"Something went wrong converting WCS to pixels "
                             f"({e}), trying wcs_world2pix next.\n")
                        xy = self._wcs.wcs_world2pix(
                            self._detections[TableColumn.RA],
                            self._detections[TableColumn.DEC], 0)
                    if TableColumn.X_CENTROID in column_names:
                        self._detections.remove_column(TableColumn.X_CENTROID)
                    if TableColumn.Y_CENTROID in column_names:
                        self._detections.remove_column(TableColumn.Y_CENTROID)
                    self._detections.add_columns(
                        xy,
                        names=[TableColumn.X_CENTROID, TableColumn.Y_CENTROID],
                        indexes=[0, 0])
                else:
                    warn("No 'RA' or 'DEC' found in AP_FILE\n")

            elif len({TableColumn.X_0, TableColumn.Y_0} & column_names) == 2:
                self._detections.rename_columns(
                    (TableColumn.X_0, TableColumn.Y_0),
                    (TableColumn.X_CENTROID, TableColumn.Y_CENTROID))
            elif (len({TableColumn.X_INIT, TableColumn.Y_INIT}
                      & column_names) == 2):
                self._detections.rename_columns(
                    (TableColumn.X_INIT, TableColumn.Y_INIT),
                    (TableColumn.X_CENTROID, TableColumn.Y_CENTROID))

            if len({TableColumn.X_CENTROID, TableColumn.Y_CENTROID} &
                   set(self._detections.colnames)) == 2:
                mask: np.ndarray = (
                    (self._detections[TableColumn.X_CENTROID] >= 0)
                    & (self._detections[TableColumn.X_CENTROID] <
                       self.main_image().shape[1])
                    & (self._detections[TableColumn.Y_CENTROID] >= 0)
                    & (self._detections[TableColumn.Y_CENTROID] <
                       self.main_image().shape[0])
                )

                # cant figure how to resolve this typing
                self._detections.remove_rows(~mask)  # noqa
                self.log(
                    "-> loaded %d sources from AP_FILE\n" %
                    len(self._detections))
            else:
                warn("Unable to determine physical coordinates"
                     " from detections table\n")
        else:
            p_error("AP_FILE='%s' does not exist\n" % f_name)

    def load_bgd_file(self, f_name: str | None = None) -> None:
        """
        Load a BGD_FILE to be used during photometry.

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

    def detect(self) -> ExitStates:
        """
        Full source detection routine. Saves the result as a table.
        self._detections

        :return: Status
        :rtype: ExitStates
        """
        main_detect: Detect = Detect()
        detections: Table | None
        detect_exit_states: ExitStates
        (detect_exit_states, detections) = main_detect.detect(
            self.log, self.main_image, self._filter, self._config, self._wcs,
            self.header)
        self._detections = detections
        return detect_exit_states

    # noinspection SpellCheckingInspection
    def aperture_photometry(self) -> ExitStates:
        """
        Executes aperture photometry processing steps.

        :return: Success or failure exit state identifier
        :rtype: ExitStates
        """
        main_aperture_photometry: AperturePhotometry = AperturePhotometry()
        detections: Table
        result: ExitStates

        (detections, result) = main_aperture_photometry.aperture_photometry(
            self.detections, self.image, self._filter, self.log, self._config,
            self.info, self._full_width_half_max, self.out_dir, self._b_name,
            self.header(), self._background, self.main_image())
        self._detections = detections
        return result

    def bgd_estimate(self) -> ExitStates:
        """
        Estimate the background of the active image.
        Saves the result as an ImageHDU self._background

        :return: The status execution code
        :rtype: ExitStates
        """
        self.log("\nEstimating Diffuse Background\n")
        status: ExitStates = ExitStates.EXIT_SUCCESS
        assert self._filter is not None
        if self._detections:
            source_list: Table = self._detections.copy()

            # noinspection DuplicatedCode
            if TableColumn.X_INIT in source_list.colnames:
                source_list.rename_column(
                    TableColumn.X_INIT, TableColumn.X_CENTROID)
            if TableColumn.Y_INIT in source_list.colnames:
                source_list.rename_column(
                    TableColumn.Y_INIT, TableColumn.Y_CENTROID)
            if TableColumn.X_DET in source_list.colnames:
                source_list.rename_column(
                    TableColumn.X_DET, TableColumn.X_CENTROID)
            if TableColumn.Y_DET in source_list.colnames:
                source_list.rename_column(
                    TableColumn.Y_DET, TableColumn.Y_CENTROID)
            if TableColumn.FLUX_DET in source_list.colnames:
                source_list.rename_column(
                    TableColumn.FLUX_DET, TableColumn.FLUX)
            mask: np.ndarray = ~(
                np.isnan(source_list[TableColumn.X_CENTROID])
                | np.isnan(source_list[TableColumn.Y_CENTROID]))

            bgd: BackGroundEstimateRoutine = BackGroundEstimateRoutine(
                source_list[mask],
                box_size=self._config.background_box_size,
                full_width_half_max=(
                    self._config.full_width_half_max_with_filter(
                        STAR_BUG_FILTERS.get(self._filter))),
                sig_sky=self._config.sigma_sky,
                bgd_r=self._config.bgd_radius,
                profile_scale=self._config.profile_scaling_factor,
                profile_slope=self._config.profile_slope,
                verbose=self._config.verbose_logs)
            header: Header = self.header()

            # Check for insanity
            if self._wcs is None:
                return ExitStates.EXIT_FAIL

            header.update(self._wcs.to_header())

            # Get image data
            image_data = bgd(
                self.main_image().data.copy(),
                output=self._config.bgd_check_file)
            assert image_data is not None

            self._background = ImageHDU(
                data=image_data.background,
                header=header)

            # Check for insanity
            if self._background is None:
                return ExitStates.EXIT_FAIL

            f_name = "%s/%s-bgd.fits" % (self._out_dir, self._b_name)
            self.log("--> %s\n" % f_name)
            self._background.writeto(f_name, overwrite=True)

        else:
            p_error("unable to estimate background, no source list loaded\n")
            status = ExitStates.EXIT_FAIL
        return status

    def bgd_subtraction(self) -> ExitStates:
        """
        Internally subtract a background array from an image array.

        :return: Success or failure exit state
        :rtype: ExitStates
        """
        self.log("Subtracting Background\n")

        if self._background is None or self._wcs is None:
            p_error("No background array loaded (-b file-bgd.fits)\n")
            return ExitStates.EXIT_FAIL

        array: np.ndarray = self.main_image().data - self._background.data
        self._residuals = array

        assert self._image is not None
        self._image[self._n_hdu].data = array
        header: Header = self.header()
        header.update(self._wcs.to_header())

        # Wrap array as Any to map properly to ImageHDU writer expectations
        hdu_output = ImageHDU(data=cast(Any, array), header=header, name="RES")
        f_name = "%s/%s-res.fits" % (self._out_dir, self._b_name)

        try:
            hdu_output.writeto(f_name, overwrite=True)
            self.log(
                "--> Background subtracted image written to %s\n" % f_name)
        except Exception as e:
            p_error(
                "Failed to write background-subtracted file: %s\n" % str(e))
            return ExitStates.EXIT_FAIL

        return ExitStates.EXIT_SUCCESS

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
        star_bug_photometry: Photometry = Photometry()
        result: ExitStates
        psf_catalogue: Table | None
        residuals: np.ndarray | None
        (result, psf_catalogue, residuals) = (
            star_bug_photometry.photometry_routine(
                self._filter, self._wcs, self._config, self.main_image(),
                self.log, self.image, self.info,  self._background,
                self.header(), self._detections, self._psf,
                self._full_width_half_max, self._ap_file,
                self._background_file, self._out_dir, self._b_name))
        self._psf_catalogue = psf_catalogue
        self._residuals = residuals
        return result


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
            self.main_image().data, slist,
            verbose=self._config.verbose_logs)
        stat: Table = sp.execute_source_props(
            full_width_half_max=STAR_BUG_FILTERS[
                self._filter].full_width_half_max,
            do_crowd=self._config.calculate_crowding_metric)

        self._source_stats = hstack((slist, stat))
        f_name: str = "%s/%s-stat.fits" % (self._out_dir, self._b_name)
        self.log("--> %s\n" % f_name)
        reindex(Table(self._source_stats))
        BinTableHDU(
            data=self._source_stats, header=self.header()).writeto(
            f_name, overwrite=True)

    # noinspection SpellCheckingInspection
    def verify(self) -> ExitStates:
        """
        This simple function verifies that everything necessary has been
        loaded properly

        :return: int where 0 on success, 1 on fail
        :rtype ExitStates
        """

        status: ExitStates = ExitStates.EXIT_SUCCESS

        self.log("Checking internal systems..\n")

        if not self._filter:
            warn("No FILTER set, please set in parameter file or "
                 "use \"-s FILTER=XXX\"\n")
            status = ExitStates.EXIT_FAIL

        d_name: str = os.path.expandvars(get_data_path())
        if not os.path.exists(d_name):
            warn("Unable to locate STARBUG_DATDIR='%s'\n" % d_name)

        if self._out_dir is not None and not os.path.exists(self._out_dir):
            warn("Unable to locate OUTPUT='%s'\n" % self._out_dir)
            status = ExitStates.EXIT_FAIL

        if self._image is None or self.main_image().data is None:
            warn("Image did not load correctly\n")
            status = ExitStates.EXIT_FAIL

        if self._ap_file and self._detections is not None:
            test = self._filter_detections()
            if not len(test):
                warn("Detection file empty or no sources overlap the image.\n")
                status = ExitStates.EXIT_FAIL

        return status

    def _filter_detections(self) -> Table:
        """
        filters the detections based on some fixed constraints.
        :return: the filtered detections
        """
        assert self._detections is not None
        detections: Table = self._detections[
            [TableColumn.X_CENTROID, TableColumn.Y_CENTROID]].copy()
        detections = detections[detections[TableColumn.X_CENTROID] >= 0]
        detections = detections[detections[TableColumn.Y_CENTROID] >= 0]
        detections = detections[
            detections[TableColumn.X_CENTROID] <
            self.main_image().header[HeaderTags.NAXIS1]]
        return detections[
            detections[TableColumn.Y_CENTROID] <
            self.main_image().header[HeaderTags.NAXIS2]]

    def __getstate__(self) -> dict[str, Any]:
        """
        extracts the inner state of this class. deleting image or/and
         background if it's there.
        :return: the internal state with those bits filtered away
        """
        assert self._image is not None
        self._image.close()
        state: dict[str, Any] = self.__dict__.copy()
        if "_image" in state:
            # Sorry but we cant have that
            del state["_image"]
            # This currently doesnt get reloaded
        if "_background" in state:
            del state["_background"]

        return state

    def __setstate__(self, state) -> None:
        self.__dict__.update(state)
        v: bool = self._config.verbose_logs
        self._config.unfreeze()
        self._config.verbose_logs = False
        self._config.freeze()
        self.load_image(self._f_name)
        self._config.unfreeze()
        self._config.verbose_logs = v
        self._config.freeze()

    def header(self) -> Header:
        """
        Construct relevant base header information for routine products

        :return:  Header file containing a series of relevant information
        :rtype: Header
        """
        head: Dict[str, str | float | None] = {
            HeaderTags.STAR_BUG: get_version(),
            HeaderTags.CALIBRATION_LV: self._stage
        }

        if self._filter:
            head[HeaderTags.FILTER] = self._filter

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
        head[VERBOSE_TAG] = self._config.verbose_logs

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
            ImageHeaderTags.FILTER, ImageHeaderTags.DETECTOR,
            ImageHeaderTags.TELESCOPE, ImageHeaderTags.INSTRUMENT,
            ImageHeaderTags.BUN_IT, ImageHeaderTags.PIXAR_A2,
            ImageHeaderTags.PIXAR_SR]
        if self._image:
            for hdu in self._image:
                out.update(
                    {(key, hdu.header[key]) for key in keys
                     if key in hdu.header})
        return out

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

        # HDU_NAME in param file
        n: str = str(self._config.hdu_name)
        if n and n in e_names:
            self._n_hdu = e_names.index(n)
            return self._image[n]

        # index?
        if isinstance(n, (int, float, np.number)):
            self._n_hdu = int(n)
            return self._image[self._n_hdu]

        # SCI, BGD, RES (common names)
        for name in (SCI, BGD, RES):
            name: str
            if name in e_names:
                self._n_hdu = e_names.index(name)
                return self._image[name]

        # First ImageHDU
        # ABS ARE WE SURE WE WANT TO LOOK FOR A INDEX WITH A ENUMERATE INDEX?
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
