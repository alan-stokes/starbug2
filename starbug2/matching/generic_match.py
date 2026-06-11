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

"""
Starbug matching functions
Primarily this is the main routines for dither/band/generic matching which are
 at the core of starbug2 and starbug2-match
"""
from typing import Any
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, Column, vstack
from starbug2.constants import (
    VERBOSE_TAG, CAT_NUM, FILTER, MATCH_THRESH, RA, DEC, SRC_GOOD, SRC_VAR,
    STD_FLUX, E_FLUX, FLUX, FLAG, NUM)
from starbug2.param import load_params
from starbug2.utils import (
    Loading, printf, remove_duplicates, p_error, fill_nan, tab2array,
    find_col_names, flux2mag)


class GenericMatch:

    @staticmethod
    def build_meta(catalogues: list[Table]) -> Table:
        """
        Extracts structural tracking headers to set up baseline combined
        outputs.
        """
        meta: dict[str, Any] = catalogues[0].meta
        base: Table = Table(None, meta=meta)
        return base

    @staticmethod
    def mask_catalogues(
            catalogues: list[Table],
            mask: list[np.ndarray | list[Any]] | np.ndarray | None) -> Table:
        """ Takes catalogues and masks and removes catalogues which don't
        match the mask.

        :param catalogues: the catalogues to match.
        :type catalogues: list (astropy.Tables)
        :param mask: the mask to apply.
        :return: an astro table with masked catalogues.
        :rtype: astropy.Table
        """
        masked: Table = Table(None)

        if mask is None or not isinstance(mask, (list, np.ndarray)):
            return masked
        if len(catalogues) != len(mask):
            return masked

        for subset, cat in zip(mask, catalogues):
            if subset is not None:
                if isinstance(subset, list):
                    subset = np.array(subset)
                if len(subset) == len(cat):
                    masked = vstack((masked, cat[~subset]))
                    cat.remove_rows(~subset)
        return masked

    @staticmethod
    def _sky_coords_not_cartesian(base: Table) -> SkyCoord:
        """
        create sky coords which are not cartesian

        :param base: the base.
        :type base: astropy.table.Table
        :return: a sky coords sing base values.
        :rtype: SkyCoord.
        """
        ra_cols: list[str] = list(
            name for name in base.colnames if RA in name)
        dec_cols: list[str] = list(
            name for name in base.colnames if DEC in name)
        ra: np.ndarray = np.nanmean(
            tab2array(base, col_names=ra_cols), axis=1)
        dec: np.ndarray = np.nanmean(
            tab2array(base, col_names=dec_cols), axis=1)
        return SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    @staticmethod
    def _sky_coords_cartesian(base: Table) -> SkyCoord:
        """
        create sky coords which are cartesian

        :param base: the base.
        :type base: astropy.table.Table
        :return: a sky coords sing base values.
        :rtype: SkyCoord.
        """
        x_cols: list[str] = list(
            name for name in base.colnames if name[0] == "x")
        y_cols: list[str] = list(
            name for name in base.colnames if name[0] == "y")
        x: np.ndarray = np.nanmean(tab2array(base, col_names=x_cols), axis=1)
        y: np.ndarray = np.nanmean(tab2array(base, col_names=y_cols), axis=1)
        return SkyCoord(
            x=x, y=y, z=np.zeros(len(x)), representation_type="cartesian")

    def __init__(
        self,
        threshold: float | None = None,
        col_names: list[str] | None = None,
        filter_string: str | None = None,
        verbose: int | None = None,
        p_file: str | None = None,
        method: str = "Generic Matching",
        load: Loading = Loading(1)
    ) -> None:
        """
        constructor for the generic match.

        :param threshold: Separation threshold in arc-seconds
        :type threshold: float or None
        :param col_names: List of str column names to include in the matching.
        Everything else will be discarded.
        :type col_names: list or None
        :param filter_string: Specifically set the filter of the catalogues
        :type filter_string: str or None
        :param verbose: Include verbose outputs
        :type verbose: int or None
        :param p_file: Parameter filename
        :type p_file: str or None
        :param load: the loading object
        :type load: Loading
        """
        options: dict[str, float | int | str] = load_params(p_file)

        self._threshold: float = options.get(MATCH_THRESH)
        self._filter: str | None = options.get(FILTER)
        self._verbose: int | None = options.get(VERBOSE_TAG)
        self.method: str = method

        if threshold is not None:
            self._threshold = threshold
        self._threshold *= u.arcsec

        if filter_string is not None:
            self._filter = filter_string
        if verbose is not None:
            self._verbose = verbose

        self._col_names: list[str] | None = col_names
        self._load: Loading = load

    def log(self, msg: str) -> None:
        """
        logs messages only when in verbose mode.
        :param msg: message to log
        :return: None
        """
        if self._verbose:
            printf(msg)

    def __str__(self) -> str:
        """
        string representation fo the generic match class.
        :return: str
        """
        s: list[str] = [
            f"{self.method}:",
            f"Filter: {self._filter}",
            f"Col names: {self._col_names}",
            f'Threshold: {self._threshold}"'
        ]
        return "\n".join(s)

    def __call__(self, *args: Any, **kwargs: Any) -> Table:
        """
        main entrance method.

        :param args: args to the matching algorithm.
        :param kwargs: extra args.
        :return: matched catalogue
        :rtype: astropy.Table
        """
        return self.match(*args, **kwargs)

    def init_catalogues(self, catalogues: list[Table]) -> list[Table]:
        # noinspection SpellCheckingInspection
        """
        This function is a bit of a "do everything" function

        It takes the input catalogues and removes any columns that aren't
        included in the self.colnames list. If this is None, then all columns
        are kept. Additionally, it initialises the loading bar with the summed
        length of all the input catalogues.
        Finally, it attempts to set the photometric filter that is being used

        :param catalogues: The input catalogues to work on
        :type catalogues: list (astropy.Tables)
        :return: The cleaned list of input catalogues
        :rtype: list (astropy.Table)
        """
        if len(catalogues) >= 2:
            self._load = Loading(
                sum(len(cat) for cat in catalogues[1:]), msg="initialising")
            if self._verbose:
                self._load.show()

        # initialise the column names if it wasn't already set
        if self._col_names is None:
            self._col_names = []
            for cat in catalogues:
                self._col_names += cat.colnames
        self._col_names = remove_duplicates(self._col_names)

        # clean out the column names not included in self._col_names
        for n, catalogue in enumerate(catalogues):
            keep: list[str] = list(
                set(catalogue.colnames) & set(self._col_names))
            keep = sorted(
                keep,
                key=lambda s: self._col_names.index(s) if
                              self._col_names else 0)
            catalogues[n] = catalogue[keep]

        if not self._filter:
            filter_string: str | None = catalogues[0].meta.get(FILTER)
            if filter_string is None:
                filter_string = "MAG"
            self._filter = filter_string

        return catalogues


    def match(
            self,
            catalogues: list[Table],
            join_type: str = "or",
            mask: list[np.ndarray | list[Any]] | np.ndarray | None = None,
            cartesian: bool = False,
            **kwargs: Any) -> Table:
        """
        This matching works as a basic match. Everything is included and the
        column names have _N appended to the end.

        :param catalogues: the catalogues to match
        :type catalogues: list (astropy.Tables)
        :param join_type: Joining method ("or" to include sources in any
                          catalogue, "and" to only include sources in all
                          catalogues)
        :type join_type: str
        :param mask: In preparation
        :type mask: list
        :param cartesian: if we should use cartesian coordinates
        :type cartesian: bool
        :param kwargs: extra args
        :type kwargs: dict
        :return: Matched catalogue
        :rtype: astropy.table.Table
        """
        catalogues = self.init_catalogues(catalogues)
        if self._col_names and CAT_NUM in self._col_names:
            self._col_names.remove(CAT_NUM)

        masked: Table = self.mask_catalogues(catalogues, mask)
        base: Table = self.build_meta(catalogues)

        if join_type == "and":
            p_error("join_type 'and' not fully implemented\n")

        # Bulk matching processes (column naming)
        for n, cat in enumerate(catalogues, 1):
            self._load.msg = "matching: %d" % n
            tmp: Table = self.inner_match(
                base, cat, join_type=join_type, cartesian=cartesian)
            tmp.rename_columns(
                tmp.colnames, [f"{name}_{n}" for name in tmp.colnames])
            base = fill_nan(hstack((base, tmp)))

        # Add in any masked bits
        if len(masked):
            masked.rename_columns(
                masked.colnames, [f"{name}_0" for name in masked.colnames])
            base = fill_nan(vstack((base, masked)))
        return base

    def inner_match(
            self, base: Table, cat: Table, join_type: str = "or",
            cartesian: bool = False) -> Table:
        """
        Base matching function between two catalogues

        :param base: The first catalogue for matching
        :type base: astropy.table.Table
        :param cat: The second catalogue for matching
        :type cat: astropy.table.Table
        :param join_type: Joining method ("or" to include all sources, "and"
                          for intersections)
        :type join_type: str
        :param cartesian: Whether to use Cartesian coordinates for matching
        :type cartesian: bool
        :return: Indices, 2D separation, and 3D separation
        :rtype: astropy.table.Table
        """
        if not len(base):
            return cat.copy()

        base = fill_nan(base.copy())
        assert self._col_names is not None
        col_names: list[str] = [
            n for n in self._col_names if n in cat.colnames]
        cat = fill_nan(cat[col_names].copy())

        sky_coords_1: SkyCoord
        sky_coords_2: SkyCoord
        if not cartesian:
            sky_coords_1 = self._sky_coords_not_cartesian(base)
            sky_coords_2 = self._sky_coords_not_cartesian(cat)
        else:
            sky_coords_1 = self._sky_coords_cartesian(base)
            sky_coords_2 = self._sky_coords_cartesian(cat)

        idx: np.ndarray
        d2d: u.Quantity
        d3d: u.Quantity
        idx, d2d, d3d = sky_coords_2.match_to_catalog_3d(sky_coords_1)

        tmp: Table = Table(
            np.full((len(base), len(col_names)), np.nan),
            names=col_names, dtype=cat[col_names].dtype)

        dist: np.ndarray | u.Quantity
        threshold: Any
        if cartesian:
            dist = d3d
            # If your threshold has an explicit unit wrapper,
            # extract its scalar magnitude
            threshold = self._threshold.value if hasattr(
                self._threshold, "value") else self._threshold
        else:
            dist = d2d
            threshold = self._threshold

        src: Any
        IDX: int
        sep: Any
        for src, IDX, sep in zip(cat.as_array(), idx, dist):
            self._load()
            if self._verbose:
                self._load.show()

            if (sep <= threshold) and (sep == min(dist[idx == IDX])):
                tmp[IDX] = src
            elif join_type == "or":
                tmp.add_row(src)

        return tmp

    def finish_matching(
        self,
        tab: Table,
        error_column: str = E_FLUX,
        num_thresh: int = -1,
        zp_mag: float = 0.0,
        col_names: list[str] | None = None) -> Table:
        """
        Averaging all the values. Combining source flags and building a NUM
        column

        :param tab: Table to work on
        :type tab: astropy.table.Table
        :param error_column: Column containing resultant photometric errors
                             (E_FLUX or STD_FLUX)
        :type error_column: str
        :param num_thresh: Minimum number of matches a source must have
                           (no cropping if <= 0)
        :type num_thresh: int
        :param zp_mag: Zero point (Magnitude) to be applied to the final
                       magnitude
        :type zp_mag: float
        :param col_names: List of column names to average; uses
                          self.col_names if None
        :type col_names: list, optional
        :return: An averaged version of the input table
        :rtype: astropy.table.Table
        """
        flags: np.ndarray = np.full(len(tab), SRC_GOOD, dtype=np.uint16)
        av: Table = Table(None)

        if col_names is None:
            col_names = self._col_names if self._col_names else []

        for ii, name in enumerate(col_names):
            if all_cols := find_col_names(tab, name):
                ar: np.ndarray = tab2array(tab, col_names=all_cols)
                col: Column

                if ar.shape[1] > 1:
                    if name == FLUX:
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                        mean: np.ndarray = np.nanmean(ar, axis=1)

                        if self._col_names and STD_FLUX not in self._col_names:
                            av.add_column(
                                Column(np.nanstd(ar, axis=1), name=STD_FLUX),
                                index=ii + 1)
                        ## if median and mean are >5% different, flag as
                        # SRC_VAR
                        flags[np.abs(mean - col) > (col / 5.0)] |= SRC_VAR
                    elif name == E_FLUX:
                        col = Column(
                            np.sqrt(np.nansum(ar * ar, axis=1)), name=name)
                    elif name == STD_FLUX:
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                    elif name == FLAG:
                        col = Column(flags, name=name)
                        f_col: np.ndarray
                        for f_col in ar.T:
                            flags |= f_col.astype(np.uint16)
                    elif name == NUM:
                        col = Column(np.nansum(ar, axis=1), name=name)
                    elif name == CAT_NUM:
                        col = Column(all_cols[0], name=name)
                    else:
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                else:
                    col = tab[all_cols[0]]
                    col.name = name
                av.add_column(col, index=ii)

        av[FLAG] = Column(flags, name=FLAG)
        if FLUX in av.colnames:
            ecol: Column | None = (
                av[error_column] if error_column in av.colnames else None)
            mag: np.ndarray
            mag_err: np.ndarray
            mag, mag_err = flux2mag(av[FLUX], flux_err=ecol)
            mag += zp_mag

            if self._filter in av.colnames:
                av.remove_column(str(self._filter))
            if f"e{self._filter}" in av.colnames:
                av.remove_column(f"e{self._filter}")
            av.add_column(mag, name=str(self._filter))
            av.add_column(mag_err, name=f"e{self._filter}")

        if NUM not in av.colnames:
            narr: np.ndarray = np.nansum(np.invert(
                np.isnan(tab2array(tab, find_col_names(tab, RA)))), axis=1)
            av.add_column(Column(narr, name=NUM))

            if num_thresh > 0:
                av.remove_rows(av[NUM] < num_thresh)
        return av

    @property
    def col_names(self) -> list[str] | None:
        return self._col_names

    @property
    def filter(self) -> str | None:
        return self._filter

    @property
    def threshold(self) -> Any:
        return self._threshold

    @property
    def verbose(self) -> int | None:
        return self._verbose