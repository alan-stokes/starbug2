"""
Starbug matching functions
Primarily this is the main routines for dither/band/generic matching which are
 at the core of starbug2 and starbug2-match
"""
from abc import ABC
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.table import Table, hstack, Column, vstack

import starbug2
from starbug2.constants import (
    VERBOSE_TAG, CAT_NUM, FILTER, MATCH_THRESH, RA, DEC, SRC_GOOD, SRC_VAR)
from starbug2.param import load_params
from starbug2.utils import (
    Loading, printf, remove_duplicates, p_error, fill_nan, tab2array,
    find_col_names, flux2mag)

# keys for catalogue fields.
# noinspection SpellCheckingInspection
_OBS = "OBSERVTN"
_VISIT = "VISIT"
_EXPOSURE = "EXPOSURE"


class GenericMatch(ABC):

    @staticmethod
    def build_meta(catalogues):
        """
        Not happy with this yet
        """
        meta = catalogues[0].meta
        base = Table(None, meta=meta)
        return base

    @staticmethod
    def mask_catalogues(catalogues, mask):
        """ takes catalogues and masks and removes catalogues which don't
        match the mask.

        :param catalogues: the catalogues to match.
        :type catalogues: list (astropy.Tables)
        :param mask: the mask to apply.
        :return: an astro table with masked catalogues.
        :rtype astropy.Table
        """
        masked = Table(None)

        if mask is None or type(mask) not in (list, np.ndarray):
            return masked
        if len(catalogues) != len(catalogues):
            return masked

        for subset, cat in zip(mask, catalogues):
            if subset is not None:
                if type(subset) == list:
                    subset=np.array(subset)
                if len(subset) == len(cat):
                    masked = vstack((masked,cat[~subset]))
                    cat.remove_rows(~subset)
        return masked


    @staticmethod
    def _sky_coords_not_cartesian(base):
        """
        create sky coords which are not cartesian

        :param base: the base.
        :type base: astropy.table.Table
        :return: a sky coords sing base values.
        :rtype: SkyCoord.
        """
        ra_cols = list(name for name in base.colnames if RA in name)
        dec_cols= list(name for name in base.colnames if DEC in name)
        ra = np.nanmean(tab2array(base, col_names=ra_cols), axis=1)
        dec = np.nanmean(tab2array(base, col_names=dec_cols), axis=1)
        return SkyCoord(ra=ra * u.deg, dec=dec * u.deg)

    @staticmethod
    def _sky_coords_cartesian(base):
        """
        create sky coords which are cartesian

        :param base: the base.
        :type base: astropy.table.Table
        :return: a sky coords sing base values.
        :rtype: SkyCoord.
        """
        x_cols = list(name for name in base.colnames if name[0] == "x")
        y_cols = list(name for name in base.colnames if name[0] == "y")
        x = np.nanmean(tab2array(base, col_names=x_cols), axis=1)
        y = np.nanmean(tab2array(base, col_names=y_cols), axis=1)
        return SkyCoord(
            x=x, y=y, z=np.zeros(len(x)), representation_type="cartesian")


    def __init__(
        self, threshold=None, col_names=None, filter_string=None,
        verbose=None, p_file=None, method="Generic Matching",
        load=Loading(1)):
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
        options = load_params(p_file)
        self._threshold = options.get(MATCH_THRESH)
        self._filter = options.get(FILTER)
        self._verbose = options.get(VERBOSE_TAG)
        self.method = method

        if threshold is not None:
            self._threshold = threshold
        self._threshold *= u.arcsec

        if filter_string is not None:
            self._filter = filter_string
        if verbose is not None:
            self._verbose = verbose

        self._col_names = col_names
        self._load = load

    def log(self, msg):
        """
        logs messages only when in verbose mode.
        :param msg: message to log
        :return: None
        """
        if self._verbose:
            printf(msg)

    def __str__(self):
        """
        string representation fo the generic match class.
        :return: str
        """
        s=[ "%s:" % self.method,
            "Filter: %s" % self._filter,
            "Col names: %s" % self._col_names,
            "Threshold: %s\"" % self._threshold]
        return "\n".join(s)

    def __call__(self, *args, **kwargs):
        """
        main entrance method.

        :param args: args to the matching algorithm.
        :param kwargs: extra args.
        :return: matched catalogue
        :rtype: astropy.Table
        """
        return self.match(*args, **kwargs)

    def init_catalogues(self, catalogues):
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

        ## Must copy here maybe?
        if len(catalogues) >= 2:
            self._load = Loading(
                sum(len(cat) for cat in catalogues[1:]), msg="initialising")
            if self._verbose:
                self._load.show()

        # initialise the column names if it wasn't already set
        if self._col_names is None:
            self._col_names = []
            for cat in catalogues:
                self._col_names += cat.col_names
        self._col_names = remove_duplicates(self._col_names)

        # clean out the column names not included in self._col_names
        for n,catalogue in enumerate(catalogues):
            keep = set(catalogue.col_names) & set(self._col_names)
            keep = sorted(keep, key= lambda s: self._col_names.index(s))
            catalogues[n] = catalogue[keep]

        # Attempt to get a value for filter if not already set
        if not self._filter:
            if (filter_string := catalogues[0].meta.get(FILTER)) is None:
                filter_string = "MAG"
            self._filter = filter_string

        return catalogues


    def match(
            self, catalogues, join_type="or", mask=None, cartesian=False,
            **kwargs):
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
        if CAT_NUM in self._col_names:
            self._col_names.remove(CAT_NUM)
        masked = self.mask_catalogues(catalogues, mask)
        base = self.build_meta(catalogues)

        if join_type == "and":
            p_error("join_type 'and' not fully implemented\n")

        # Bulk matching processes (column naming)
        for n, cat in enumerate(catalogues, 1):
            self._load.msg = "matching: %d" % n
            tmp = self.inner_match(
                base, cat, join_type=join_type, cartesian=cartesian)
            tmp.rename_columns(
                tmp.colnames, ["%s_%d"%(name,n) for name in tmp.colnames] )
            base = fill_nan(hstack((base, tmp)))

        # Add in any masked bits
        if len(masked):
            masked.rename_columns(
                masked.colnames, ["%s_0" % n for n in masked.colnames])
            base = fill_nan(vstack((base, masked)))
        return base


    def inner_match(self, base, cat, join_type="or", cartesian=False):
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
        if not len(base): return cat.copy()

        base = fill_nan(base.copy())
        col_names = [n for n in self._col_names if n in cat.colnames]
        cat = fill_nan(cat[col_names].copy())

        if not cartesian:
            sky_coords_1 = self._sky_coords_not_cartesian(base)
            sky_coords_2 = self._sky_coords_not_cartesian(cat)
        else:
            sky_coords_1 = self._sky_coords_cartesian(base)
            sky_coords_2 =  self._sky_coords_cartesian(cat)

        #######################
        # The actual Matching #
        #######################
        idx, d2d, d3d = sky_coords_2.match_to_catalog_3d(sky_coords_1)
        tmp = Table(
            np.full((len(base), len(col_names)), np.nan),
            names=col_names, dtype=cat[col_names].dtype)

        if cartesian:
            dist = d3d
            threshold = self._threshold.value
        else:
            dist = d2d
            threshold = self._threshold

        for src, IDX, sep in zip(cat, idx, dist):
            self._load()
            if self._verbose:
                self._load.show()

            ##GOODMATCH
            if (sep <= threshold) and (sep == min(dist[idx == IDX])):
                tmp[IDX] = src

            ## Append a source
            elif join_type == "or":
                tmp.add_row(src)
        return tmp


    def finish_matching(
        self, tab, error_column="eflux", num_thresh=-1, zp_mag=0,
        col_names=None):
        """
        Averaging all the values. Combining source flags and building a NUM
        column

        :param tab: Table to work on
        :type tab: astropy.table.Table
        :param error_column: Column containing resultant photometric errors
                             ("eflux" or "stdflux")
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
        flags = np.full(len(tab), SRC_GOOD, dtype=np.uint16)
        av = Table(None)

        if col_names is None:
            col_names = self._col_names
        for ii, name in enumerate(col_names):
            if all_cols := find_col_names(tab, name):
                ar = tab2array(tab, col_names=all_cols)
                if ar.shape[1] > 1:
                    if name == "flux":
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                        mean = np.nanmean(ar, axis=1)
                        if "stdflux" not in self._col_names:
                            av.add_column(
                                Column(np.nanstd(ar, axis=1),
                                       name="stdflux"),
                                index=ii + 1)
                        ## if median and mean are >5% different, flag as SRC_VAR
                        flags[np.abs(mean-col)>(col/5.0)] |= SRC_VAR
                    elif name == "eflux":
                        col = Column(np.sqrt(np.nansum(ar * ar, axis=1)),
                                     name=name)
                    elif name == "stdflux":
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                    elif name == "flag":
                        col = Column(flags, name=name)
                        for f_col in ar.T:
                            flags |= f_col.astype(np.uint16)
                    elif name == "NUM":
                        col = Column(np.nansum(ar, axis=1), name=name)
                    elif name == CAT_NUM:
                        col = Column(all_cols[0], name=name)
                    else:
                        col = Column(np.nanmedian(ar, axis=1), name=name)
                else:
                    col = tab[all_cols[0]]
                    col.name = name
                av.add_column(col,index=ii)

        av["flag"] = Column(flags, name="flag")
        if "flux" in av.colnames:
            ecol = av[error_column] if error_column in av.colnames else None
            mag, mag_err = flux2mag(av["flux"], flux_err=ecol)
            mag += zp_mag

            if self._filter in av.colnames:
                av.remove_column(self._filter)
            if "e%s" % self._filter in av.colnames:
                av.remove_column("e%s" % self._filter)
            av.add_column(mag, name=self._filter)
            av.add_column(mag_err, name="e%s" % self._filter)

        if "NUM" not in av.colnames:
            narr = np.nansum(np.invert(
                np.isnan(tab2array(tab, find_col_names(tab, RA)))), axis=1)
            av.add_column(Column(narr, name="NUM"))

            if num_thresh > 0:
                av.remove_rows( av["NUM"] < num_thresh)
        return av