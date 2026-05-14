"""
Starbug matching functions
Primarily this is the main routines for dither/band/generic matching which are
 at the core of starbug2 and starbug2-match
"""
from typing import override

import numpy as np
import astropy.units as u
from astropy.table import Table, hstack
from starbug2.constants import FILTER, RA, DEC
from starbug2.filters import filters
from starbug2.matching.generic_match import GenericMatch
from starbug2.utils import (
    Loading, printf, p_error, fill_nan, find_col_names, warn, puts)

# keys for catalogue fields.
# noinspection SpellCheckingInspection
_OBS = "OBSERVTN"
_VISIT = "VISIT"
_EXPOSURE = "EXPOSURE"


class BandMatch(GenericMatch):
    # filter flag for kwargs.
    FILTER = "fltr"
    THRESHOLD = "threshold"

    # match methods
    _FIRST = "first"
    _LAST = "last"
    _BOOT_STRAP = "bootstrap"

    # warning messages
    _WRONG_THRESHOLD = (
        "Threshold values must be scalar or list with length 1 less than the "
        "catalogue list. The final element is being ignored.\n")


    def __init__(self, **kwargs):
        if BandMatch.FILTER in kwargs:
            if not isinstance(kwargs[BandMatch.FILTER], list):
                warn("{} input should be a list, "
                     "there may be unexpected behaviour\n", self.FILTER)

        if BandMatch.THRESHOLD in kwargs:
            if isinstance(kwargs[BandMatch.THRESHOLD], list):
                kwargs[BandMatch.THRESHOLD] = (
                    np.array(kwargs[BandMatch.THRESHOLD]))

        super().__init__(**kwargs, method="Band Matching")

    def order_catalogues(self,catalogues):
        """
        Reorder catalogue list into increasing wavelength size.
        This only works for JWST bands. Unrecognised filters will be left
        unchanged. The function should also set the self.FILTER variable if 
        possible.

        :param catalogues: List of astropy tables with meta keys 'FILTER'
        :type catalogues: list[astropy.table.Table]
        :return: The same list reordered by wavelength
        :rtype: list[astropy.table.Table]
        """

        status = -1
        _ii = None
        sorters = [
            ## META in JWST filters
            lambda t: list(filters.keys()).index(t.meta.get(FILTER)),

            ## col_names in JWST filters
            lambda t: list(filters.keys()).index(
                (set(t.col_names) & set(filters.keys())).pop()),

            ## META in self.filters
            lambda t: self.FILTER.index( t.meta.get(FILTER)),

            ## col_names in JWST filters
            lambda t: self.FILTER.index(
                (set(t.col_names) & set(self.FILTER)).pop())
        ]

        for n, fn in enumerate(sorters):
            try:
                catalogues.sort(key=fn)
                _ii = map(fn, catalogues)
                status = n
                break
            except (KeyError, AttributeError, TypeError, ValueError):
                pass

        if status < 0:
            p_error(
                "Unable to reorder catalogues, leaving input order"
                " untouched.\n")
        ## JWST filters
        elif status <= 1 and (_ii is not None):
            self._filter = [list(filters.keys())[i] for i in _ii]

        self._load = Loading(sum(len(c) for c in catalogues[1:]))

        return catalogues

    def jwst_order(self,catalogues):
        pass

    @override
    def match(self, catalogues, method="first", **kwargs):
        """
        Given a list of catalogues, it will reorder them into increasing
        wavelength or to match the fltr= keyword in the initializer.
        The matching then uses the shortest wavelength available position.
        I.e If F115W, F444W, F770W are input, the F115W centroid positions will
        be taken as "correct". If a source is not resolved in this band, the
        next most astrometric ally accurate position is taken, i.e. F444W

        :param catalogues: List of `astropy.table.Table` objects containing
        the meta item "FILTER=XXX"
        :type catalogues: List of `astropy.table.Table`
        :param method: Centroid method
            "first" -   Use the position corresponding to the earliest
                        appearance of the source
            "last"  -   Use the position corresponding to the latest
                        appearance of the source
            "bootsrap"- ..
            "average" - ..
        :type method: str
        :param kwargs:
        :return: Matched catalogue
        :rtype: astropy.Table
        """
        catalogues = self.order_catalogues(catalogues)

        if (isinstance(self._filter, list) and 
                len(self._filter) == len(catalogues)):
            printf("Bands: %s\n"%', '.join(self._filter))
        else:
            printf("Bands: Unknown\n")

        if type(self._threshold.value) in (list,np.ndarray):
            if len(self._threshold) != (len(catalogues) - 1):
                warn(self._WRONG_THRESHOLD)
                self._threshold = self._threshold[:-1]
        else:
            self._threshold = (
                np.full(len(catalogues) - 1, self._threshold) * u.arcsec)

        printf("Thresholds: %s\n"%", ".join(
            ["%g\""%g for g in self._threshold.value]))

        if self._col_names is None:
            self._col_names = [
                RA, DEC, "flag", "NUM",
                *self._filter, *["e%s" % f for f in self._filter]]

        printf("Columns: %s\n"%", ".join(self._col_names))

        if method not in (self._FIRST, self._LAST, self._BOOT_STRAP):
            method = self._FIRST

        #########
        # Begin #
        #########

        base = self.build_meta(catalogues)
        _threshold = self._threshold.copy()
        for n, tab in enumerate(catalogues):
            ## Temporarily recast threshold
            self._threshold = _threshold[n - 1]
            self._load.msg = "%s (%g\")" % (
                self._filter[n], self._threshold.value)
            col_names = [
                name for name in self._col_names if name in tab.colnames]

            tmp = self.inner_match(base, tab, join_type="or")

            col_names.remove(RA)
            col_names.remove(DEC)
            base = fill_nan(hstack((base, tmp[col_names])))
            base.rename_columns(
                col_names, ["%s_%d" % (name, n + 1) for name in col_names])

            if RA not in base.col_names:
                base = fill_nan(hstack((tmp[RA, DEC], base)))
            elif method == self._FIRST:
                _mask = np.logical_and(np.isnan(base[RA]), tmp[RA] != np.nan)
                base[RA][_mask] = tmp[RA][_mask]
                base[DEC][_mask] = tmp[DEC][_mask]
            elif method == self._LAST:
                _mask = ~np.isnan(tmp[RA])
                base[RA][_mask] = tmp[RA][_mask]
                base[DEC][_mask] = tmp[DEC][_mask]
            elif method == self._BOOT_STRAP:
                _mask = ~np.isnan(tmp[RA])
                base.rename_columns((RA, DEC), ("_RA_%d"%n, "_DEC_%d" % n))
                base = hstack((base, tmp[[RA,DEC]]))

        # Set threshold back at the end
        self._threshold= _threshold

        ####################
        # Fix column names #
        for name in self._col_names:
            all_cols = find_col_names(base, name)
            if len(all_cols) == 1:
                base.rename_column(all_cols.pop(), name)

        ################################
        # Finalise NUM and flag column #
        tmp = self.finish_matching(base, col_names=["NUM", "flag"])
        base.remove_columns(
            (*find_col_names(base, "NUM"), *find_col_names(base, "flag")))
        base.add_column(tmp["NUM"], index=2)
        base.add_column(tmp["flag"], index=3)

        return base


    def band_match(self, catalogues, col_names=(RA,DEC)):
        """
        Given a list of catalogues (with filter names in the metadata), match
        them in order of decreasing astrometric accuracy. If F115W, F444W, F770W
        are input, the F115W centroid positions will be taken as "correct". If a
        source is not resolved in this band, the next most astrometric ally
        accurate position is taken, i.e. F444W

        :param catalogues:
        :param col_names:
        :return:
        """

        ### ORDER the tables into the correct order (increasing wavelength)
        tables = np.full(len(filters), None)
        mask = np.full(len(filters), False)
        for tab in catalogues:
            if FILTER in tab.meta.keys():
                if tab.meta[FILTER] in filters:
                    ii = list(filters.keys()).index(tab.meta[FILTER])
                    tables[ii] = tab
                    mask[ii] = True
                else:
                    p_error(
                        "Unknown filter '%s' (skipping)..\n" % tab.meta[FILTER])
            elif _tmp := set(filters.keys()) & set(tab.col_names):
                ii = list(filters.keys()).index(_tmp.pop())
                tables[ii] = tab
                mask[ii] = True
            else:
                p_error("Cannot find 'FILTER' in table meta (skipping)..\n")

        # document bands
        s = "Bands: "
        for filter_string, tab in zip(filters.keys(),tables):
            if tab: s += "%5s "% filter_string
        puts(s)

        ### Match in increasing wavelength order
        base = Table(None)
        load = Loading(
            sum([len(t) for t in tables[mask][1:]]), "matching", res=100)
        for filter_string, tab in zip(filters.keys(), tables):
            if not tab:
                continue

            ## removing empty magnitude rows
            tab.remove_rows(np.isnan(tab[filter_string]))
            load.msg = "matching:%s" % filter_string
            _col_names = (
                list(name for name in tab.col_names if name in col_names))
            if not len(base):
                tmp = tab[_col_names].copy()
            else:
                idx, d2d, _ = self.inner_match(base, tab)
                tmp = Table(
                    np.full( (len(base),len(_col_names)), np.nan),
                    names = _col_names)

                ###################################
                # Hard coding separations for now #
                separation = 0.06
                f_id = list(filters.keys()).index(filter_string)
                if f_id >= list(filters.keys()).index("F277W"):
                    separation = 0.10
                if f_id >= list(filters.keys()).index("F560W"):
                    separation = 0.15
                if f_id >= list(filters.keys()).index("F1000W"):
                    separation = 0.20
                if f_id >= list(filters.keys()).index("F1500W"):
                    separation = 0.25

                for ii, (src, IDX, sep) in enumerate(zip(tab, idx, d2d)):
                    load.msg = (
                        "matching:%s(%.2g\")" % (filter_string, separation))
                    load()
                    load.show()

                    if ((sep <= separation * u.arcsec)
                        and (sep == min(d2d[idx == IDX]))):
                        for name in _col_names: tmp[IDX][name] = src[name]
                    else:
                        tmp.add_row(src[_col_names])

            tmp.rename_column("flag", "flag_%s" % filter_string)
            base = hstack((
                base, tmp[[filter_string, "e%s" % filter_string,
                           "flag_%s" % filter_string]]
            ))
            base = Table(
                base, dtype=[float] * len(base.col_names)).filled(np.nan)

            ### Only keep the most astronomically correct position
            if RA not in base.col_names:
                base = hstack((tmp[[RA, DEC]], base))
            else:
                _mask = np.logical_and( np.isnan(base[RA]), tmp[RA] != np.nan)
                base[RA][_mask] = tmp[RA][_mask]
                base[DEC][_mask] = tmp[DEC][_mask]

        ## Sort out flags
        flag = np.zeros(len(base),dtype=np.uint16)
        for f_col in find_col_names(base, "flag"):
            flag |= base[f_col].value.astype(np.uint16)
            base.remove_column(f_col)
        base.add_column(flag,name="flag")

        return base.filled(np.nan)