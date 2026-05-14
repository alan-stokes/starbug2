"""
Starbug matching functions
Primarily this is the main routines for dither/band/generic matching which are
 at the core of starbug2 and starbug2-match
"""
from typing import override

import numpy as np
from astropy.table import Table, hstack, vstack

from starbug2.constants import CAT_NUM
from starbug2.matching.generic_match import GenericMatch
from starbug2.utils import p_error, fill_nan

# keys for catalogue fields.
# noinspection SpellCheckingInspection
_OBS = "OBSERVTN"
_VISIT = "VISIT"
_EXPOSURE = "EXPOSURE"

class ExactValueMatch(GenericMatch):
    """
    Match catalogues based on a single value in one of their columns
    The normal use case is matching sources based on their *Catalogue_Number*
    value.
    """

    def __init__(self, value=CAT_NUM, **kwargs):
        """
        setup method.

        :param value: Column name to take exact values from
        :type value: str
        :param kwargs:
        """
        self.value = value
        super().__init__(**kwargs, method="Exact Value Matching")

        if "colnames" in kwargs:
            p_error("Colnames not implemented in %s\n" % self.method)

    def __str__(self):
        s=[ "%s:" % self.method,
            "Value: \"%s\"" % self.value,
            "Colnames: %s" % self._col_names,
            ]

        return "\n".join(s)

    @override
    def inner_match(self, base, cat, join_type="or", cartesian=False):
        """
         The low level matching function.

        :param base: Table onto which to match *cat*
        :type base: astropy.Table
        :param cat: the Table to match to *base*
        :type cat:  astropy.Table
        :param join_type: Joining method ("or" to include all sources, "and"
                          for intersections)
        :type join_type: str
        :param cartesian: Whether to use Cartesian coordinates for matching
        :type cartesian: bool
        :return: A new catalogue, it is a reordered version of *cat*, in the
            correct sorting to be h-stacked with *base*
        :rtype:  astropy.Table
        """

        tmp = Table(
            np.full((len(base), len(cat.col_names)), np.nan),
            names=cat.col_names, dtype=cat.dtype, masked=True)
        for col in tmp.columns.values():
            col.mask |= True

        if not len(base):
            return vstack([tmp,cat])

        for src in cat:
            if self._verbose:
                self._load()
                self._load.show()
            ii = np.where(base[self.value] == src[self.value])[0]
            if len(ii):
                tmp[ii] = src
            else:
                tmp.add_row(src)
        return tmp

    @override
    def match(self, catalogues, **kwargs):
        """
        Core matching function

        :param catalogues: List of astropy Tables to match together
        :type catalogues: list of astropy.Table
        :param kwargs:
        :return: Full matched catalogue
        :rtype: astropy.Table
        """

        catalogues = self.init_catalogues(catalogues)
        base = self.build_meta(catalogues)

        if self.value not in self._col_names:
            p_error("Exact value '%s' not in column names.\n" % self.value)
            return None

        for n, cat in enumerate(catalogues, 1):
            self._load.msg = "matching: %d"%n
            tmp = self.inner_match(base, cat)
            tmp.rename_columns(
                tmp.col_names, ["%s_%d" % (name, n) for name in tmp.col_names])
            base = hstack([base,tmp])

            if n > 1:
                ii = (base[self.value].mask
                      & ~base["%s_%d" % (self.value, n)].mask)
                base["%s" % self.value][ii] = (
                    base["%s_%d" % (self.value, n)][ii])
                base.remove_column("%s_%d" % (self.value, n))

            else: base.rename_column("%s_1" % self.value, self.value)

        return fill_nan(base)