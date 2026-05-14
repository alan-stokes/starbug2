from typing import override

from starbug2.constants import CAT_NUM
from starbug2.matching.generic_match import GenericMatch
from starbug2.utils import h_cascade, fill_nan


class CascadeMatch(GenericMatch):
    """
    A simple advancement on "Generic Matching" where the number
    of columns are not preserved. At the end of each sub match, the
    table is left justified, to reduce the total number of columns needed.
    """
    def __init__(self, **kwargs):
        super().__init__(**kwargs, method="Cascade Matching")

    @override
    def match(self, catalogues, **kwargs):
        """
        Match a list of catalogues with RA and DEC columns

        :param catalogues: The input catalogues to work on
        :type catalogues: list (astropy.Tables)
        :param kwargs: Keyword arguments passed to GenericMatch.match
        :return: A left aligned catalogue of all the matched values
        :rtype: astropy.table.Table
        """
        catalogues = self.init_catalogues(catalogues)
        if CAT_NUM in self._col_names:
            self._col_names.remove(CAT_NUM)
        base = self.build_meta(catalogues)

        for n, cat in enumerate(catalogues, 1):
            self._load.msg = "matching: %d" % n
            tmp = self.inner_match(base, cat, join_type="or")
            tmp.rename_columns(
                tmp.colnames, ["%s_%d" % (name, n) for name in tmp.colnames] )
            base = h_cascade([base, tmp], col_names=self._col_names)
        base = fill_nan(base)
        return base