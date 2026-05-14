from typing import override
from starbug2.matching.generic_match import GenericMatch


class DitherMatch(GenericMatch):
    """
    The same as Generic Matching

    """

    def __init__(self, catalogues, p_file=None):
        super().__init__(catalogues, p_file, method="DitherMatch")

    @override
    def match(self, **kwargs):
        return None