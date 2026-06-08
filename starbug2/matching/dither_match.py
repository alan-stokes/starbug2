from typing import override, Any
from astropy.table import Table
from starbug2.matching.generic_match import GenericMatch


class DitherMatch(GenericMatch):
    """
    The same as Generic Matching
    """

    def __init__(
            self, catalogues: list[Table], p_file: str | None = None) -> None:
        # Note: Routed using keywords to align correctly with
        # GenericMatch.__init__
        super().__init__(p_file=p_file, method="DitherMatch")
        # Ensure the instance sets up column structures using your tracking
        # list
        self.init_catalogues(catalogues)

    @override
    def match(self, **kwargs: Any) -> Table | None:
        """
        Match pipeline implementation placeholder.
        """
        return None