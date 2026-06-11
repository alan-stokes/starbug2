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