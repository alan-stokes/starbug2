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
import numpy as np
from PyQt6.QtWidgets import QWidget
from astropy.table import Table
from photutils.psf import EPSFBuildResult
from scipy.spatial import KDTree

from constants import TableColumn
from custom_psf_gui.star_grid_panel import StarGridPanel
from main_components.custom_psf import CustomPSF
from star_bug_config import StarBugMainConfig

# the list of scale options
SCALE_LIST: list[str] = [
    "Linear", "Log", "Power", "Sqrt", "Squared", "AsinH", "SinH",
    "Histogram"]
DEFAULT_SCALE_LIST_SELECTED = SCALE_LIST.index("AsinH")

# list of scale options mut.
SCALE_MUT_LIST: list[str] = ["Min max", "Z scale"]
DEFAULT_SCALE_LIST_MUT_SELECTED = SCALE_MUT_LIST.index("Z scale")


def _find_star_matches(
    psf_stars: Table, filtered_stars: Table) -> Table:
    """
    find stars based off closest x adn y coordinates.

    :param psf_stars: the stars found by the psf finder.
    :param filtered_stars: the stars selected by the user.
    :return: The matched stars and a bool stating if all were matched
    :rtype: Table
    """
    filtered_coords: np.ndarray = np.column_stack((
        filtered_stars[TableColumn.X_CENTROID],
        filtered_stars[TableColumn.Y_CENTROID]
    ))
    psf_coords: np.ndarray = np.column_stack((
        psf_stars[TableColumn.X_CENTROID],
        psf_stars[TableColumn.Y_CENTROID]
    ))

    tree = KDTree(psf_coords)
    max_distance_pixels = 2.0
    distances, matching_indices = tree.query(
        filtered_coords, distance_upper_bound=max_distance_pixels)

    valid_matches_mask = distances < max_distance_pixels
    matched_dao_indices = matching_indices[valid_matches_mask]

    matched_dao_sources: Table = psf_stars[matched_dao_indices].copy()
    return matched_dao_sources


def generate_epsf_and_view(
        selected_stars: list[str], detected_stars: Table,
        image_data: np.ndarray, config: StarBugMainConfig, parent: QWidget,
        scaling_list_row: int,  scaling_list_mut_row:int) -> None:
    """
    generate epsf and views result.

    :param selected_stars: the selected stars from the GUI.
    :type selected_stars: list[str]
    :param detected_stars: the detected stars from starbug
    :type detected_stars: Table
    :param image_data: the image data
    :type image_data: np.ndarray
    :param config:the main config for starbug
    :type config: StarBugMainConfig
    :param parent: the parent gui component
    :type parent: QWidget
    :param scaling_list_row: the row of the scaling list to select
    :type scaling_list_row: int
    :param scaling_list_mut_row: the row of the scaling list mut to select
    :type scaling_list_mut_row: int
    :return: None
    """
    # build the list of stars.
    clean_selected_ids = {
        s.replace("Star_", "") for s in selected_stars}
    table_ids = np.char.strip(
        detected_stars[TableColumn.CAT_NUM].astype(str))
    mask = np.isin(table_ids, list(clean_selected_ids))
    filtered_table = detected_stars[mask].copy()

    # get equivalent psf generator table.
    psf_sources: Table | None = CustomPSF.get_psf_sources(
        image_data.copy(), config,
        config.full_width_half_max)

    assert psf_sources is not None
    safe_sources: Table = _find_star_matches(
        psf_sources, filtered_table)
    safe_sources = filtered_table

    # run the custom psf process.
    result: EPSFBuildResult = CustomPSF.generate_epsf(
        safe_sources, image_data.copy(), config)
    epsf = result.epsf

    # open the grid panel with the epsf so the user can inspect.
    dialog = StarGridPanel(
        parent=parent, images=[("Star_PSF", epsf.data)], sole_ui=False,
        scale_selected_row=scaling_list_row, config=config,
        scale_selected_mut_row=scaling_list_mut_row,
        detected_stars=detected_stars, image_data=image_data)
    dialog.exec()
