"""Copyright (C) 2026 UKATC

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT any WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>."""
from typing import Tuple
import numpy as np
from constants import TableColumn
from astropy.table import Table
from scipy.spatial import KDTree


def find_stars_to_select(
        image: np.ndarray, detections: Table, stars_to_select: int,
        min_separation_px: float, saturation_limit: float,
        sharp_range_min: float, sharp_range_max: float, grid_bin_x: int,
        grid_bin_y: int, edge_buffer: int) -> Tuple[Table | None, str | None]:
    """
    automatically detects the best stars to place into a psf.

    :param image: the image data
    :type image: np.array
    :param detections: the list of detections.
    :type detections: Table
    :param stars_to_select: the number of stars to select.
    :type stars_to_select: int
    :param min_separation_px: the min separation of pixels between stars.
    :type min_separation_px: float
    :param saturation_limit: the max limit of saturation.
    :type saturation_limit: float
    :param sharp_range_min: the min sharp range.
    :type sharp_range_min: float
    :param sharp_range_max: the max sharp range.
    :type sharp_range_max: float
    :param grid_bin_x: the size of the spacial grid in pixels for x-axis.
    :type grid_bin_x: int
    :param grid_bin_y: the size of the spacial grid in pixels for y-axis.
    :type grid_bin_y: int
    :param edge_buffer: the number of pixels away from the edge of the picture
                        for valid stars.
    :type edge_buffer: int
    :return: the filtered list of detections for use in custom psf generation,
             or an error string if not enough stars. So a Table if correct
             with None for str. or None for table and a str for error.
             or table and error when spacial returns less, opening the option
             to review them anyhow
    :rtype: tuple[Table | None, str | None]
    """
    filtered_table: Table = detections.copy()

    # handle saturation by dropping stars above that limit
    if saturation_limit is not None:
        filtered_table = filtered_table[
            filtered_table[TableColumn.PEAK] < saturation_limit]

    # remove stars outside sharpness level.
    filtered_table = filtered_table[
        (filtered_table[TableColumn.SHARPNESS] >= sharp_range_min) &
        (filtered_table[TableColumn.SHARPNESS] <= sharp_range_max)]

    # Remove stars too close to array borders.
    shape_y: int
    shape_x: int
    shape_y, shape_x = image.shape
    x: np.ndarray = filtered_table[TableColumn.X_CENTROID]
    y: np.ndarray = filtered_table[TableColumn.Y_CENTROID]
    valid_bounds: np.ndarray = (
        (x >= edge_buffer) & (x < shape_x - edge_buffer) &
        (y >= edge_buffer) & (y < shape_y - edge_buffer)
    )
    filtered_table = filtered_table[valid_bounds]

    # if not enough stars, return what's left.
    if len(filtered_table) == 0 or len(filtered_table) < stars_to_select:
        return None, "Not enough stars meet the selection criteria."

    # Sort by Brightness / Flux / MAG
    filtered_table.sort(TableColumn.FLUX, reverse=True)

    # Isolation Check using KDTree (Drop stars with close neighbours)
    coords: np.ndarray = np.column_stack((
        filtered_table[TableColumn.X_CENTROID],
        filtered_table[TableColumn.Y_CENTROID]))
    tree: KDTree = KDTree(coords)

    # Query all pairs within min_separation_px
    isolated_mask: np.ndarray = np.ones(len(filtered_table), dtype=bool)
    coord_index: int
    point: np.ndarray
    for coord_index, point in enumerate(coords):
        if not isolated_mask[coord_index]:
            # Already marked as a neighbour to a brighter star
            continue

        # Find all neighbours within distance radius
        neighbours: list[int] = tree.query_ball_point(
            point, r=min_separation_px)
        neighbour_index: int
        for neighbour_index in neighbours:
            # Mark fainter neighbour for removal
            if neighbour_index != coord_index:
                isolated_mask[neighbour_index] = False

    # mask out isolated stars
    isolated_stars: Table = filtered_table[isolated_mask]

    if len(isolated_stars) == 0:
        return None, "No isolated stars remain after separation filtering."

    # ensure stars within a certain grid size are considered separately.
    # ensuring spatial grid coverage.
    stars_per_bin: int = max(1, stars_to_select // (grid_bin_y * grid_bin_x))

    # create stores
    x_bins: np.ndarray = np.linspace(0, shape_x, grid_bin_x + 1)
    y_bins: np.ndarray = np.linspace(0, shape_y, grid_bin_y + 1)
    selected_indices: list[int] = []
    iso_x: Table = isolated_stars[TableColumn.X_CENTROID]
    iso_y: Table = isolated_stars[TableColumn.Y_CENTROID]

    # select over the bins
    x_bin_index: int
    for x_bin_index in range(grid_bin_x):
        for y_bin_index in range(grid_bin_y):
            in_cell: np.ndarray = (
                (iso_x >= x_bins[x_bin_index])
                & (iso_x < x_bins[x_bin_index + 1])
                & (iso_y >= y_bins[y_bin_index])
                & (iso_y < y_bins[y_bin_index + 1]))
            cell_indices: np.ndarray = np.where(in_cell)[0]

            # Take top stars from this bin
            star_index: int
            for star_index in cell_indices[:stars_per_bin]:
                if star_index not in selected_indices:
                    selected_indices.append(int(star_index))

    # If grid binning yielded fewer than stars_to_select return error.
    found_count: int = len(selected_indices)
    if found_count < stars_to_select:
        return (
            isolated_stars[selected_indices[:stars_to_select]],
            (f"Only {found_count} optimal spatial PSF stars were available "
             f"(requested {stars_to_select})."))

    # extract requested number of stars (top N isolated candidates)
    final_selected: Table = isolated_stars[selected_indices[:stars_to_select]]
    return final_selected, None
