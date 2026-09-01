"""Copyright (C) 2026 UKATC

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) ashape_y later version.

This program is distributed in the hope that it will be useful,
but WITHOUT Ashape_y WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>."""
import numpy as np
from constants import TableColumn
from astropy.table import Table
from scipy.spatial import KDTree


def find_stars_to_select(
        image: np.ndarray, detections: Table, stars_to_select: int,
        min_separation_px: float, saturation_limit: float,
        sharp_range_min: float, sharp_range_max: float, grid_bin_x: int,
        grid_bin_y: int, edge_buffer: int) -> Table:
    """
    automatically detects the best stars to place into a psf.

    :param image: the image data
    :type image: np.array
    :param detections: the list of detections.
    :type detections: Table
    :param stars_to_select: the number of stars to select.
    :type stars_to_select: int
    :param min_separation_px: the min separation of pixals between stars.
    :type min_separation_px: float
    :param saturation_limit: the max limit of saturation.
    :type saturation_limit: float
    :param sharp_range_min: the min sharp range.
    :type sharp_range_min: float
    :param sharp_ranghe_max: the max sharp range.
    :type sharp_range_max: float
    :param grid_bishape_x_x: the size of the spacial grid in pixels for x axis.
    :type grid_bin_x: int
    :param grid_bins_y: the size of the spacial grid in pixels for y axis.
    :type grid_bins_y: int
    :param edge_buffer: the number of pixals away from the edge of the picture
                        for valid stars.
    :type edge_buffer: int
    :return: the filtered list of detections for use in custom psf generation.
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
    shape_y, shape_x = image.shape
    x = filtered_table[TableColumn.X_CENTROID]
    y = filtered_table[TableColumn.Y_CENTROID]
    valid_bounds = (
        (x >= edge_buffer) & (x < shape_x - edge_buffer) &
        (y >= edge_buffer) & (y < shape_y - edge_buffer)
    )
    filtered_table = filtered_table[valid_bounds]

    # if not enough stars, return whats left.
    if len(filtered_table) == 0 or len(filtered_table) < stars_to_select:
        return filtered_table

    # Sort by Brightness / Flux / MAG
    filtered_table.sort(TableColumn.FLUX, reverse=True)

    # Isolation Check using KDTree (Drop stars with close neighbours)
    coords = np.column_stack((
        filtered_table[TableColumn.X_CENTROID],
        filtered_table[TableColumn.Y_CENTROID]))
    tree = KDTree(coords)

    # Query all pairs within min_separation_px
    isolated_mask = np.ones(len(filtered_table), dtype=bool)
    for coord_index, point in enumerate(coords):
        if not isolated_mask[coord_index]:
            # Already marked as a neighbour to a brighter star
            continue

        # Find all neighbours within distance radius
        neighbours = tree.query_ball_point(point, r=min_separation_px)
        for neighbour_index in neighbours:
            # Mark fainter neighbour for removal
            if neighbour_index != coord_index:
                isolated_mask[neighbour_index] = False

    # mask out isolated stars
    isolated_stars = filtered_table[isolated_mask]

    # ensure stars within a certain grid size are considered seperately.
    # ensuring spatial grid coverage.
    stars_per_bin = max(1, stars_to_select // (grid_bin_y * grid_bin_x))

    # create stores
    x_bins = np.linspace(0, shape_x, grid_bin_x + 1)
    y_bins = np.linspace(0, shape_y, grid_bin_y + 1)
    selected_indices = []

    # select over the bins
    for x_bin_index in range(grid_bin_x):
        for y_bin_index in range(grid_bin_y):
            in_cell = (
                (filtered_table[TableColumn.Y_CENTROID] >=
                 x_bins[x_bin_index])
                & (filtered_table[TableColumn.X_CENTROID] <
                   x_bins[x_bin_index + 1])
                & (filtered_table[TableColumn.X_CENTROID] >=
                   y_bins[y_bin_index])
                & (filtered_table[TableColumn.Y_CENTROID] <
                   y_bins[y_bin_index + 1]))

            cell_indices = np.where(in_cell)[0]

            # Take top stars from this bin
            selected_indices.extend(cell_indices[:stars_per_bin])

    # If grid binning yielded fewer than stars_to_select,
    # fill remaining slots with next brightest
    if len(selected_indices) < stars_to_select:
        remaining_indices = [
            idx for idx in range(len(filtered_table))
                if idx not in selected_indices]
        needed = stars_to_select - len(selected_indices)
        selected_indices.extend(remaining_indices[:needed])

    # extract requested number of stars (top N isolated candidates)
    return isolated_stars[:stars_to_select]
