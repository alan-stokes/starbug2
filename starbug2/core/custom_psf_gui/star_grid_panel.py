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
from typing import List, Tuple

import numpy as np
from PyQt6.QtWidgets import (
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget, QListWidget,
)
import pyqtgraph as pg
from pyqtgraph import ImageItem

from custom_psf_gui.scale_elements import ScaleElements


class StarGridPanel(QDialog):
    """Pop-up panel / solo gui containing scale parameters and a grid of
       astronomical images with the ability to do PSF generation is solo."""

    def __init__(
            self, images: List[Tuple[str, np.ndarray]], sole_ui: bool,
            parent=None):
        """

        :param images:
        :param sole_ui:
        :param parent:
        """
        super().__init__(parent)
        self.setWindowTitle("Image Inspection & PSF generation")
        self.resize(800, 600)
        self._images: List[Tuple[str, np.ndarray]] = images
        self._image_items: List[ImageItem] = list()
        self._solo: bool = sole_ui

        # scale bits
        self._scaling_list: QListWidget | None = None
        self._scaling_list_mut: QListWidget | None = None
        self._scale_builder: ScaleElements | None = None

        # create the UI.
        self._create_components()

    def _create_image_viewer(self) -> QScrollArea:
        scroll_area: QScrollArea = QScrollArea(self)
        scroll_area.setWidgetResizable(True)

        grid_widget: QWidget = QWidget(self)
        grid_layout: QGridLayout = QGridLayout(grid_widget)

        # Populate grid (e.g., 3 columns)
        cols: int = 3
        for idx, (star_id, img_data) in enumerate(self._images):
            row = idx // cols
            col = idx % cols

            # Create a PyQtGraph GraphicsLayoutWidget for each cell
            win = pg.GraphicsLayoutWidget()
            view = win.addViewBox()
            view.setAspectLocked(True)

            img_item: ImageItem = pg.ImageItem(img_data.T)
            view.addItem(img_item)
            self._image_items.append(img_item)

            grid_layout.addWidget(win, row, col)

        scroll_area.setWidget(grid_widget)
        return scroll_area

    def _create_controls_layout(self) -> QHBoxLayout:
        # --- Top Section: Scaling Controls Panel ---
        controls_layout = QHBoxLayout()
        controls_layout.addWidget(QLabel("Control Panel"))

        # extract just the image data.
        image_data = []
        for star_id, img_data in self._images:
            image_data.append(img_data)

        # add scaling component
        self._scale_builder: ScaleElements = ScaleElements(
            image_data, self._image_items)
        assert self._scale_builder is not None
        scaling_group, self._scaling_list, self._scaling_list_mut = (
            self._scale_builder.create_scaling_group(self))
        controls_layout.addWidget(scaling_group)

        # Example action button inside the GUI if in solo mode.
        if self._solo:
            apply_btn = QPushButton("Execute PSF generation")
            controls_layout.addWidget(apply_btn)
        else:
            apply_btn = QPushButton("update PSF star selection")
            controls_layout.addWidget(apply_btn)

        return controls_layout

    def _create_components(self) -> None:
        """
        builds the UI components.
        :return: None
        """

        # Main Layout
        main_layout = QVBoxLayout(self)

        # --- Bottom Section: Scrollable Grid of Images ---
        scroll_area: QScrollArea = self._create_image_viewer()
        controls_layout: QHBoxLayout = self._create_controls_layout()

        main_layout.addLayout(controls_layout)
        main_layout.addWidget(scroll_area)