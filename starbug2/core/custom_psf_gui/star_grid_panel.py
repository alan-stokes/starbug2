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
    QWidget,
)
import pyqtgraph as pg


class StarGridPanel(QDialog):
    """Pop-up panel containing scale parameters and a grid of astronomical images."""

    def __init__(self, images: List[Tuple[str, np.ndarray]], parent=None):
        super().__init__(parent)
        self.setWindowTitle("Image Inspection & Scaler")
        self.resize(800, 600)
        self._images = images or []

        # Main Layout
        main_layout = QVBoxLayout(self)

        # --- Top Section: Scaling Controls Panel ---
        controls_layout = QHBoxLayout()
        controls_layout.addWidget(
            QLabel("Grid Parameters / Controls Panel Here")
        )

        # Example action button inside the dialog
        apply_btn = QPushButton("Apply to Grid")
        controls_layout.addWidget(apply_btn)
        main_layout.addLayout(controls_layout)

        # --- Bottom Section: Scrollable Grid of Images ---
        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)

        grid_widget = QWidget(self)
        grid_layout = QGridLayout(grid_widget)

        # Populate grid (e.g., 3 columns)
        cols = 3
        for idx, (star_id, img_data) in enumerate(self._images):
            row = idx // cols
            col = idx % cols

            # Create a PyQtGraph GraphicsLayoutWidget for each cell
            win = pg.GraphicsLayoutWidget()
            view = win.addViewBox()
            view.setAspectLocked(True)

            img_item = pg.ImageItem(img_data.T)
            view.addItem(img_item)

            grid_layout.addWidget(win, row, col)

        scroll_area.setWidget(grid_widget)
        main_layout.addWidget(scroll_area)