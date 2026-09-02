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
from PyQt6 import QtCore
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QVBoxLayout,
    QWidget, QListWidget, QCheckBox,
)
import pyqtgraph as pg
from astropy.table import Table
from pyqtgraph import ImageItem

from custom_psf_gui import common_code
from custom_psf_gui.scale_elements import ScaleElements
from star_bug_config import StarBugMainConfig


class StarGridPanel(QDialog):
    """Pop-up panel / solo gui containing scale parameters and a grid of
       astronomical images with the ability to do PSF generation is solo."""

    def __init__(
            self, images: List[Tuple[str, np.ndarray]], sole_ui: bool,
            scale_selected_row: int, scale_selected_mut_row: int,
            image_data: np.ndarray, config: StarBugMainConfig,
            detected_stars: Table, parent=None):
        """

        :param images: the images needed to present for finer selection
        :type images: List[Tuple[str, np.ndarray]]
        :param sole_ui: if this is a main ui or a dialogue box.
        :type sole_ui: bool
        :param parent: the parent ui.
        :type parent: QWidget
        :param scale_selected_row: the selected row number
        :type scale_selected_row: int
        :param scale_selected_mut_row: the selected mut row number
        :type scale_selected_mut_row: int
        :param image_data: the image data
        :type image_data: np.ndarray
        :param config: the config object
        :type config: StarBugMainConfig
        :param detected_stars: the detected stars
        :type detected_stars: Table
        """
        super().__init__(parent)
        self.setWindowTitle("Image Inspection & PSF generation")
        self._images: List[Tuple[str, np.ndarray]] = images
        self._image_items: List[ImageItem] = list()
        self._solo: bool = sole_ui
        self._selected_stars: list[str] = list()

        # update selected stars to all given
        for star_id, _ in self._images:
            self._selected_stars.append(star_id)

        # scale bits
        self._scaling_list: QListWidget | None = None
        self._scaling_list_mut: QListWidget | None = None
        self._scale_builder: ScaleElements | None = None

        # bits for psf.
        self._image_data: np.ndarray = image_data
        self._config: StarBugMainConfig = config
        self._detected_stars: Table = detected_stars

        # create the UI.
        self._create_components(scale_selected_row, scale_selected_mut_row)

    def _create_image_viewer(self) -> QScrollArea:
        """
        builds the image viewer.
        :return: The image viewer area.
        :rtype: QScrollArea
        """
        scroll_area: QScrollArea = QScrollArea(self)
        scroll_area.setWidgetResizable(True)

        grid_widget: QWidget = QWidget(self)
        grid_layout: QGridLayout = QGridLayout(grid_widget)
        grid_layout.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignTop |
            QtCore.Qt.AlignmentFlag.AlignLeft)

        # Populate grid (e.g., 3 columns)
        cols: int = 3
        square_size: int = 200
        for idx, (star_id, img_data) in enumerate(self._images):
            row = idx // cols
            col = idx % cols

            cell_widget = QWidget(grid_widget)
            cell_layout = QVBoxLayout(cell_widget)
            cell_layout.setContentsMargins(0, 0, 0, 0)
            cell_layout.setSpacing(0)

            # Create a PyQtGraph GraphicsLayoutWidget for each cell
            win = pg.GraphicsLayoutWidget()
            win.setFixedSize(square_size, square_size)

            # add image
            view = win.addViewBox()
            view.setAspectLocked(True)

            # disable zoom and pan
            view.setMouseEnabled(x=False, y=False)
            # Disables the right-click context menu
            view.setMenuEnabled(False)
            view.enableAutoRange(axis=pg.ViewBox.XYAxes, enable=True)

            img_item: ImageItem = pg.ImageItem(img_data.T)
            view.addItem(img_item)
            self._image_items.append(img_item)

            # add check box
            select_cb = QCheckBox(win)
            select_cb.setChecked(True)

            # position in top right corner of image.
            cb_size = select_cb.sizeHint()
            select_cb.move(square_size - cb_size.width() - 6, 6)

            # ensure clean background.
            select_cb.setStyleSheet("""
                QCheckBox {
                    background-color: #ffffff;
                    border: 1px solid #cccccc;
                    border-radius: 3px;
                    padding: 2px;
                }
                QCheckBox::indicator {
                    width: 14px;
                    height: 14px;
                }
            """)

            # noinspection PyUnresolvedReferences
            select_cb.clicked.connect(
                lambda checked, s_id=star_id:
                    self._on_select_cb(s_id, checked))

            # 2. Add click handler to the image canvas
            def make_mouse_press_handler(cb=select_cb):
                def mouse_press_event(event):
                    # Ignore right clicks if needed, or toggle on any left
                    # click
                    if event.button() == Qt.MouseButton.LeftButton:
                        # Toggles checked state
                        cb.toggle()
                        # Triggers the callback
                        cb.clicked.emit(cb.isChecked())
                    # Call original base event processing
                    pg.GraphicsLayoutWidget.mousePressEvent(win, event)
                return mouse_press_event

            # add event handler to ensure tiick matches when clciking on image.
            win.mousePressEvent = make_mouse_press_handler()

            # Assemble cell layout
            cell_layout.addWidget(win)

            grid_layout.addWidget(cell_widget, row, col)

        # ensures the squares are not stretched to rectangles
        grid_layout.setRowStretch(grid_layout.rowCount(), 1)

        scroll_area.setWidget(grid_widget)
        return scroll_area

    def _create_controls_layout(
            self, scale_selected_row: int,
            scale_selected_mut_row: int) -> QVBoxLayout:
        """
        creates the control layout.
        :param scale_selected_row: the selected row number
        :type scale_selected_row: int
        :param scale_selected_mut_row: the selected mut row number
        :type scale_selected_mut_row: int
        :return: the control setup.
        :rtype: QHBoxLayout
        """

        # --- Top Section: Scaling Controls Panel ---
        controls_layout = QVBoxLayout()
        controls_layout.addWidget(QLabel("Control Panel"))

        # extract just the image data.
        image_data = []
        for _, img_data in self._images:
            image_data.append(img_data)

        # add scaling component
        self._scale_builder: ScaleElements = ScaleElements(
            image_data, self._image_items)
        assert self._scale_builder is not None
        scaling_group, self._scaling_list, self._scaling_list_mut = (
            self._scale_builder.create_scaling_group(
                self, scale_selected_row, scale_selected_mut_row))
        controls_layout.addWidget(scaling_group)

        # Example action button inside the GUI if in solo mode.
        if self._solo:
            apply_btn = QPushButton("Execute PSF generation")
            controls_layout.addWidget(apply_btn, stretch=1)
            # noinspection PyUnresolvedReferences
            apply_btn.clicked.connect(self._solo_execute_psf_generation)
        else:
            apply_btn = QPushButton("update PSF star selection")
            controls_layout.addWidget(apply_btn, stretch=1)
            # noinspection PyUnresolvedReferences
            apply_btn.clicked.connect(self._solo_update_psf_selection)

        # push everything to the top
        controls_layout.addStretch()

        # trigger scaling
        self._scale_builder.on_scaling_item_clicked()
        return controls_layout

    def _create_components(
            self, scale_selected_row: int,
            scale_selected_mut_row: int) -> None:
        """
        builds the UI components.
        :param scale_selected_row: the selected row number
        :type scale_selected_row: int
        :param scale_selected_mut_row: the selected mut row number
        :type scale_selected_mut_row: int
        :return: None
        """

        # Main Layout
        main_layout = QHBoxLayout(self)
        self.setStyleSheet("""
            QWidget#main_window_frame {
                border: 2px solid #555555;
                border-radius: 6px;
            }
            #main_window_frame > QWidget {
                background-color: transparent;
            }
        """)

        # handle the main views.
        controls_widget = QWidget(self)
        controls_widget.setFixedWidth(280)

        # create image grid
        image_grid_area: QScrollArea = self._create_image_viewer()

        # create control panel
        controls_layout: QVBoxLayout = self._create_controls_layout(
            scale_selected_row, scale_selected_mut_row)

        controls_widget.setLayout(controls_layout)
        main_layout.addWidget(controls_widget, stretch=0)
        main_layout.addWidget(image_grid_area, stretch=1)

    def _on_select_cb(self, star_id: str, checked: bool) -> None:
        """
        updates the selection for a given star.
        :param star_id: the star id to update.
        :param checked: the value to update to.
        :return: None
        """
        print(f"Star {star_id} selection state: {checked}")
        if checked:
            self._selected_stars.append(star_id)
        else:
            self._selected_stars.remove(star_id)

    def _solo_execute_psf_generation(self) -> None:
        """
        execute the PSF generation.
        :return: None
        """
        assert self._selected_stars is not None
        assert self._scaling_list is not None
        assert self._scaling_list_mut is not None
        common_code.generate_epsf_and_view(
            self._selected_stars, self._detected_stars,
            self._image_data.copy(), self._config, self,
            self._scaling_list.currentRow(),
            self._scaling_list_mut.currentRow()
        )

    def _solo_update_psf_selection(self) -> None:
        """
        update the selection. and hands back to the initial GUI.
        :return: None
        """
        self.accept()

    @property
    def selected_stars(self) -> List[str]:
        return self._selected_stars
