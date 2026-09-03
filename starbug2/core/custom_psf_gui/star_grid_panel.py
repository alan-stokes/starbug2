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
import os
from typing import List, Tuple, cast, Any

import numpy as np
from PyQt6 import QtCore
from PyQt6.QtCore import Qt, QSize
from PyQt6.QtWidgets import (
    QDialog, QGridLayout, QHBoxLayout, QLabel, QPushButton, QScrollArea,
    QVBoxLayout, QWidget, QListWidget, QCheckBox, QGroupBox, QMessageBox,
)
from astropy.io.fits import ImageHDU, Header
from astropy.table import Table, Column
from photutils.psf import ImagePSF
from pyqtgraph import ImageItem, GraphicsLayoutWidget, ViewBox

from constants import ExitStates, TableColumn, FileExtensions
from custom_psf_gui.background_generator import BackgroundGenerator
from custom_psf_gui.common_gui_code import (
    detect_stars, create_gui_instance_with_icon,
    run_starbug_for_image_and_ap_file, STAR_IMAGE_SIZE)
from custom_psf_gui.psf_star_selector import find_stars_to_select
from custom_psf_gui.scale_elements import ScaleElements
from star_bug_config import StarBugMainConfig
from utilities.utils import printf, export_table

# size of an image in pixels when taking multiple into account.
DEFAULT_SQUARE_SIZE = 200


class StarGridPanel(QDialog):
    """Pop-up panel / solo gui containing scale parameters and a grid of
       astronomical images with the ability to do PSF generation is solo."""


    @staticmethod
    def boot_up(config: StarBugMainConfig) -> ExitStates:
        image_data: np.ndarray | None
        detections: Table | None
        exit_state: ExitStates
        if config.ap_file is None:
            image_data, detections, exit_state = detect_stars(config)
            if exit_state != ExitStates.EXIT_SUCCESS:
                return exit_state
        else:
            image_data, detections, exit_state = (
                run_starbug_for_image_and_ap_file(config))
            if exit_state != ExitStates.EXIT_SUCCESS:
                return exit_state

        assert image_data is not None
        assert detections is not None

        # select stars for psf.
        selected_stars, error = find_stars_to_select(
            image_data, detections, config.psf_generator_stars_to_select,
            config.psf_generator_min_separation,
            config.psf_generator_saturation_limit,
            config.sharp_cutoff_low, config.sharp_cutoff_high,
            config.psf_generator_grid_bin_x, config.psf_generator_grid_bin_y,
            config.psf_generator_edge_buffer)


        if error is not None and selected_stars is None:
            printf(f"Automatic selection failed with error: {error}. "
                   f"Shutting down")
            return ExitStates.EXIT_FAIL

        if selected_stars is not None and error is not None:
            if len(selected_stars) < config.psf_min_allowed_stars:
                printf(
                    f"Automatic selection failed, as the amount of stars "
                    f"detected ({len(selected_stars)}) was less than the "
                    f"minimum allowed of {config.psf_min_allowed_stars}. "
                    f"Shutting down. ")
                return ExitStates.EXIT_FAIL
            else:
                printf(
                    f"Automatic selection failed to reach your requested "
                    f"number of stars ({config.psf_generator_stars_to_select}"
                    f") but as it reached the threshold of minimum allowed "
                    f"stars ({config.psf_min_allowed_stars}) "
                    f"We will continue.")

        # populate images from selected stars.
        images: List[Tuple[str, np.ndarray]] = []
        assert selected_stars is not None

        # get coords and iterate over them to build the images to display.
        x_coords: Column = selected_stars[TableColumn.X_CENTROID]
        y_coords: Column = selected_stars[TableColumn.Y_CENTROID]
        cat_id: Column = selected_stars[TableColumn.CAT_NUM]
        x: float
        y: float
        for x, y, cat_id in zip(x_coords, y_coords, cat_id):
            x_center: int = int(round(x))
            y_center: int = int(round(y))
            img_height, img_width = image_data.shape

            x_min: int = max(0, x_center - STAR_IMAGE_SIZE)
            x_max: int = min(img_width, x_center + STAR_IMAGE_SIZE)
            y_min: int = max(0, y_center - STAR_IMAGE_SIZE)
            y_max: int = min(img_height, y_center + STAR_IMAGE_SIZE)

            star_data: np.ndarray = image_data[y_min:y_max, x_min:x_max]
            images.append((cat_id, star_data))

        # this allows us to debug whilst in test mode as well as working
        # via command line.
        # noinspection PyArgumentList
        gui_app, app_icon = create_gui_instance_with_icon()
        custom_psf_gui = StarGridPanel(
            images, True, ScaleElements.DEFAULT_SCALE_LIST_SELECTED,
            ScaleElements.DEFAULT_SCALE_LIST_MUT_SELECTED, image_data,
            config, detections)
        custom_psf_gui.setWindowTitle("starbug2 PSF generator")

        if not app_icon.isNull():
            custom_psf_gui.setWindowIcon(app_icon)

        custom_psf_gui.show()

        try:
            return_code = gui_app.exec()
            if return_code == 0:
                return ExitStates.EXIT_SUCCESS
            else:
                print(f"GUI closed with non-zero exit code: {return_code}")
                return ExitStates.EXIT_FAIL
        except Exception as e:
            print(f" failed to run custom PFS GUI due to: {e}")
            return ExitStates.EXIT_FAIL



    @staticmethod
    def window_clicked(
        event, window: GraphicsLayoutWidget,
        select_cb: QCheckBox) -> None:
        # Ignore right clicks if needed, or toggle on any left
        # click
        if event.button() == Qt.MouseButton.LeftButton:
            # Toggles checked state
            select_cb.toggle()

            # Triggers the callback
            # noinspection PyUnresolvedReferences
            select_cb.clicked.emit(select_cb.isChecked())

        # Call original base event processing
        GraphicsLayoutWidget.mousePressEvent(window, event)

    def __init__(
            self, images: List[Tuple[str, np.ndarray]], sole_ui: bool,
            scale_selected_row: int, scale_selected_mut_row: int,
            image_data: np.ndarray, config: StarBugMainConfig,
            detected_stars: Table, parent=None,
            original_selected_stars: list[str] | None=None):
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
        :param original_selected_stars: the original selected stars, used
               when showing the generated psf.
        :type original_selected_stars: list[str] | None
        """
        super().__init__(parent)
        self.setWindowTitle("Image Inspection & PSF generation")
        self._images: List[Tuple[str, np.ndarray]] = images
        self._original_selected_stars: list[str] | None = (
            original_selected_stars)
        self._image_items: List[ImageItem] = list()
        self._solo: bool = sole_ui
        self._selected_stars: list[str] = list()

        # Track grid structures for dynamic reflowing
        self._cell_widgets: List[QWidget] = []
        self._grid_layout: QGridLayout | None = None
        self._current_cols: int = 3

        # update selected stars to all given
        star_id: str
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

        # selection elements.
        self._de_select_all: QPushButton | None = None
        self._select_all: QPushButton | None = None
        self._check_boxes: list[QCheckBox] = []

        # apply button.
        self._apply_button: QPushButton | None = None

        # save button
        self._save_button: QPushButton | None = None

        # create the UI.
        self._create_components(scale_selected_row, scale_selected_mut_row)

    def _create_checkbox(
            self, window: GraphicsLayoutWidget, star_id: str) -> None:
        """
        creates a checkbox.
        :param window: thee window to apply the checkbox to.
        :type window: GraphicsLayoutWidget
        :param star_id: the id of the star
        :type star_id: str
        :return: None
        """
        # add check box
        the_select_box: QCheckBox = QCheckBox(window)
        the_select_box.setChecked(True)
        self._check_boxes.append(the_select_box)

        # position in top right corner of image.
        cb_size: QSize = the_select_box.sizeHint()
        the_select_box.move(DEFAULT_SQUARE_SIZE - cb_size.width() - 6, 6)

        # ensure clean background.
        the_select_box.setStyleSheet("""
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
        the_select_box.clicked.connect(
            lambda checked, s_id=star_id: self._on_select_cb(s_id, checked))

        # Add event handler to ensure tick matches when clicking on image.
        window.mousePressEvent = lambda event: StarGridPanel.window_clicked(
            event, window, the_select_box)

    def _create_cell(
            self, grid_widget: QWidget, star_id: str,
            img_data: np.ndarray) -> QWidget:
        """
        creates a cell for the grid.

        :param grid_widget: the grid holder.
        :return: the cell
        :rtype: QWidget
        """
        cell_widget: QWidget = QWidget(grid_widget)
        cell_layout: QVBoxLayout = QVBoxLayout(cell_widget)
        cell_layout.setContentsMargins(0, 0, 0, 0)
        cell_layout.setSpacing(0)

        # Create a PyQtGraph GraphicsLayoutWidget for each cell
        window: GraphicsLayoutWidget = GraphicsLayoutWidget()
        window.setFixedSize(DEFAULT_SQUARE_SIZE, DEFAULT_SQUARE_SIZE)

        # add image
        view: ViewBox = window.addViewBox()
        view.setAspectLocked(True)

        # disable zoom and pan
        view.setMouseEnabled(x=False, y=False)

        # Disables the right-click context menu
        view.setMenuEnabled(False)
        view.enableAutoRange(axis=ViewBox.XYAxes, enable=True)

        # add image data to image / window
        img_item: ImageItem = ImageItem(img_data.T)
        view.addItem(img_item)
        self._image_items.append(img_item)

        if len(self._selected_stars) != 1:
            self._create_checkbox(window, star_id)

        # Assemble cell layout
        cell_layout.addWidget(window)
        return cell_widget


    def _create_image_viewer(self) -> QScrollArea:
        """
        builds the image viewer.
        :return: The image viewer area.
        :rtype: QScrollArea
        """
        scroll_area: QScrollArea = QScrollArea(self)
        scroll_area.setWidgetResizable(True)

        # create grid wigit
        grid_widget: QWidget = QWidget(self)
        grid_layout: QGridLayout = QGridLayout(grid_widget)
        grid_layout.setAlignment(
            QtCore.Qt.AlignmentFlag.AlignTop |
            QtCore.Qt.AlignmentFlag.AlignLeft)

        # Populate grid (e.g., 3 columns)
        cols: int = 3
        for idx, (star_id, img_data) in enumerate(self._images):
            row = idx // cols
            col = idx % cols
            cell_widget: QWidget = self._create_cell(
                grid_widget, star_id, img_data)
            grid_layout.addWidget(cell_widget, row, col)

        # ensures the squares are not stretched to rectangles
        grid_layout.setRowStretch(grid_layout.rowCount(), 1)
        scroll_area.setWidget(grid_widget)
        return scroll_area

    def _update_selects(self, new_state: bool) -> None:
        """
        updates all the checkboxes based off the new state.
        :param new_state: the new state of the checkboxes
        :type new_state: bool
        :return: None
        """
        check_box: QCheckBox
        for check_box in self._check_boxes:
            check_box.setChecked(new_state)

        # ensure the internal state works.
        self._selected_stars.clear()
        if new_state:
            star_id: str
            for star_id, _ in self._images:
                self._selected_stars.append(star_id)


    def _build_select_buttons(self, controls_layout: QVBoxLayout) -> None:
        """
        builds the select all and unselect all buttons.
        :param controls_layout: the layout to store the buttons in
        :type controls_layout: QVBoxLayout
        :return: None
        """
        self._de_select_all: QPushButton = QPushButton("Remove all", self)
        # noinspection PyUnresolvedReferences
        self._de_select_all.clicked.connect(lambda: self._update_selects(False))

        self._select_all: QPushButton = QPushButton("Select All", self)
        # noinspection PyUnresolvedReferences
        self._select_all.clicked.connect(lambda: self._update_selects(True))

        # add to widget
        controls_layout.addWidget(self._select_all)
        controls_layout.addWidget(self._de_select_all)

    def _build_scale_components(
            self, controls_layout: QVBoxLayout, scale_selected_row: int,
            scale_selected_mut_row: int) -> None:
        """
        builds the scale components.
        :param controls_layout: the control layout to instill them.
        :type controls_layout: QVBoxLayout
        :param scale_selected_row: the selected row number
        :type scale_selected_row: int
        :param scale_selected_mut_row: the selected mut row number
        :type scale_selected_mut_row: int
        :return: None
        """
        # extract just the image data.
        image_data: list[np.ndarray] = []
        img_data: np.ndarray
        for _, img_data in self._images:
            image_data.append(img_data)

        # add scaling component
        self._scale_builder: ScaleElements = ScaleElements(
            image_data, self._image_items)
        assert self._scale_builder is not None
        scaling_group: QGroupBox
        scaling_group, self._scaling_list, self._scaling_list_mut = (
            self._scale_builder.create_scaling_group(
                self, scale_selected_row, scale_selected_mut_row))
        controls_layout.addWidget(scaling_group)
        # trigger scaling
        self._scale_builder.on_scaling_item_clicked()

    def _build_trigger_buttons(self, controls_layout: QVBoxLayout) -> None:
        # Example action button inside the GUI if in solo mode.
        if self._solo:
            self._apply_button: QPushButton = QPushButton(
                "Execute PSF generation")
            controls_layout.addWidget(self._apply_button, stretch=1)
            # noinspection PyUnresolvedReferences
            self._apply_button.clicked.connect(self._solo_execute_psf_generation)
            self._update_apply_enable_state()
        else:
            self._apply_button = QPushButton("update PSF star selection")
            controls_layout.addWidget(self._apply_button, stretch=1)
            # noinspection PyUnresolvedReferences
            self._apply_button.clicked.connect(self._solo_update_psf_selection)

    def _update_apply_enable_state(self) -> None:
        """
        updates the apply button enable state based off x selected.
        :return: None
        """
        assert self._apply_button is not None
        if not self._solo:
            return

        if len(self._selected_stars) > self._config.psf_min_allowed_stars:
            self._apply_button.setEnabled(True)
        else:
            self._apply_button.setEnabled(False)

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
        controls_layout: QVBoxLayout = QVBoxLayout()
        controls_layout.addWidget(QLabel("Control Panel"))

        self._build_scale_components(
            controls_layout, scale_selected_row, scale_selected_mut_row)

        # add select buttons
        if len(self._selected_stars) != 1:
            self._build_select_buttons(controls_layout)

            # build trigger buttons
            self._build_trigger_buttons(controls_layout)
        else:
            self._build_save_button(controls_layout)

        # push everything to the top
        controls_layout.addStretch()
        return controls_layout

    def _build_save_button(self, controls_layout: QVBoxLayout) -> None:
        """
        creates the save button.
        :param controls_layout: the layout to instill them.
        :type controls_layout: QVBoxLayout
        :return: None
        """
        self._save_button: QPushButton = QPushButton("Save")
        controls_layout.addWidget(self._save_button, stretch=1)
        # noinspection PyUnresolvedReferences
        self._save_button.clicked.connect(self._save_state)

    def _save_state(self) -> None:
        """
        saves the selected stars as a fits file for users to read as well as
        the psf file.
        :return: None
        """
        # filter detections to just selected stars
        catalog_ids: np.ndarray = np.asarray(
            self._detected_stars[TableColumn.CAT_NUM], dtype=str
        )
        selected_ids: np.ndarray = np.asarray(
            self._original_selected_stars, dtype=str)
        mask: np.ndarray = np.isin(catalog_ids, selected_ids)
        filtered_table: Table = self._detected_stars[mask]

        # write to files
        output_dir: str | None = self._config.output_file
        assert output_dir is not None
        file_name: str = os.path.join(
            str(self._config.output_file),
            f"psf_selected_stars{FileExtensions.CUSTOM_PSF}")
        export_table(filtered_table, file_name, header=None)
        file_name_psf: str = os.path.join(
            str(self._config.output_file),
            f"custom{FileExtensions.CUSTOM_PSF}")
        ImageHDU(data=cast(Any, self._images[0][1]), header=Header()).writeto(
            file_name_psf, overwrite=True)

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
        main_layout: QHBoxLayout = QHBoxLayout(self)
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
        controls_widget: QWidget = QWidget(self)
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

    def _handle_failure(self, error: str) -> None:
        """
        handles failure of epsf.
        :param error: the error message.
        :return: None
        """
        assert self._apply_button is not None
        self._apply_button.setEnabled(True)
        msg_box = QMessageBox(self)
        msg_box.setIcon(QMessageBox.Icon.Critical)
        msg_box.setWindowTitle("Error")
        msg_box.setText(f"Failed to generate PSF:\n{error}")
        msg_box.exec()

    def _handle_success(self, epsf: ImagePSF) -> None:
        """
        handles a successful psf generation.
        :return: None
        """
        assert self._apply_button is not None
        self._apply_button.setEnabled(True)

        assert self._scaling_list is not None
        assert self._scaling_list_mut is not None

        dialog = StarGridPanel(
            parent=self, images=[("Star_PSF", epsf.data)], sole_ui=False,
            scale_selected_row=self._scaling_list.currentRow(),
            config=self._config,
            scale_selected_mut_row=self._scaling_list_mut.currentRow(),
            detected_stars=self._detected_stars, image_data=self._image_data,
            original_selected_stars=self._selected_stars)
        dialog.exec()

    def _solo_execute_psf_generation(self) -> None:
        """
        execute the PSF generation.
        :return: None
        """
        assert self._selected_stars is not None
        assert self._scaling_list is not None
        assert self._scaling_list_mut is not None
        assert self._apply_button is not None
        self._apply_button.setEnabled(False)
        self._background_generator = BackgroundGenerator()
        assert self._background_generator is not None
        self._background_generator.generate_epsf_and_view(
            self._selected_stars, self._detected_stars,
            self._image_data.copy(), self._config, self,
            self._handle_failure, self._handle_success
        )

    def _solo_update_psf_selection(self) -> None:
        """
        update the selection. and hands back to the initial GUI.
        :return: None
        """
        self.accept()

    @property
    def selected_stars(self) -> List[str]:
        """
        hands over the selected stars to parent if needed.
        NOTE: used with the main GUI.
        :return: the selected stars
        :type: List[str]
        """
        return self._selected_stars
