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
from typing import Tuple

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QLabel,
    QPushButton, QGroupBox, QFormLayout, QDoubleSpinBox, QCheckBox, QSpinBox,
    QMessageBox, QScrollArea, QComboBox, QListWidget, QAbstractItemView,
    QDialog)
from astropy.table import Table
from photutils.psf import ImagePSF
from pyqtgraph import ImageItem, GraphicsLayoutWidget, ViewBox

import numpy as np

from starbug2.core.custom_psf_gui.background_generator import (
    BackgroundGenerator)
from starbug2.core.custom_psf_gui.common_gui_code import (
    detect_stars, update_config, create_gui_instance_with_icon, STAR_IMAGE_SIZE)
from custom_psf_gui.scale_elements import ScaleElements
from starbug2.core.custom_psf_gui.star_grid_panel import StarGridPanel
from starbug2.constants import ExitStates, TableColumn
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.custom_psf_gui.clickable_circle_overlay import (
    ClickableCircleOverlay)
from starbug2.core.custom_psf_gui.psf_star_selector import (
    find_stars_to_select)

# Search box radius in pixels
RADIUS: float = 2.0

# the geometry of the UI.
GEOMETRY: Tuple[int, int, int, int] = (100, 100, 1200, 700)


class CustomPSFGui(QMainWindow):
    """
    main entrance class to the custom PSF generation GUI.
    """

    @staticmethod
    def execute_gui(config: StarBugMainConfig) -> ExitStates:
        """
        generates a GUI via the command line arguments interface.

        :param config: the main config file.
        :type config: StarBugMainConfig
        :return: exit state, successful if the gui runs and exits cleanly,
                 fail otherwise.
        :rtype: ExitStates
        """
        image_data: np.ndarray | None
        detections: Table | None
        exit_state: ExitStates
        image_data, detections, exit_state = detect_stars(config)
        if exit_state != ExitStates.EXIT_SUCCESS:
            return exit_state

        assert image_data is not None
        assert detections is not None

        # this allows us to debug whilst in test mode as well as working
        # via command line.
        # noinspection PyArgumentList
        gui_app, app_icon = create_gui_instance_with_icon()

        custom_psf_gui = CustomPSFGui(image_data, detections, config)
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

    def __init__(self, image_data: np.ndarray, detections: Table,
                 config: StarBugMainConfig):
        """
        builds the GUI.

        :param image_data: the image data.
        :type image_data: np.ndarray
        :param detections: the detected detections.
        :type detections: Table
        :param config: the main config.
        :type config: StarBugMainConfig
        """
        super().__init__(parent=None)

        # Data storage
        self._image_data: np.ndarray = image_data
        self._circles_for_psf_generation: (
            dict[str, ClickableCircleOverlay]) = dict()
        self._selected_stars: list[str] = []
        self._detected_stars: Table = detections
        self._config = config

        # info
        self._info_label: QLabel = QLabel()

        # detection form elements
        self._full_width_half_max_spin: QDoubleSpinBox
        self._sig_sky: QDoubleSpinBox
        self._sig_source: QDoubleSpinBox
        self._sharp_lo: QDoubleSpinBox
        self._sharp_hi: QDoubleSpinBox
        self._round1_hi: QDoubleSpinBox
        self._round2_hi: QDoubleSpinBox
        self._smooth_lo: QDoubleSpinBox
        self._smooth_hi: QDoubleSpinBox
        self._ricker_r: QDoubleSpinBox
        self._do_bkg: QCheckBox
        self._do_convolution: QCheckBox
        self._clean_sources: QCheckBox

        # psf form elements
        self._detected_list: QComboBox
        self._selected_list: QComboBox

        # psf automatic params
        self._stars_to_select: QSpinBox
        self._min_separation: QDoubleSpinBox
        self._saturation_limit: QDoubleSpinBox
        self._grid_bin_x: QSpinBox
        self._grid_bin_y: QSpinBox
        self._star_finder_sharp_min: QDoubleSpinBox
        self._star_finder_sharp_max: QDoubleSpinBox
        self._edge_buffer: QDoubleSpinBox

        # image
        self._img_item: ImageItem | None = None

        # scale form elements
        self._scaling_list: QListWidget | None = None
        self._scaling_list_mut: QListWidget | None = None
        self._scale_builder: ScaleElements | None = None

        # buttons
        self._automatic_psf_star_selection_btn: QPushButton | None = None
        self._do_custom_psf_btn: QPushButton | None = None
        self._redo_detection_btn: QPushButton | None = None

        # background generator to stop garbage collector.
        self._background_generator: BackgroundGenerator | None = None

        # Set up UI components (Side-by-side layout)
        self._set_up_components()

        # Initialise image content and NaN overlays
        self._create_nan_layer()

        # Add circles for selection from detections
        self._populate_star_circles()
        self._populate_star_combos()

        # update to correct visualisation scale.
        assert self._scale_builder is not None
        self._scale_builder.on_scaling_item_clicked()

    def _setup_central_widget(self) -> None:
        """
        builds the central widget.
        :return: None
        """
        self._central_widget: QWidget = QWidget(self)
        self.setCentralWidget(self._central_widget)
        self._main_layout: QHBoxLayout = QHBoxLayout(self._central_widget)
        self._central_widget.setObjectName("main_window_frame")
        self._central_widget.setStyleSheet("""
            QWidget#main_window_frame {
                border: 2px solid #555555;
                border-radius: 6px;
            }
            #main_window_frame > QWidget {
                background-color: transparent;
            }
        """)

        # add some padding
        self._main_layout.setContentsMargins(8, 8, 8, 8)

    def _setup_right_panel(self) -> None:
        # Right Image View Panel
        self._graphics_layout: GraphicsLayoutWidget = GraphicsLayoutWidget()
        self._view_box: ViewBox = self._graphics_layout.addViewBox()
        self._view_box.setAspectLocked(True)

        # Image items attached to ViewBox
        self._img_item: ImageItem = ImageItem()
        self._nan_overlay_item: ImageItem = ImageItem()
        self._view_box.addItem(self._img_item)
        self._view_box.addItem(self._nan_overlay_item)

        # Set zoom limits
        height, width = self._image_data.shape
        self._view_box.setRange(
            xRange=(0, width), yRange=(0, height), padding=0)

    def _set_up_components(self) -> None:
        """
        sets up the GUI components
        :return: None
        """
        self.setWindowTitle("StarbugII - Custom e-PSF Selector")
        self.setGeometry(*GEOMETRY)

        # Main Central Widget & Horizontal Layout
        self._setup_central_widget()
        self._setup_right_panel()

        # Left Control Panel
        self._control_panel: QWidget = self._create_control_panel()

        # add widgets
        self._main_layout.addWidget(self._control_panel, stretch=1)
        self._main_layout.addWidget(self._graphics_layout, stretch=1)

    def _add_detection_form_elements(self, param_form: QFormLayout) -> None:
        """
        adds the detection params to the form.

        :param param_form: the form to add the params to.
        :type param_form: QFormLayout
        :return: None
        """
        self._full_width_half_max_spin = QDoubleSpinBox(self)
        self._full_width_half_max_spin.setRange(0.5, 20.0)
        self._full_width_half_max_spin.setValue(
            self._config.full_width_half_max)
        param_form.addRow("FWHM (pixels):", self._full_width_half_max_spin)

        self._sig_sky = QDoubleSpinBox(self)
        self._sig_sky.setRange(0.5, 20.0)
        self._sig_sky.setValue(self._config.sigma_sky)
        param_form.addRow("sigma sky:", self._sig_sky)

        self._sig_source = QDoubleSpinBox(self)
        self._sig_source.setRange(0.5, 20.0)
        self._sig_source.setValue(self._config.sigma_source)
        param_form.addRow("sigma source:", self._sig_source)

        self._sharp_lo = QDoubleSpinBox(self)
        self._sharp_lo.setRange(0.5, 20.0)
        self._sharp_lo.setValue(self._config.sharp_cutoff_low)
        param_form.addRow("sharp low:", self._sharp_lo)

        self._sharp_hi = QDoubleSpinBox(self)
        self._sharp_hi.setRange(0.5, 20.0)
        self._sharp_hi.setValue(self._config.sharp_cutoff_high)
        param_form.addRow("sharp high:", self._sharp_hi)

        self._round1_hi = QDoubleSpinBox(self)
        self._round1_hi.setRange(0.5, 20.0)
        self._round1_hi.setValue(self._config.round1_cutoff_high)
        param_form.addRow("round 1 high:", self._round1_hi)

        self._round2_hi = QDoubleSpinBox(self)
        self._round2_hi.setRange(0.5, 20.0)
        self._round2_hi.setValue(self._config.round2_cutoff_high)
        param_form.addRow("round 2 high:", self._round2_hi)

        self._smooth_lo = QDoubleSpinBox(self)
        self._smooth_lo.setRange(0.5, 20.0)
        self._smooth_lo.setValue(self._config.smooth_low)
        param_form.addRow("sharp low:", self._smooth_lo)

        self._smooth_hi = QDoubleSpinBox(self)
        self._smooth_hi.setRange(0.5, 20.0)
        self._smooth_hi.setValue(self._config.smooth_high)
        param_form.addRow("sharp high:", self._smooth_hi)

        self._ricker_r = QDoubleSpinBox(self)
        self._ricker_r.setRange(0.5, 20.0)
        self._ricker_r.setValue(self._config.ricker_wavelet_radius)
        param_form.addRow("sharp high:", self._ricker_r)

        # ensure the checkboxes are on same line
        checkbox_container = QWidget(self)
        checkbox_layout = QHBoxLayout(checkbox_container)
        checkbox_layout.setContentsMargins(0, 0, 0, 0)

        self._do_bkg = QCheckBox("background 2d", self)
        self._do_bkg.setChecked(self._config.do_bgd_2d)
        checkbox_layout.addWidget(self._do_bkg)

        self._do_convolution = QCheckBox("do convolution", self)
        self._do_convolution.setChecked(self._config.do_convolution)
        checkbox_layout.addWidget(self._do_convolution)

        self._clean_sources = QCheckBox("clean sources", self)
        self._clean_sources.setChecked(self._config.clean_sources)
        checkbox_layout.addWidget(self._clean_sources)
        param_form.addRow("Options:", checkbox_container)

    def _create_detection_param_group(self) -> QGroupBox:
        """
        creates the detection param group.
        :return: the detection param group
        :rtype: QGroupBox
        """
        detection_param_group = QGroupBox("Detection Parameters", self)

        # support collapsing this detection param group.
        detection_param_group.setCheckable(True)
        detection_param_group.setChecked(True)

        # Toggling the title checkbox hides/shows the form parameters inside
        # noinspection PyUnresolvedReferences
        detection_param_group.toggled.connect(
            lambda checked: [
                param_form.itemAt(i).widget().setVisible(checked)
                for i in range(param_form.count())
                if param_form.itemAt(i).widget()
            ]
        )

        # populate the detection form elements.
        param_form = QFormLayout(detection_param_group)
        self._add_detection_form_elements(param_form)
        return detection_param_group

    def _create_buttons(self) -> None:
        """
        generates the buttons
        :return: tuple of the redo detection button, and the do custom psf
        button
        :rtype: None
        """
        self._redo_detection_btn = QPushButton("Redo Detection", self)
        # noinspection PyUnresolvedReferences
        self._redo_detection_btn.clicked.connect(self.on_redo_detection)

        self._do_custom_psf_btn = QPushButton("Generate Custom PSF", self)
        assert self._do_custom_psf_btn is not None
        self._do_custom_psf_btn.setStyleSheet("font-weight: bold;")
        # noinspection PyUnresolvedReferences
        self._do_custom_psf_btn.clicked.connect(self.on_generate_custom_psf)

    def _create_psf_stars_group(self):
        psf_stars_group = QGroupBox("custom PSF selections", self)

        # support collapsing this detection param group.
        psf_stars_group.setCheckable(True)
        psf_stars_group.setChecked(True)

        def toggle_group_contents(checked: bool) -> None:
            """
            ensures disabled aspects for the psf star panel.
            NOTE. the same work for the detections panel doesn't work, as the
            group box doesn't like it and is more inclined to disable instead
            of going invisible.
            :param checked: if its to be visible or invisible.
            :return: None
            """
            widgets = [
                self._detected_list,
                self._selected_list,
                transfer_stars_to_selected,
                transfer_stars_to_detections,
                params_group,
                review_selected,
                self._automatic_psf_star_selection_btn
            ]
            for w in widgets:
                assert w is not None
                w.setVisible(checked)
                # Overrides Qt's default grey-out behaviour
                w.setEnabled(True)

        # Toggling the title checkbox hides/shows the form parameters inside
        # noinspection PyUnresolvedReferences
        psf_stars_group.toggled.connect(toggle_group_contents)

        # create holder for the star selection
        group_layout = QVBoxLayout(psf_stars_group)
        dropdown_layout = QHBoxLayout()

        # create the multi select drops downs.
        self._detected_list = QListWidget(self)
        self._selected_list = QListWidget(self)

        # put them in multi select mode
        self._detected_list.setSelectionMode(
            QAbstractItemView.SelectionMode.MultiSelection
        )
        self._selected_list.setSelectionMode(
            QAbstractItemView.SelectionMode.MultiSelection
        )

        # create automatic selection process.
        self._automatic_psf_star_selection_btn = QPushButton(
            "Pick PSF Stars", self)
        # noinspection PyUnresolvedReferences
        self._automatic_psf_star_selection_btn.clicked.connect(
            self.on_automatic)

        # add button
        transfer_stars_to_selected = QPushButton("Move to selected ->", self)
        # noinspection PyUnresolvedReferences
        transfer_stars_to_selected.clicked.connect(
            self.on_psf_add_star_selected)

        transfer_stars_to_detections = QPushButton("Move back <-", self)
        # noinspection PyUnresolvedReferences
        transfer_stars_to_detections.clicked.connect(
            self.on_psf_remove_star_selected)

        # review button
        review_selected = QPushButton("Review selection", self)
        # noinspection PyUnresolvedReferences
        review_selected.clicked.connect(self.do_review_stars)

        # add automatic params
        self._stars_to_select = QSpinBox(self)
        self._stars_to_select.setRange(1, 1000)
        self._stars_to_select.setValue(
            self._config.psf_generator_stars_to_select)

        self._min_separation = QDoubleSpinBox(self)
        self._min_separation.setRange(1, 100)
        self._min_separation.setValue(
            self._config.psf_generator_min_separation)

        self._saturation_limit = QDoubleSpinBox(self)
        self._saturation_limit.setRange(1, 1000000000)
        self._saturation_limit.setValue(
            self._config.psf_generator_saturation_limit)

        self._grid_bin_x = QSpinBox(self)
        self._grid_bin_x.setRange(1, 100)
        self._grid_bin_x.setValue(
            self._config.psf_generator_grid_bin_x)

        self._grid_bin_y = QSpinBox(self)
        self._grid_bin_y.setRange(1, 100)
        self._grid_bin_y.setValue(
            self._config.psf_generator_grid_bin_y)

        self._star_finder_sharp_min = QDoubleSpinBox(self)
        self._star_finder_sharp_min.setRange(0, 100)
        self._star_finder_sharp_min.setValue(
            self._config.sharp_cutoff_low)

        self._star_finder_sharp_max = QDoubleSpinBox(self)
        self._star_finder_sharp_max.setRange(0, 100)
        self._star_finder_sharp_max.setValue(
            self._config.sharp_cutoff_high)

        self._edge_buffer = QDoubleSpinBox(self)
        self._edge_buffer.setRange(1, 100)
        self._edge_buffer.setValue(
            self._config.psf_generator_edge_buffer)

        # add widgets in order.
        params_group = QGroupBox("PSF Parameters", self)
        param_form = QFormLayout(params_group)
        param_form.addRow("N stars", self._stars_to_select)
        param_form.addRow("Min separation", self._min_separation)
        param_form.addRow("Saturation Limit", self._saturation_limit)
        param_form.addRow("Spacial grid x", self._grid_bin_x)
        param_form.addRow("Spacial grid y", self._grid_bin_y)
        param_form.addRow("Sharp min", self._star_finder_sharp_min)
        param_form.addRow("Sharp max", self._star_finder_sharp_max)
        param_form.addRow("Edge buffer", self._edge_buffer)
        group_layout.addWidget(params_group)
        group_layout.addWidget(self._automatic_psf_star_selection_btn)
        dropdown_layout.addWidget(self._detected_list)
        dropdown_layout.addWidget(self._selected_list)
        group_layout.addLayout(dropdown_layout)
        group_layout.addWidget(transfer_stars_to_selected)
        group_layout.addWidget(transfer_stars_to_detections)
        group_layout.addWidget(review_selected)

        return psf_stars_group

    def _create_control_panel(self) -> QWidget:
        """
        creates the control panel
        :return: the control panel QWidget
        :rtype QWidget
        """

        scroll_area = QScrollArea(self)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        scroll_area.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )

        control_panel: QWidget = QWidget(self)
        control_layout: QVBoxLayout = QVBoxLayout(control_panel)
        control_panel.setObjectName("control_panel")

        # create detection field.
        detection_param_group = self._create_detection_param_group()

        # create scaling field.
        assert self._img_item is not None
        self._scale_builder: ScaleElements = ScaleElements(
            [self._image_data], [self._img_item])
        assert self._scale_builder is not None
        scaling_group, self._scaling_list, self._scaling_list_mut = (
            self._scale_builder.create_scaling_group(
                self, ScaleElements.DEFAULT_SCALE_LIST_SELECTED,
                ScaleElements.DEFAULT_SCALE_LIST_MUT_SELECTED))
        self._scale_builder.on_scaling_item_clicked()

        # create psf stars field.
        psf_stars_group = self._create_psf_stars_group()

        # Action Buttons
        self._create_buttons()

        # add to the control panel layout.
        control_layout.addWidget(detection_param_group)
        control_layout.addWidget(scaling_group)
        control_layout.addWidget(psf_stars_group)
        control_layout.addWidget(self._redo_detection_btn)
        control_layout.addWidget(self._do_custom_psf_btn)

        # put everything at top
        control_layout.addStretch()

        # add the control panel to the scroll area.
        scroll_area.setWidget(control_panel)

        # Status Label
        self._info_label: QLabel = QLabel(
            "Click a star marker to select.", self)
        self._info_label.setWordWrap(True)
        control_layout.addWidget(self._info_label)

        return scroll_area

    def _create_nan_layer(self) -> None:
        """ creates the nan graphical layer.
        :return: None
        """
        """Prepares background data and loads main/NaN overlays safely."""
        display_data = np.copy(self._image_data)
        nan_mask = np.isnan(display_data)

        # Clean Nans for main image rendering
        bg_median = float(np.nanmedian(display_data))
        display_data[nan_mask] = bg_median

        assert self._img_item is not None
        self._img_item.setImage(display_data.T, autoLevels=True)

        # Red NaN Overlay LUT
        nan_lut = np.array([
            # Valid pixels -> transparent
            [0, 0, 0, 0],
            # Nans -> solid red
            [255, 0, 0, 255]
        ], dtype=np.uint8)

        self._nan_overlay_item.setImage(
            nan_mask.T.astype(np.uint8), levels=[0, 1], lut=nan_lut,
            autoLevels=False)

    def _populate_star_combos(self) -> None:
        """Populates dropdowns based on current detection state."""
        self._detected_list.clear()
        self._selected_list.clear()

        # Populate detected dropdown with stars that are not currently selected
        for star_id in self._circles_for_psf_generation.keys():
            if star_id not in self._selected_stars:
                self._detected_list.addItem(star_id)

        # Populate selected dropdown
        for star_id in self._selected_stars:
            self._selected_list.addItem(star_id)

    def _populate_star_circles(self) -> None:
        """
        populate circles from detections.
        :return: None
        """
        # Clear existing overlays
        for star_circle in self._circles_for_psf_generation.values():
            self._view_box.removeItem(star_circle)
        self._circles_for_psf_generation.clear()

        # Draw interactive circles
        detection_needed: Table = self._detected_stars[
            TableColumn.CAT_NUM, TableColumn.X_CENTROID,
            TableColumn.Y_CENTROID]
        for idx, x, y in detection_needed:
            star_id: str = f"Star_{idx}"
            star_circle = ClickableCircleOverlay(
                pos=(x, y),
                radius=RADIUS,
                star_id=star_id,
                callback=self.on_star_clicked
            )
            self._view_box.addItem(star_circle)
            self._circles_for_psf_generation[star_id] = star_circle

    def _populate_star_lists(self) -> None:
        """
        updates the 2-star lists for the psf selection.
        :return: None
        """
        self._detected_list.clear()
        self._selected_list.clear()

        # Populate detected stars list
        star_id: str
        for star_id in self._circles_for_psf_generation.keys():
            if star_id not in self._selected_stars:
                self._detected_list.addItem(star_id)

        # Populate selected stars list
        for star_id in self._selected_stars:
            self._selected_list.addItem(star_id)

    def on_star_clicked(self, star_id: str, pos: tuple[float, float]) -> None:
        """
        callback when a star circle is selected.

        :param star_id: the star id
        :type star_id: str
        :param pos: the location.
        :type pos: tuple[float, float]
        :return: None
        """
        if star_id in self._selected_stars:
            self._circles_for_psf_generation[star_id].turn_off()
            self._selected_stars.remove(star_id)
            self._info_label.setText(
                f"Unselected: {star_id} at Pixel Coordinates (X={pos[0]:.2f}, "
                f"Y={pos[1]:.2f})")
        else:
            self._circles_for_psf_generation[star_id].turn_on()
            self._selected_stars.append(star_id)
            self._info_label.setText(
                f"Selected: {star_id} at Pixel Coordinates (X={pos[0]:.2f}, "
                f"Y={pos[1]:.2f})")
        self._populate_star_lists()

    def on_redo_detection(self) -> None:
        """
        executes when redoing detection
        :return: None
        """
        assert self._redo_detection_btn is not None
        self._redo_detection_btn.setEnabled(False)

        config_copy = update_config(self._config)
        # add the values from the form.
        config_copy.unfreeze()
        config_copy.full_width_half_max = (
            self._full_width_half_max_spin.value())
        config_copy.sigma_sky = self._sig_sky.value()
        config_copy.sigma_source = self._sig_source.value()
        config_copy.sharp_cutoff_low = self._sharp_lo.value()
        config_copy.sharp_cutoff_high = self._sharp_hi.value()
        config_copy.round1_cutoff_high = self._round1_hi.value()
        config_copy.round2_cutoff_high = self._round2_hi.value()
        config_copy.smooth_low = self._smooth_lo.value()
        config_copy.smooth_high = self._smooth_hi.value()
        config_copy.ricker_wavelet_radius = self._ricker_r.value()
        config_copy.do_bgd_2d = self._do_bkg.isChecked()
        config_copy.do_convolution = self._do_convolution.isChecked()
        config_copy.clean_sources = self._clean_sources.isChecked()
        config_copy.freeze()

        # run and get new detections.
        self._info_label.setText("updating detected stars.")
        _, detections, exit_state = CustomPSFGui._detect_stars(config_copy)

        if exit_state != ExitStates.EXIT_SUCCESS:
            self._info_label.setText("failed to detect stars.")
        else:
            self._info_label.setText("located stars.")

        assert detections is not None
        self._detected_stars: Table = detections
        self._populate_star_circles()
        self._selected_stars.clear()
        self._populate_star_lists()

        self._info_label.setText(
            "Detection redo complete. Click a star marker to select.")
        self._redo_detection_btn.setEnabled(True)

    def _handle_failure(self, error: str) -> None:
        """
        handles failure of epsf.
        :param error: the error message.
        :return: None
        """
        self._info_label.setText("PSF Generation Failed.")
        assert self._do_custom_psf_btn is not None
        self._do_custom_psf_btn.setEnabled(True)
        msg_box = QMessageBox(self)
        msg_box.setIcon(QMessageBox.Icon.Critical)
        msg_box.setWindowTitle("Error")
        msg_box.setText(f"Failed to generate PSF:\n{error}")
        msg_box.exec()

    def _handle_success(self, epsf: ImagePSF) -> None:
        """
        handles a successful psf generation.
        :param epsf: the epsf.
        :type epsf: ImagePSF
        :return: None
        """
        self._info_label.setText("PSF Generation Succeeded.")
        assert self._do_custom_psf_btn is not None
        self._do_custom_psf_btn.setEnabled(True)

        assert self._scaling_list is not None
        assert self._scaling_list_mut is not None
        dialog = StarGridPanel(
            parent=self, images=[("Star_PSF", epsf.data)], sole_ui=False,
            scale_selected_row=self._scaling_list.currentRow(),
            config=self._config,
            scale_selected_mut_row=self._scaling_list_mut.currentRow(),
            detected_stars=self._detected_stars, image_data=self._image_data)
        dialog.exec()

    def on_generate_custom_psf(self) -> None:
        """
        executes when generating the custom psf.
        :return: None
        """
        assert self._do_custom_psf_btn is not None
        self._do_custom_psf_btn.setEnabled(False)
        assert self._selected_stars is not None
        assert self._scaling_list is not None
        assert self._scaling_list_mut is not None
        self._info_label.setText("Generating custom PSF. Please wait...")
        self._background_generator = BackgroundGenerator()
        assert self._background_generator is not None
        self._background_generator.generate_epsf_and_view(
            self._selected_stars, self._detected_stars,
            self._image_data.copy(), self._config, self,
            self._handle_failure, self._handle_success
        )

    def on_psf_add_star_selected(self) -> None:
        """
        moves stars from the detections into the psf selection list.
        :return: None
        """
        selected_items: list = self._detected_list.selectedItems()
        if not selected_items:
            return

        for item in selected_items:
            star_id: str = item.text()
            if star_id not in self._selected_stars:
                self._selected_stars.append(star_id)

                # Highlight circle overlay on canvas
                if star_id in self._circles_for_psf_generation:
                    self._circles_for_psf_generation[star_id].turn_on()

        self._populate_star_lists()

    def on_psf_remove_star_selected(self) -> None:
        """
        removes stars selected from the psf remove and puts them back in the
        detections list.
        :return: None
        """
        selected_items: list = self._selected_list.selectedItems()
        if not selected_items:
            return

        for item in selected_items:
            star_id: str = item.text()
            if star_id in self._selected_stars:
                self._selected_stars.remove(star_id)

                # Turn off highlight circle overlay on canvas
                if star_id in self._circles_for_psf_generation:
                    self._circles_for_psf_generation[star_id].turn_off()

        self._populate_star_lists()

    def on_automatic(self) -> None:
        """
        execute an automatic detection of stars for the psf.
        :return: None
        """
        assert self._automatic_psf_star_selection_btn is not None
        self._automatic_psf_star_selection_btn.setEnabled(False)

        selected_stars_from_algorithm: Table | None
        error: str | None
        selected_stars_from_algorithm, error = find_stars_to_select(
            self._image_data, self._detected_stars,
            self._stars_to_select.value(), self._min_separation.value(),
            self._saturation_limit.value(),
            self._star_finder_sharp_min.value(),
            self._star_finder_sharp_max.value(), self._grid_bin_x.value(),
            self._grid_bin_y.value(), self._edge_buffer.value())

        # handle states
        if error is not None and selected_stars_from_algorithm is None:
            self._info_label.setText(error)
            msg_box = QMessageBox(self)
            msg_box.setIcon(QMessageBox.Icon.Critical)
            msg_box.setWindowTitle("Selection Failed")
            msg_box.setText(
                f"{error}\n\nPlease adjust your filter parameters and try "
                f"again.")
            msg_box.exec()
            self._automatic_psf_star_selection_btn.setEnabled(True)
            return
        if selected_stars_from_algorithm is not None and error is not None:
            self._info_label.setText(error)
            reply = QMessageBox.warning(
                self, "Spatial Grid Selection Warning",
                f"{error}\n\nDo you want to continue with the available "
                f"stars?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                # Default focused button
                QMessageBox.StandardButton.No
            )
            if reply == QMessageBox.StandardButton.No:
                self._automatic_psf_star_selection_btn.setEnabled(True)
                return

        # update selected stars and the ui.
        self._info_label.setText(
            "Populating image with selected Stars. Please wait....")
        self._selected_stars.clear()
        assert selected_stars_from_algorithm is not None
        for star_id in selected_stars_from_algorithm[TableColumn.CAT_NUM]:
            self._selected_stars.append(f"Star_{star_id}")
        self._update_ui_elements()
        self._automatic_psf_star_selection_btn.setEnabled(True)
        self._info_label.setText("Completed Pick PSF Stars.")

    def _update_ui_elements(self):
        # update tables and plot.
        self._populate_star_circles()
        self._populate_star_lists()

        # remove and reapply circles as needed
        for star in self._circles_for_psf_generation.values():
            star.turn_off()
        for star_id in self._selected_stars:
            self._circles_for_psf_generation[star_id].turn_on()

    def do_review_stars(self) -> None:
        """
        opens up the review panel which allows fine tune selection of stars
        :return: None
        """
        selected_stars_arrays: list[Tuple[str, np.ndarray]] = []

        # verify data exists
        if len(self._selected_stars) == 0:
            return

        # extract mini images for each selected star
        for item in self._selected_stars:
            star_id: str = item.split("Star_")[1]
            row_mask = self._detected_stars[TableColumn.CAT_NUM] == star_id
            matched_row = self._detected_stars[row_mask]
            x_coord: float = float(matched_row[TableColumn.X_CENTROID][0])
            y_coord: float = float(matched_row[TableColumn.Y_CENTROID][0])

            x_center: int = int(round(x_coord))
            y_center: int = int(round(y_coord))
            img_height, img_width = self._image_data.shape

            x_min: int = max(0, x_center - STAR_IMAGE_SIZE)
            x_max: int = min(img_width, x_center + STAR_IMAGE_SIZE)
            y_min: int = max(0, y_center - STAR_IMAGE_SIZE)
            y_max: int = min(img_height, y_center + STAR_IMAGE_SIZE)

            star_data: np.ndarray = self._image_data[y_min:y_max, x_min:x_max]
            selected_stars_arrays.append((star_id, star_data))

        assert self._scaling_list is not None
        assert self._scaling_list_mut is not None
        dialog = StarGridPanel(
            parent=self, images=selected_stars_arrays, sole_ui=False,
            scale_selected_row=self._scaling_list.currentRow(),
            scale_selected_mut_row=self._scaling_list_mut.currentRow(),
            image_data=self._image_data, config=self._config,
            detected_stars=self._detected_stars)
        if dialog.exec() == QDialog.DialogCode.Accepted:

            # update selected stars
            updated_stars: list[str] = dialog.selected_stars
            self._selected_stars.clear()
            for star_id in updated_stars:
                self._selected_stars.append(f"Star_{star_id}")

            # update tables and plot.
            self._update_ui_elements()
