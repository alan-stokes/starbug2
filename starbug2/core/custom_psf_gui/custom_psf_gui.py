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
import copy
import sys
from typing import Tuple

import numpy as np
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QLabel,
    QPushButton, QGroupBox, QFormLayout, QDoubleSpinBox)
from astropy.table import Table
from pyqtgraph import ImageItem, GraphicsLayoutWidget, ViewBox

from starbug2.constants import ExitStates, TableColumn
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug2.core.starbug_main import StarbugBase
from starbug2.core.custom_psf_gui.clickable_circle_overlay import (
    ClickableCircleOverlay)

# Search box radius in pixels
RADIUS = 2.0

# the geometry of the UI.
GEOMETRY = (100, 100, 1200, 700)


class CustomPSFGui(QMainWindow):
    """
    main entrance class to the custom PSF generation GUI.
    """

    @staticmethod
    def _update_config(config: StarBugMainConfig) -> StarBugMainConfig:
        """
        updates a config to just do star detection and aperture photometry.
        :param config: the main config.
        :type config: StarBugMainConfig
        :return: a modified config.
        :rtype: StarBugMainConfig
        """
        config_copy: StarBugMainConfig = copy.deepcopy(config)
        config_copy.unfreeze()
        config_copy.do_star_detection = True
        config_copy.do_aperture_photometry = True
        config_copy.do_photometry_routine = False
        config_copy.do_bgd_estimate = False

        # turn off any other functionality.
        config_copy.execute_jwst_initialisation = False
        config_copy.update_param = False
        config_copy.generate_local_param_file = False
        config_copy.generate_region = False
        config_copy.do_artificial_star_test = False
        config_copy.generate_psf = False
        config_copy.do_custom_psf_gui = False
        config_copy.freeze()
        return config_copy

    @staticmethod
    def _run_starbug_for_detection(config: StarBugMainConfig) ->  (
            Tuple[np.ndarray | None, Table | None, ExitStates]):
        """
        runs basic starbug to generate the image and detections.
        :param config: the config
        :type config: StarBugMainConfig
        :return: a tuple containing the main image, and the detections, and
                 exit state. If not successful, main image and detections will
                 be None.
        :rtype: Tuple[np.ndarray | None, Table | None, ExitStates]
        """
        # create new base and execute
        starbug_base: StarbugBase = StarbugBase(
            config.fits_images[0], config, ap_file=None,
            bkg_file=None)
        result: ExitStates = starbug_base.run_starbug(config)

        # if not successful, pass upwards.
        if result != ExitStates.EXIT_SUCCESS:
            return None, None, result
        else:
            return (starbug_base.main_image().data, starbug_base.detections,
                    result)

    @staticmethod
    def _detect_stars(
            config: StarBugMainConfig) -> (
                Tuple[np.ndarray | None, Table | None, ExitStates]):
        """
        runs basic starbug to generate the image and detections.
        :param config: the main config which will contain detection params.
        :type config: StarBugMainConfig
        :return: a tuple containing the main image, and the detections, and
                 exit state. If not successful, main image and detections will
                 be None.
        :rtype: Tuple[np.ndarray | None, Table | None, ExitStates]
        """
        config_copy = CustomPSFGui._update_config(config)
        return CustomPSFGui._run_starbug_for_detection(config_copy)

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
        image_data, detections, exit_state = CustomPSFGui._detect_stars(config)
        if exit_state != ExitStates.EXIT_SUCCESS:
            return exit_state

        assert image_data is not None
        assert detections is not None

        # this allows us to debug whilst in test mode as well as working
        # via command line.
        # noinspection PyArgumentList
        gui_app = QApplication.instance() or QApplication(sys.argv)
        custom_psf_gui = CustomPSFGui(image_data, detections, config)
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

        # Set up UI components (Side-by-side layout)
        self._set_up_components()

        # Initialise image content and NaN overlays
        self._create_nan_layer()

        # Add circles for selection from detections
        self._populate_star_circles()

    def _set_up_components(self) -> None:
        """
        sets up the GUI components
        :return: None
        """
        self.setWindowTitle("StarbugII - Custom ePSF Selector")
        self.setGeometry(*GEOMETRY)

        # Main Central Widget & Horizontal Layout
        self._central_widget: QWidget = QWidget(self)
        self.setCentralWidget(self._central_widget)
        self._main_layout: QHBoxLayout = QHBoxLayout(self._central_widget)

        # Left Control Panel (Stretch = 1)
        self._control_panel: QWidget = self._create_control_panel()
        self._main_layout.addWidget(self._control_panel, stretch=1)

        # Right Image View Panel (Stretch = 1)
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
        #self._view_box.setLimits(
        #    xMin=0, xMax=width, yMin=0, yMax=height, minXRange=5, minYRange=5
        #)
        self._view_box.setRange(
            xRange=(0, width), yRange=(0, height), padding=0)

        self._main_layout.addWidget(self._graphics_layout, stretch=1)

    def _create_control_panel(self) -> QWidget:
        """
        creates the control panel
        :return: the control panel QWidget
        :rtype QWidget
        """

        control_panel: QWidget = QWidget(self)
        control_layout: QVBoxLayout = QVBoxLayout(control_panel)

        # Status Label
        self._info_label: QLabel = QLabel(
            "Click a star marker to select.", self)
        self._info_label.setWordWrap(True)
        control_layout.addWidget(self._info_label)

        # Group Box: Detection Parameters
        param_group = QGroupBox("Detection Parameters", self)
        param_form = QFormLayout(param_group)

        full_width_half_max_spin = QDoubleSpinBox(self)
        full_width_half_max_spin.setRange(0.5, 20.0)
        full_width_half_max_spin.setValue(self._config.full_width_half_max)
        param_form.addRow("FWHM (pixels):", full_width_half_max_spin)

        sig_sky = QDoubleSpinBox(self)
        sig_sky.setRange(0.5, 20.0)
        sig_sky.setValue(self._config.sigma_sky)
        param_form.addRow("sigma sky:", sig_sky)

        sig_source = QDoubleSpinBox(self)
        sig_source.setRange(0.5, 20.0)
        sig_source.setValue(self._config.sigma_source)
        param_form.addRow("sigma source:", sig_source)

        sharp_lo = QDoubleSpinBox(self)
        sharp_lo.setRange(0.5, 20.0)
        sharp_lo.setValue(self._config.sharp_cutoff_low)
        param_form.addRow("sharp low:", sharp_lo)

        control_layout.addWidget(param_group)

        # Action Buttons
        redo_detection_btn = QPushButton("Redo Detection", self)
        # noinspection PyUnresolvedReferences
        redo_detection_btn.clicked.connect(self.on_redo_detection)
        control_layout.addWidget(redo_detection_btn)

        do_custom_psf_btn = QPushButton("Generate Custom PSF", self)
        do_custom_psf_btn.setStyleSheet("font-weight: bold;")
        # noinspection PyUnresolvedReferences
        do_custom_psf_btn.clicked.connect(self.on_generate_custom_psf)
        control_layout.addWidget(do_custom_psf_btn)

        control_layout.addStretch()
        return control_panel

    def _create_nan_layer(self) -> None:
        """ creates the nan graphical layer.
        :return: None
        """
        """Prepares background data and loads main/NaN overlays safely."""
        display_data = np.copy(self._image_data)
        nan_mask = np.isnan(display_data)

        # Clean NaNs for main image rendering
        bg_median = float(np.nanmedian(display_data))
        display_data[nan_mask] = bg_median

        self._img_item.setImage(display_data.T, autoLevels=True)

        # Red NaN Overlay LUT
        nan_lut = np.array([
            # Valid pixels -> transparent
            [0, 0, 0, 0],
            # NaNs -> solid red
            [255, 0, 0, 255]
        ], dtype=np.uint8)

        self._nan_overlay_item.setImage(
            nan_mask.T.astype(np.uint8), levels=[0, 1], lut=nan_lut,
            autoLevels=False)

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

    def on_redo_detection(self) -> None:
        """
        executes when redoing detection
        :return: None
        """
        config_copy = CustomPSFGui._update_config(self._config)
        # add the values from the form.
        config_copy.unfreeze()

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

        self._info_label.setText(
            "Detection redo complete. Click a star marker to select.")


    def on_generate_custom_psf(self) -> None:
        """
        executes when generating the custom psf.
        :return: None
        """
        # build the list of stars.

        # run the custom psf process.

        # open new viewer for the epsfs.
        pass