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
from PyQt6.QtCore import QObject, pyqtSignal
from astropy.table import Table
from photutils.psf import EPSFBuildResult, ImagePSF

from main_components.custom_psf import CustomPSF
from star_bug_config import StarBugMainConfig


class PSFWorker(QObject):
    # Signals to communicate results/errors back to the main GUI thread safely
    # Emits the generated e-PSF result object
    finished = pyqtSignal(object, name="success")

    # Emits error messages if execution fails
    error = pyqtSignal(str, name="error")

    def __init__(
            self, filtered_table: Table, image_data: np.ndarray,
            config: StarBugMainConfig, parent: QObject | None = None):
        """
        builds the psf worker.
        :param filtered_table: the table of stars to use for the psf
        generation.
        :type filtered_table: astropy.table.Table
        :param image_data: the image data
        :type image_data: np.ndarray
        :param config: the main starbug config
        :type config: StarBugMainConfig
        :param parent: the parent widget
        :type parent: QObject
        """
        super().__init__(parent)
        self._image_data: np.ndarray = image_data
        self._config: StarBugMainConfig = config
        self._filtered_table: Table = filtered_table

    def run(self) -> None:
        """
        runs the psf worker, which does the heavy PSF processing in a
        background thread. This should stop the UI from locking up and
        resulting in the os thinking its frozen.

        :return: None
        """
        try:
            # Perform calculation here (
            # avoid touching UI components directly inside run)
            result: EPSFBuildResult = CustomPSF.generate_epsf(
                self._filtered_table, self._image_data, self._config)
            epsf: ImagePSF = result.epsf
            self.finished.emit(epsf)
        except Exception as e:
            self.error.emit(str(e))