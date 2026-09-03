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
from PyQt6.QtCore import QThread
from PyQt6.QtWidgets import QWidget
from astropy.table import Table
from collections.abc import Callable

from photutils.psf import ImagePSF

from constants import TableColumn
from custom_psf_gui.psf_worker import PSFWorker
from star_bug_config import StarBugMainConfig


class BackgroundGenerator:

    def __init__(self):
        self._psf_thread: QThread | None = None
        self._worker: PSFWorker | None = None

    def generate_epsf_and_view(
            self, selected_stars: list[str], detected_stars: Table,
            image_data: np.ndarray, config: StarBugMainConfig, parent: QWidget,
            failure_callback: Callable[[str], None],
            success_callback: Callable[[ImagePSF], None]) -> None:
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
        :param failure_callback: the function to call when a failure occurs.
        :type failure_callback: Callable[[str], None]
        :param success_callback: the callback when it is successful.
        :type success_callback: Callable[None]
        :return: None
        """
        # build the list of stars.
        clean_selected_ids = {
            s.replace("Star_", "") for s in selected_stars}
        table_ids = np.char.strip(
            detected_stars[TableColumn.CAT_NUM].astype(str))
        mask = np.isin(table_ids, list(clean_selected_ids))
        filtered_table = detected_stars[mask].copy()

        # run the custom psf process.
        self._psf_thread = QThread(parent)
        assert self._psf_thread is not None
        self._psf_thread.setObjectName("psf_thread")
        self._worker = PSFWorker(filtered_table, image_data.copy(), config)
        assert self._worker is not None
        self._worker.moveToThread(self._psf_thread)

        # connect the signals.
        # noinspection PyUnresolvedReferences
        self._psf_thread.started.connect(self._worker.run)
        self._worker.finished.connect(lambda epsf: success_callback(epsf))
        self._worker.error.connect(failure_callback)

        # cleanup
        self._worker.finished.connect(self._psf_thread.quit)
        self._worker.finished.connect(self._worker.deleteLater)
        # noinspection PyUnresolvedReferences
        self._psf_thread.finished.connect(self._psf_thread.deleteLater)

        # start the process.
        self._psf_thread.start()
