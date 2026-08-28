from typing import Tuple

import numpy as np
from PyQt6.QtWidgets import (
    QGroupBox, QListWidget, QFormLayout, QMainWindow, QDialog)
from pyqtgraph import ImageItem
from astropy.visualization import (
    AsinhStretch,
    HistEqStretch,
    ImageNormalize,
    LinearStretch,
    LogStretch,
    MinMaxInterval,
    PowerStretch,
    SinhStretch,
    SqrtStretch,
    SquaredStretch,
    ZScaleInterval,
)


class ScaleElements:

    # Mapping UI labels to Astropy stretch objects
    STRETCH_MAP = {
        "Linear": LinearStretch(),
        "Log": LogStretch(a=1000.0),
        "Power": PowerStretch(a=2.0),
        "Sqrt": SqrtStretch(),
        "Squared": SquaredStretch(),
        "AsinH": AsinhStretch(a=0.1),
        "SinH": SinhStretch(a=0.5),
        "Histogram": None
    }

    # Mapping UI labels to Astropy interval objects
    INTERVAL_MAP = {
        "Min max": MinMaxInterval(),
        "Z scale": ZScaleInterval(),
    }

    def __init__(
            self, images: list[np.ndarray],
            image_items: list[ImageItem]) -> None:
        """
        constructor
        :param images: the list of images
        :param image_items: the list of image items.
        """
        self._image_data = images
        self._img_items = image_items
        self._scaling_list: QListWidget | None = None
        self._scaling_list_mut: QListWidget | None = None

    @staticmethod
    def scale_astronomy_image(
            image_data: np.ndarray, stretch_name: str,
            interval_name: str) -> np.ndarray:
        """Scales a 2D numpy image array using Astropy visualisation based on
           UI selections.

        :param image_data: Raw 2D numpy array of the astronomy image.
        :param stretch_name: Selected item from scaling_list (
                             e.g., 'Linear', 'AsinH').
        :param interval_name: Selected item from scaling_list_mut (
                              'Min max' or 'Z scale').
        :return: 2D numpy array normalised to [0.0, 1.0] ready for
                 QImage/display.
        """
        # 1. Select Interval (Min max vs Z scale)
        if interval_name == "Z scale":
            interval = ZScaleInterval()
        else:
            interval = MinMaxInterval()

        # 2. Select Stretch
        if stretch_name == "Log":
            stretch = LogStretch(a=1000.0)
        elif stretch_name == "Power":
            stretch = PowerStretch(a=2.0)
        elif stretch_name == "Sqrt":
            stretch = SqrtStretch()
        elif stretch_name == "Squared":
            stretch = SquaredStretch()
        elif stretch_name == "AsinH":
            stretch = AsinhStretch(a=0.1)
        elif stretch_name == "SinH":
            stretch = SinhStretch(a=0.5)
        elif stretch_name == "Histogram":
            stretch = HistEqStretch(image_data)
        else:
            stretch = LinearStretch()

        # 3. Create Normalisation wrapper
        norm = ImageNormalize(
            image_data, interval=interval, stretch=stretch, clip=True
        )

        # Return normalised floats [0.0, 1.0]
        return norm(image_data)

    @staticmethod
    def adjust_list_widget_height(list_widget: QListWidget) -> None:
        """
        Resizes the QListWidget to fit its content items exactly.
        :param list_widget: the wigit to resize.
        :return: None
        """
        # Force layout recalculation to ensure valid row heights
        list_widget.doItemsLayout()

        # Sum total height of all items
        total_height = 0
        for i in range(list_widget.count()):
            total_height += list_widget.sizeHintForRow(i)

        # Add frame margin padding (top/bottom borders)
        total_height += list_widget.frameWidth() * 2

        # Set maximum height so it shrinks to content
        list_widget.setMaximumHeight(total_height)

    def create_scaling_group(
            self, parent: QMainWindow | QDialog, default_selection_list: int,
            default_selection_mul: int) -> Tuple[
            QGroupBox, QListWidget, QListWidget]:
        """
        Creates the scaling group
        :param parent: the parent
        :type parent: QMainWindow
        :param default_selection_list: which item to select by default.
        :type default_selection_list: int
        :param default_selection_mul: which item to select by default.
        :type default_selection_mul: int
        :return: the scaling group
        """
        scale_param_group = QGroupBox("Scales Parameters", parent)

        # support collapsing this detection param group.
        scale_param_group.setCheckable(True)
        scale_param_group.setChecked(True)

        # Toggling the title checkbox hides/shows the form parameters inside
        # noinspection PyUnresolvedReferences
        scale_param_group.toggled.connect(
            lambda checked: [
                param_form.itemAt(i).widget().setVisible(checked)
                for i in range(param_form.count())
                if param_form.itemAt(i).widget()
            ]
        )

        self._scaling_list = QListWidget(parent)
        assert self._scaling_list is not None
        self._scaling_list.addItems([
            "Linear", "Log", "Power", "Sqrt", "Squared", "AsinH", "SinH",
            "Histogram"])
        self._scaling_list.setCurrentRow(default_selection_list)
        # noinspection PyUnresolvedReferences
        self._scaling_list.itemClicked.connect(self.on_scaling_item_clicked)
        self.adjust_list_widget_height(self._scaling_list)

        self._scaling_list_mut = QListWidget(parent)
        assert self._scaling_list_mut is not None
        self._scaling_list_mut.addItems(["Min max", "Z scale"])
        self._scaling_list_mut.setCurrentRow(default_selection_mul)
        # noinspection PyUnresolvedReferences
        self._scaling_list_mut.itemClicked.connect(
            self.on_scaling_item_clicked)
        self.adjust_list_widget_height(self._scaling_list_mut)

        # fix max heights
        self._scaling_list.setMaximumHeight(
            self._scaling_list.sizeHintForRow(0) *
            self._scaling_list.count() + 10)
        self._scaling_list_mut.setMaximumHeight(
            self._scaling_list_mut.sizeHintForRow(0) *
            self._scaling_list_mut.count() + 10)
        self._scaling_list.setMinimumHeight(
            self._scaling_list.sizeHintForRow(0) *
            self._scaling_list.count() + 10)
        self._scaling_list_mut.setMinimumHeight(
            self._scaling_list_mut.sizeHintForRow(0) *
            self._scaling_list_mut.count() + 10)

        # add to the form
        param_form = QFormLayout(scale_param_group)
        param_form.addRow(self._scaling_list)
        param_form.addRow(self._scaling_list_mut)
        return scale_param_group, self._scaling_list, self._scaling_list_mut

    def on_scaling_item_clicked(self) -> None:
        """
        execute a scaling adjustment.
        :return: None
        """
        assert self._scaling_list is not None
        assert self._scaling_list_mut is not None
        selected_stretch_items = self._scaling_list.selectedItems()
        selected_interval_items = self._scaling_list_mut.selectedItems()

        if not selected_stretch_items or not selected_interval_items:
            return

        stretch_choice = selected_stretch_items[0].text()
        interval_choice = selected_interval_items[0].text()

        # Process image data using Astropy
        for image_data, img_item in zip(self._image_data, self._img_items):
            display_data = ScaleElements.scale_astronomy_image(
                image_data.copy(), stretch_choice, interval_choice
            )

            # Pass display_data (scaled 0-1) to your Matplotlib canvas or
            # QImage renderer
            img_item.setImage(
                display_data.T, autoLevels=False, levels=(0.0, 1.0))
