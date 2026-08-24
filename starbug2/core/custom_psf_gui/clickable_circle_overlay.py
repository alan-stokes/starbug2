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
import pyqtgraph as pg


class ClickableCircleOverlay(pg.ROI):
    """Custom interactive circle overlay that highlights on click and reports
       its location."""

    def __init__(self, pos, radius, star_id, callback):
        # Centre the ROI bounding box around the star position
        pos_offset = [pos[0] - radius, pos[1] - radius]
        super().__init__(
            pos=pos_offset, size=[radius * 2, radius * 2], movable=False,
            removable=False)

        self._star_id = star_id
        self._center_pos = pos
        self._callback = callback

        # Make the ROI shape circular visually
        self._aspectLocked = True

        # Style default appearance
        self.turn_off()

    def mouseClickEvent(self, ev) -> None:
        """
        event handler for mouse click on circle.
        :param ev: the event.
        :return: None
        """
        if ev.button() == pg.QtCore.Qt.MouseButton.LeftButton:
            self._callback(self._star_id, self._center_pos)
            ev.accept()
        else:
            super().mouseClickEvent(ev)

    def turn_off(self) -> None:
        """
        turns the circle to not selected mode
        :return: None
        """
        self.setPen(pg.mkPen(color='white', width=1.5))

    def turn_on(self) -> None:
        """
        turns the circle to selected mode.
        :return: None
        """
        self.setPen(pg.mkPen(color='yellow', width=2.5))
