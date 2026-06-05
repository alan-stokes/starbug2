from __future__ import annotations
from typing import List, Optional, Tuple
import getopt
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import Polygon
from astropy.table import Table

from starbug2.utils import tab2array, colour_index, fill_nan

class Mask(object):

    @staticmethod
    def from_file(f_name: str) -> Mask:
        """
        makes a mask object from a file. The file must have only 1 line
        in it, which must match the format defined in Mask.from_string.

        :param f_name: the file name for the mask string.
        :return: A Mask instance.
        :rtype: Mask
        """
        with open(f_name) as fp:
            return Mask.from_string(fp.readline())


    @staticmethod
    def from_string(string: str) -> Mask:
        # noinspection SpellCheckingInspection
        """
        method to create a mask object from a string. the string should be
        in the following format:

        [-x XCOL] [-y YCOL] [-l Label] : x1 y1 x2 y2 x3 y3 ...

        :param string: the string to create the mask of.
        :return: the constructed mask.
        """
        label: Optional[str] = None
        keys: List[Optional[str]] = [None, None]
        colour: str = 'k'

        # type definitions
        opts_str: str
        coords: str
        opts: List[Tuple[str, str]]
        args: List[str]

        opts_str, coords = string.split(':')
        opts, args = getopt.getopt(opts_str.split(' '), "c:l:x:y:")
        for opt, opt_arg in opts:
            if opt == "-x":
                keys[0] = opt_arg
            if opt == "-y":
                keys[1] = opt_arg
            if opt == "-l":
                label = opt_arg.replace('_', ' ')
            if opt == "-c":
                colour = opt_arg
        strip_coords: List[str] = coords.strip().rstrip().split(' ')
        points: np.ndarray = (
            np.array(strip_coords, dtype=float).reshape(
                (int(len(strip_coords) / 2), 2)))
        return Mask(points, keys, label=label, colour=colour)


    def __init__(self, bounds, keys, label=None, colour="k") -> None:
        """
        mask constructor

        :param bounds: The path vertices, as an array, masked array or
                       sequence of pairs.
        :param keys: iterable array with 2 elements.
        :param label: the mask label
        :param colour: the colour of the mask.
        """

        self._path: Path = Path(bounds)
        if len(keys) == 2:
            self._keys: List[str] = keys
        else:
            raise Exception
        self._label: Optional[str] = label

        self._colour: str = colour


    def apply(self, data_table) -> np.ndarray:
        """
        applies a data table based off the masks keys.
        :param data_table: the table to apply the mask on.
        :return: length-N bool array
        :rtype: ndarray
        """
        d: Table = fill_nan(colour_index(data_table, self._keys))
        return self._path.contains_points(tab2array(d))


    def plot(self, plot_axis, **kwargs) -> None:
        """
        plots a polygon onto the axis.
        :param plot_axis: the axis to plot the polygon onto.
        :param kwargs: arbitrary polygon parameters.
        :return: None
        """
        patch: Polygon = Polygon(
            self._path.vertices,
            label=self._label.replace('_', ' ') if self._label else None,
            fill=False, edgecolor=self._colour, **kwargs)
        plot_axis.add_patch(patch)
