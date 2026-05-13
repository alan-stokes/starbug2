import getopt
import numpy as np
from matplotlib.path import Path
from matplotlib.patches import Polygon
from astropy.table import Table

from starbug2.constants import MASK_MAIN_TABLE_PATH
from starbug2.utils import tab2array, colour_index, fill_nan

class Mask(object):
    colour = 'k'

    @staticmethod
    def from_file(f_name):
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
    def from_string(string):
        """
        method to create a mask object from a string. the string should be
        in the following format:

        [-x XCOL] [-y YCOL] [-l Label] : x1 y1 x2 y2 x3 y3 ...

        :param string: the string to create the mask of.
        :return: the constructed mask.
        """
        label = None
        keys = [None, None]
        colour = 'k'
        opts, coords = string.split(':')
        opts, args = getopt.getopt(opts.split(' '), "c:l:x:y:")
        for opt, opt_arg in opts:
            if opt == "-x":
                keys[0] = opt_arg
            if opt == "-y":
                keys[1] = opt_arg
            if opt == "-l":
                label = opt_arg.replace('_', ' ')
            if opt == "-c":
                colour = opt_arg
        strip_coords = coords.strip().rstrip().split(' ')
        points = (
            np.array(strip_coords, dtype=float).reshape(
                (int(len(strip_coords) / 2), 2)))
        return Mask(points, keys, label=label, colour=colour)


    def __init__(self, bounds, keys, label=None, **kwargs):
        """
        mask constructor

        :param bounds: The path vertices, as an array, masked array or
                       sequence of pairs.
        :param keys: iterable array with 2 elements.
        :param label: the mask label
        :param kwargs: kwargs!
        """

        self.path = Path(bounds)
        if len(keys) == 2:
            self.keys = keys
        else:
            raise Exception
        self.label = label
        
        if "colour" in kwargs:
            self.colour = kwargs.get("colour")


    def apply(self, data_table):
        """
        applies a data table based off the masks keys.
        :param data_table: the table to apply the mask on.
        :return: length-N bool array
        :rtype: ndarray
        """
        d = fill_nan(colour_index(data_table, self.keys))
        return self.path.contains_points(tab2array(d))


    def plot(self, axis, **kwargs):
        """
        plots a polygon onto the axis.
        :param axis: the axis to plot the polygon onto.
        :param kwargs: arbitrary polygon parameters.
        :return: None
        """
        patch = Polygon(
            self.path.vertices,
            label=self.label.replace('_', ' ') if self.label else None,
            fill=False, edgecolor=self.colour, **kwargs)
        axis.add_patch(patch)
        



if __name__== "__main__":
    """
    main method if you ran mask object.
    """
    mask_string = "-yF115W -xF115W-F200W -lTestCut 0 20 1 21 1 24 0 24"
    table = Table.read(MASK_MAIN_TABLE_PATH, format="fits").filled(np.nan)
    mask = Mask.from_string(mask_string)
    masked_table = mask.apply(table)
    import matplotlib.pyplot as plt
    tt = colour_index(table, ("F115W-F200W", "F115W"))
    plt.scatter(tt["F115W-F200W"], tt["F115W"], c='k', lw=0, s=1)
    mask.plot(plt.gca(), fill=False, edgecolor="blue", label="test")
    plt.legend()
    plt.show()

