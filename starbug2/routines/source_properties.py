"""
Core routines for StarbugII.
"""
import numpy as np
from astropy.table import Table, QTable, hstack
from photutils.detection import DAOStarFinder

from starbug2.constants import X_CENTROID, Y_CENTROID
from starbug2.utils import Loading, printf, p_error



class SourceProperties:
    status = 0
    def __init__(self, image, source_list, verbose=1):
        """
        source properties.

        :param image: the image
        :type image: numpy.ndarray
        :param source_list: the source list
        :type source_list:  astropy.Table
        :param verbose: int for verbose
        :type verbose: int
        """
        self._image = image
        self._source_list = None
        self._verbose = verbose

        if source_list and type(source_list) in (Table, QTable):
            if len({X_CENTROID, Y_CENTROID} & set(source_list.col_names)) == 2:
                self.source_list = Table(source_list[[X_CENTROID, Y_CENTROID]])
            elif len({"x_0", "y_0"} & set(source_list.col_names)) == 2:
                self.source_list = Table(source_list[["x_0", "y_0"]])
                self.source_list.rename_columns(
                    ("x_0", "y_0"), (X_CENTROID, Y_CENTROID))
            else:
                p_error("no positional columns in source list\n")
        else:
            p_error("bad source list type: %s\n" % type(source_list))


    def __call__(self, do_crowd=1, **kwargs):
        """
        trigger source properties

        :param do_crowd: int
        :type do_crowdL int
        :param kwargs: extra args
        :type kwargs: dict.
        """

        out = Table()

        ## This can be slow
        if do_crowd:
            out = hstack(
                (out, Table([self.calculate_crowding(**kwargs)],
                            names=["crowding"])))

        out = hstack((out, self.calculate_geometry(**kwargs)))
        return out

    def calculate_crowding(self, n_closest_sources=10, **kwargs):
        """
        Crowding Index: Sum of magnitude of separation of n closest sources

        :param n_closest_sources: the number of closest sources.
        :type n_closest_sources: int
        :param kwargs: extra arguments.
        :type kwargs: extra args.
        """
        if self.source_list is None:
            p_error("no source list\n")
            return None

        crowd = np.zeros(len(self.source_list))
        load = Loading(
            len(self.source_list), msg="calculating crowding", res=10)

        for i, src in enumerate(self._source_list):
            dist = np.sqrt(
                (src[X_CENTROID] - self.source_list[X_CENTROID]) ** 2
                + (src[Y_CENTROID] - self.source_list[Y_CENTROID]) ** 2)
            dist.sort()
            crowd[i]= sum( dist[1 : n_closest_sources])
            load()
            if self._verbose:
                load.show()
        return crowd

    def calculate_geometry(self, full_width_half_max=2.0, **kwargs):
        """
        calculate geometry

        :param full_width_half_max: the full width half max.
        :type full_width_half_max: float
        :param kwargs: extra args
        :type kwargs: dict.
        """
        if self.source_list is None:
            p_error("no source list\n")
            return None
        if self._verbose:
            printf("-> measuring source geometry\n")
        xy_coords = np.array(
            (self.source_list[X_CENTROID], self.source_list[Y_CENTROID])).T

        dao_find = DAOStarFinder(
            -np.inf, full_width_half_max, sharplo=-np.inf, sharphi=np.inf,
            roundlo=-np.inf, roundhi=np.inf, xycoords=xy_coords,
            peakmax=np.inf)

        # ABS protected access. yuck
        return dao_find._get_raw_catalog(self._image).to_table()

