"""
Core routines for StarbugII.
"""
from typing import Optional

import numpy as np
from astropy.table import Table, QTable, hstack
from photutils.detection import DAOStarFinder

from starbug2.constants import X_CENTROID, Y_CENTROID
from starbug2.utils import Loading, printf, p_error



class SourceProperties:
    status: int = 0

    def __init__(self, image: Optional[np.ndarray],
                 source_list: Optional[Table], verbose: int | bool=1) -> None:
        """
        source properties.

        :param image: the image
        :type image: numpy.ndarray or None
        :param source_list: the source list
        :type source_list:  astropy.Table or None
        :param verbose: int for verbose
        :type verbose: int
        """
        self._image: Optional[np.ndarray] = image
        self._source_list: Optional[Table] = None
        self._verbose: int | bool = verbose

        if source_list and type(source_list) in (Table, QTable):
            if len({X_CENTROID, Y_CENTROID} & set(source_list.colnames)) == 2:
                self._source_list = (
                    Table(source_list[[X_CENTROID, Y_CENTROID]]))
            elif len({"x_0", "y_0"} & set(source_list.colnames)) == 2:
                self._source_list = Table(source_list[["x_0", "y_0"]])
                assert self._source_list is not None
                self._source_list.rename_columns(
                    ("x_0", "y_0"), (X_CENTROID, Y_CENTROID))
            else:
                p_error("no positional columns in source list\n")
        else:
            p_error("bad source list type: %s\n" % type(source_list))


    def __call__(self, do_crowd: int=1, n_closest_sources: int = 10,
                 full_width_half_max: float = 2.0) -> Table:
        """
        trigger source properties

        :param do_crowd: int check for doing crowd
        :type do_crowd: int
        :param n_closest_sources: the number of closest sources.
        :type n_closest_sources: int
        :param full_width_half_max: the full width half max.
        :type full_width_half_max: float
        """

        out: Table = Table()

        ## This can be slow
        if do_crowd:
            out = hstack(
                (out, Table([self.calculate_crowding(n_closest_sources)],
                            names=["crowding"])))

        out = hstack((out, self.calculate_geometry(full_width_half_max)))
        return out

    def calculate_crowding(
            self, n_closest_sources: int = 10) -> np.ndarray | None:
        """
        Crowding Index: Sum of magnitude of separation of n closest sources

        :param n_closest_sources: the number of closest sources.
        :type n_closest_sources: int
        """
        if self._source_list is None:
            p_error("no source list\n")
            return None

        crowd: np.ndarray = np.zeros(len(self._source_list))
        load: Loading = Loading(
            len(self._source_list), msg="calculating crowding", res=10)

        for i, src in enumerate([self._source_list]):
            i: int
            src: Table
            dist: np.ndarray = np.sqrt(
                (src[X_CENTROID] - self._source_list[X_CENTROID]) ** 2
                + (src[Y_CENTROID] - self._source_list[Y_CENTROID]) ** 2)
            dist.sort()
            crowd[i]= sum( dist[1 : n_closest_sources])
            load()
            if self._verbose:
                load.show()
        return crowd

    def calculate_geometry(
            self, full_width_half_max: float=2.0) -> Table | None:
        """
        calculate geometry

        :param full_width_half_max: the full width half max.
        :type full_width_half_max: float
        :return the result of geometry
        :rtype Table
        """
        if self._source_list is None:
            p_error("no source list\n")
            return None
        if self._verbose:
            printf("-> measuring source geometry\n")
        xy_coords: np.ndarray = np.array(
            (self._source_list[X_CENTROID], self._source_list[Y_CENTROID])).T

        dao_find: DAOStarFinder = DAOStarFinder(
            -np.inf, full_width_half_max, sharplo=-np.inf, sharphi=np.inf,
            roundlo=-np.inf, roundhi=np.inf, xycoords=xy_coords,
            peakmax=np.inf)

        # ABS protected access. yuck
        return dao_find._get_raw_catalog(self._image).to_table()

