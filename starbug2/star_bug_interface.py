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
from abc import ABC, abstractmethod
from typing import Optional, Tuple
from astropy.table import Table
from astropy.io.fits import PrimaryHDU, ImageHDU, HDUList, Header

import numpy as np


class StarBugInterface(ABC):

    @abstractmethod
    def log(self, msg: str) -> None:
        """
        Print message if in verbose mode

        :param msg: Message to print out
        :type msg: str
        :return: None
        """
        pass

    @abstractmethod
    def load_image(self, f_name: str) -> None:
        """
        Given f_name, load the image into starbug to be worked on.

        :param f_name: Filename of fits image (with any number of extensions).
            If using a non-standard HDU index, set the name or index of the
            extension with "HDU_NAME=XXX" in the parameter file.
        :type f_name: str
        :return: None
        """
        pass

    @abstractmethod
    def load_ap_file(self, f_name:Optional[str]=None) -> None:
        """
        Load an AP_FILE to be used during photometry

        :param f_name: Filename for fits table containing source coordinates.
         These coordinates can be x-centroid / y-centroid, x_init / y_init,
          x_0, y_0 or RA/DEC. The latter is used if starbug gets "USE_WCS=1"
          in the parameter file.
        :type f_name: str
        :return: None
        """
        pass

    @abstractmethod
    def load_bgd_file(self, f_name: Optional[str]=None) -> None:
        """
        Load a BGD_FILE to be used during photometry

        :param f_name: Filename of fits image the same dimensions as the
                       main image
        :type f_name: str
        :return: None
        """
        pass

    @abstractmethod
    def load_psf(self, f_name: Optional[str]=None) -> int:
        """
        Load a PSF_FILE to be used during photometry

        :param f_name: Filename of a PSF fits image
        :type f_name: str
        :return: the status
        :rtype int
        """
        pass

    @abstractmethod
    def prepare_image_arrays(self) -> (
        Tuple[np.ndarray , np.ndarray, np.ndarray | None, np.ndarray]):
        """
        Make a copy of the original image, and prepare the other image arrays

        :return: tuple of image, error, bgd, mask
        :rtype: tuple of int /float, np.array, np.array or None, int
        """
        pass


    @abstractmethod
    def detect(self) -> int:
        """
        Full source detection routine. Saves the result as a table
        self._detections
        :return: status
        :rtype: int
        """
        pass


    @abstractmethod
    def aperture_photometry(self) -> int:
        """
        executes aperture photometry
        :return: 0 for success 1 for failure
        :rtype int
        """
        pass


    @abstractmethod
    def bgd_estimate(self) -> int:
        """
        Estimate the background of the active image
        Saves the result as an ImageHDU self._background
        :return: the status.
        :rtype: int
        """
        pass


    @abstractmethod
    def bgd_subtraction(self) -> int:
        """
        Internally subtract a background array from an image array
        :return: 0 for success, 1 otherwise
        :rtype int.
        """
        pass


    @abstractmethod
    def photometry_routine(self) -> int:
        """
        Full photometry routine
        Saves the result as a table self._psf_catalogue,
        Additionally it appends a residual Image onto the
        self._residuals HDUList

        :return: 0 for success, 1 otherwise
        :rtype int
        """
        pass


    @abstractmethod
    def source_geometry(self) -> None:
        """
        Calculate source geometry stats for a given image and source list
        :return: None
        """
        pass

    @abstractmethod
    def verify(self) -> int:
        """
        This simple function verifies that everything necessary has been
        loaded properly

        :return: int where 0 on success, 1 on fail
        :rtype int
        """
        pass


    @property
    @abstractmethod
    def header(self) -> Header:
        """
        Construct relevant base header information for routine products

        :return:  Header file containing a series of relevant information
        :rtype: Header
        """
        pass


    @property
    @abstractmethod
    def info(self) -> dict[str, str]:
        """
        Get some useful information from the image header file.

        :return: extracted keys and elements from the image header.
        :rtype: dict of str, to str.
        """
        pass


    @property
    @abstractmethod
    def main_image(self) -> ImageHDU | PrimaryHDU:
        # noinspection SpellCheckingInspection
        """
        automagically find the main image array to use
        Order of importance is:
        > self._nHDU (if set)
        > param[ HDUNAME ]
        > SCI, BGD, RES
        > first ImageHDU
        > first ImageHDU
        > image[0]

        :return: the main image array.
        :rtype: HDUList
        """
        pass


    @property
    @abstractmethod
    def filter(self) -> str | None:
        pass


    @property
    @abstractmethod
    def n_hdu(self) -> int:
        pass


    @property
    @abstractmethod
    def image(self) -> HDUList | None:
        pass


    @image.setter
    @abstractmethod
    def image(self, new_image: HDUList) -> None:
        pass


    @property
    @abstractmethod
    def psf_catalogue(self) -> Table | None:
        pass


    @property
    @abstractmethod
    def psf(self)-> np.ndarray | None:
        pass


    @property
    @abstractmethod
    def f_name(self) -> str | None:
        pass


    @property
    @abstractmethod
    def detections(self)-> Table | None:
        pass


    @detections.setter
    @abstractmethod
    def detections(self, new_detections: Table) -> None:
        pass


    @property
    @abstractmethod
    def out_dir(self) -> str | None:
        pass
