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
        Tuple[np.array, np.array, np.array or None, np.array]):
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

    @abstractmethod
    def artificial_stars(self) -> int:
        # noinspection SpellCheckingInspection
        """
        Execute the automated artificial star testing and completeness
        pipeline.

        This routine injects synthetic point spread function (PSF) source
        profiles into the active observation framework across multiple
        configuration slices to empirically estimate target detection
        completeness thresholds, stellar recovery fractions, and photometric
        parameter variability.

        . note::
            * Flux calculations are normalized automatically into Jansky units
              if the primary FITS image headers track surface brightness
              profiles in Mega-Janskys per steradian (MJy/sr).
            * Background matrices must be explicitly calculated and bound to
              `self.background.data` prior to execution to handle
              background-subtracted PSF fitting accurately.

        :return: Execution status code (0 for clean completion).
        :rtype: int
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
    def options(self) -> dict[str, int | float | str]:
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
