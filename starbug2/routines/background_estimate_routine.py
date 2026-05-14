"""
Core routines for StarbugII.
"""
import numpy as np
from photutils.background import Background2D, BackgroundBase
from photutils.aperture import CircularAperture

from starbug2.constants import X_CENTROID, Y_CENTROID
from starbug2.utils import Loading, printf,  warn

class BackGroundEstimateRoutine(BackgroundBase):
    def __init__(
            self, source_list, box_size=2, full_width_half_max=2, sig_sky=2,
            bgd_r=-1, profile_scale=1, profile_slope=0.5, verbose=0, bgd=None):
        """
        Diffuse background emission estimator run by starbug.

        :param source_list: List of sources in the image in a table
                            containing X_CENTROID Y_CENTROID
        :type source_list: astropy.table.Table
        :param box_size: Size of the kernel to pass over the image in
                         un-sharp masking
        :type box_size: int
        :param full_width_half_max: Source full width half maximum in the image
        :type full_width_half_max: float
        :param sig_sky: Sigma threshold for background clipping
        :type sig_sky: float
        :param bgd_r: Fixed aperture mask radius around each source
        :type bgd_r: float
        :param profile_scale: Scaling factor for the aperture mask radius
                              profile
        :type profile_scale: float
        :param profile_slope: Slope of the aperture mask radius profile
        :type profile_slope: float
        :param verbose: Set whether to print verbose output information
        :type verbose: bool
        """
        self._source_list = source_list
        self._box_size = box_size
        self._full_width_half_max = full_width_half_max
        self._sig_sky = sig_sky
        self._bgd_r = bgd_r
        self._a = profile_scale
        self._b = profile_slope
        self._verbose = verbose
        self._bgd = bgd
        super().__init__()

    def calc_peaks(self, im):
        """
        Determine peak pixel value for each source in xy
        :param im: ??????
        :return: peaks
        :rtype: np.array
        """
        x = self._source_list[X_CENTROID]
        y = self._source_list[Y_CENTROID]
        apertures = CircularAperture(np.array((x,y)).T, 2).to_mask()
        peaks = np.full(len(x), np.nan)
        for i, mask in enumerate(apertures):
            peaks[i] = np.nanmax(mask.multiply(im))
        return peaks

    def log(self, msg):
        """
        log this message.
        :param msg: the message to log
        :return: None
        """
        if self._verbose:
            printf(msg)

    def __call__(self, data, axis=None, masked=False, output=None):
        """
        does background estimation routine.

        :param data: the data to process
        :param axis: the axis
        :param masked: if the data is masked
        :param output: if the data should be outputted
        :return: the new background 2D object.
        :rtype: Background2D
        """
        if self._source_list is None or data is None:
            return self._bgd
        _data = np.copy(data)
        x_grid, y_grid = np.ogrid[:data.shape[1], :data.shape[0]]

        default_r = 2 * self._full_width_half_max

        if self._bgd_r and self._bgd_r>0:
            self.log(
                "-> using BGD_R=%g masking aperture radii\n" % self._bgd_r)
            rlist= self._bgd_r * np.ones(len(self._source_list))

        else:
            if "flux" in self._source_list.colnames:
                self.log("-> calculating source aperture mask radii\n")
                sky = (
                    self._source_list["sky"]
                    if "sky" in self._source_list.colnames else 1.0)
                rlist = (
                    self._a * self._full_width_half_max * (
                        np.log(self._source_list["flux"] / sky))
                    ** self._b)
                rlist[np.isnan(rlist)] = default_r

                if output:
                    with open(output,'w') as fp:
                        for i in range(len(rlist)):
                            fp.write(
                                "circle %f %f %f #color=green;" % (
                                    1 + self._source_list[i][X_CENTROID],
                                    1 + self._source_list[i][Y_CENTROID],
                                    rlist[i]))
                            fp.write(
                                "annulus %f %f %f %f #color=white;" % (
                                    1 + self._source_list[i][X_CENTROID],
                                    1 + self._source_list[i][Y_CENTROID],
                                    1.5 * rlist[i],
                                    1.5 * rlist[i] + 1))
                    self.log("-> exporting check file \"%s\"\n" % output)
            else:
                warn("Unable to calculate aperture mask sizes, "
                     "add '-A' to starbug command.\n")
                rlist = default_r * np.ones(len(self._source_list))

        dimension = 50
        load = Loading(len(self._source_list), msg="masking sources", res=10)
        for r,src in zip(rlist, self._source_list):

            rin = 1.5 * r
            rout = rin + 1

            x = round(src[X_CENTROID])
            y = round(src[Y_CENTROID])
            _X = x_grid[
                 max( x - dimension, 0) : min(x + dimension, data.shape[1])]
            _Y = y_grid[
                 :,max( y - dimension, 0) : min(y + dimension, data.shape[0])]

            radius = np.sqrt(
                (_X-src[X_CENTROID]) ** 2 + (_Y -src[Y_CENTROID]) ** 2)

            mask = (radius < r)
            annuli_mask = ((radius > rin) & (radius < rout))

            tmp = _data[_Y, _X]
            tmp[mask] = np.median(data[_Y, _X][annuli_mask])
            _data[_Y,_X] = tmp

            load()

            ## This will slow the thing down quite a lot
            if self._verbose:
                load.show()
        if self._verbose:
            printf("-> estimating bgd2d\n")
        self._bgd = Background2D(_data, self._box_size).background
        return self._bgd

    def calc_background(self,data, axis=None, masked=None):
        if self._bgd is None:
            self.__call__(data)
        return self._bgd