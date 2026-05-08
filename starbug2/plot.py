"""
A collection of plotting functions
"""
import os
import numpy as np
from astropy.visualization import ZScaleInterval
from scipy.interpolate import RegularGridInterpolator
from multiprocessing import Pool

from starbug2.constants import CAT_NUM, URL_DOCS, FILTER, PLOT_MAIN_TABLE_PATH
import matplotlib.image as mpimg
from matplotlib.colors import LinearSegmentedColormap

from astropy.wcs import WCS

import starbug2
from starbug2 import utils

# try to import pyplot as plt.
try: import matplotlib.pyplot as plt
except ImportError:
    from matplotlib import use; use("TkAgg")
    import matplotlib.pyplot as plt


def load_style(f_name):
    """
    Load a pyplot style sheet

    :param f_name: Filename of the style sheet
    :type f_name: str
    :return: None
    """
    if os.path.exists(f_name):
        plt.style.use(f_name)
    else:
        utils.p_error("Unable to load style sheet \"%s\"\n" % f_name)


def get_point_density(x, y, bins=30):
    """
    get point density

    :param x: x coord
    :param y: y coord
    :param bins: the bins.
    :return: densities.
    :rtype: ?????
    """
    f = _generate_regular_grid_interpolator(x, y, bins)

    with Pool(processes=8) as pool:
        dens = pool.map(f, zip(x, y))
    return dens


def _generate_regular_grid_interpolator(x, y, bins):
    """
    generates a regular grid interpolator

    :param x: x coord
    :param y: y coord
    :param bins: the bins.
    :return: the configured RegularGridInterpolator
    :rtype: a RegularGridInterpolator
    """
    hist, _x, _y = np.histogram2d(x, y, bins=bins)
    xx = np.linspace(min(x), max(x), bins)
    yy = np.linspace(min(y), max(y), bins)
    return RegularGridInterpolator((xx, yy), hist)


def plot_test(axes):
    """
    Just plot the starbug image

    :param axes: Ax to plot into
    :type axes: plt.axes
    :return: The working axes
    :rtype: plt.axes
    """
    if not utils.wget(URL_DOCS, f_name="/tmp/sb.png"):
        axes.imshow(mpimg.imread("/tmp/sb.png"))
        axes.set_title("starbug2 v%s" % utils.get_version())
        axes.set_axis_off()
    return axes


def plot_cmd(
        tab, colour, mag, axis=None, col=None, hess=True, x_lim=None,
        y_lim=None, **kwargs):
    """
    plot command.

    :param tab: ???
    :param colour: ???
    :param mag: ???
    :param axis: axis
    :param col: ???
    :param hess: ???
    :param x_lim: ???
    :param y_lim: ???
    :param kwargs: ???
    :return: axes.
    :rtype: plt.axes
    """
    tt = utils.colour_index(tab, (colour,mag))
    mask =~ (tt[colour].mask | tt[mag].mask)
    cc = tt[colour][mask]
    mm = tt[mag][mask]

    if not axis:
        plt.subplots(1)

    if x_lim is None:
        x_lim = (np.nanmin(cc),np.nanmax(cc))
    if y_lim is None:
        y_lim = (np.nanmin(mm),np.nanmax(mm))

    mask = (
        (cc >= x_lim[0]) & (cc <= x_lim[1]) &
        (mm >= y_lim[0]) & (mm <= y_lim[1]))
    cc = cc[mask]
    mm = mm[mask]

    if col is None:
        col = plt.rcParams["axes.prop_cycle"].by_key()["color"][0]
    cmap = LinearSegmentedColormap.from_list(
        "", [plt.rcParams["axes.prop_cycle"].by_key()["color"][0], col])

    if hess:
        bins = 100
        f = _generate_regular_grid_interpolator(cc, mm, bins)
        col = [f([X,Y]) for X,Y in zip(cc,mm)]
    pyplot_kw = {"lw":0,"s":3}
    pyplot_kw.update(kwargs)
    ax.scatter(cc, mm, c=col, cmap=cmap, **pyplot_kw)

    ax.set_xlabel(colour)
    ax.set_ylabel(mag)
    ax.set_xlim(x_lim)
    ax.set_ylim(*y_lim[::-1])
    return ax


def plot_inspect_source(src, images):
    """
    Show a source in an array of images

    :param src: Input source to look at
    :type src: astropy.Table Row
    :param images: List of fits images to inspect
    :type images: list of fits.Image
    :return: The figure
    :rtype: plt.figure
    """

    n = len(images)
    figure, axs = plt.subplots(1, n, figsize=(1.7 * n, 2))
    if n == 1:
        axs=[axs]
    images = sorted(
        images, key=lambda a:
            list(starbug2.filters.keys()).index(a.header[FILTER]))

    #arcsec?
    size = 0.1
    for n, (im, axis) in enumerate(zip(images,axs)):
        wcs = WCS(im)
        x, y = wcs.all_world2pix(
            np.ones(2) * src["RA"],
            np.array([1, 1 + (size / 3600)]) * src["DEC"], 0)
        dp = np.sqrt((x[1] - x[0]) ** 2 + (y[1] - y[0]) ** 2)

        x_min = max(0, int(np.floor(x[0] - dp)))
        x_max = min(im.data.shape[1] - 1, int(np.ceil(x[0] + dp)))
        y_min = max(0, int(np.floor(y[0] - dp)))
        y_max = min(im.data.shape[0] - 1, int(np.ceil(y[0] + dp)))

        dat = im.data[
              min(y_min, y_max): max(y_min, y_max),
              min(x_min, x_max): max(x_min, x_max)]
        if all(dat.shape):
            ax.imshow(ZScaleInterval()(dat), cmap="Greys_r", origin="lower")
            ax.text(0, 0, im.header.get(FILTER), c="white")

        ax.set_axis_off()
        figure.suptitle(src[CAT_NUM][0])
    figure.tight_layout()
        
    return figure


if __name__=="__main__":
    from astropy.table import Table
    fig, ax = plt.subplots(1)
    t = Table().read(PLOT_MAIN_TABLE_PATH)
    plot_cmd(t, "F115W-F200W","F200W", ax=ax)
    plt.show()
