import os
from typing import List, Optional, Any, Final

from starbug2.constants import (
    JWST_MIRI_APCORR_0010_FITS_URL, JWST_NIRCAM_APCORR_0004_FITS_URL,
    JWST_MIRI_ABVEGA_OFFSET_URL, JWST_NIRCAM_ABVEGA_OFFSET_URL, NIRCAM,
    WEBBPSF_PATH_ENV_VAR, DetectorLengths, STARBUG_DATA_DIR)
from starbug2.constants import STAR_BUG_MIRI
from starbug2.filters import STAR_BUG_FILTERS, FilterStruct
from astropy.io import fits
from starbug2.starbug import StarbugBase
from starbug2.utils import printf, wget, puts, Loading, p_error
import stpsf

# noinspection SpellCheckingInspection
# the detector labels for nircam short length detectors
NIRCAM_SHORT_DETECTORS: Final[list[str]] = [
    "NRCA1","NRCA2","NRCA3","NRCA4","NRCB1","NRCB2",
    "NRCB3","NRCB4"]

# noinspection SpellCheckingInspection
# the detector labels for nircam long length detectors.
NIRCAM_LONG_DETECTORS: Final[list[str]] = ["NRCA5","NRCB5"]

##########################
# One time run functions #
##########################
def init_starbug_for_jwst() -> None:
    """
    Initialise Starbug for jwst.
        - generate PSFs
        - download crds files
    INPUT:
        data_name : data directory
    """
    printf("Initialising StarbugII\n")

    data_name: str = StarbugBase.get_data_path()

    # noinspection SpellCheckingInspection
    printf("-> using %s=%s\n" % (
        STARBUG_DATA_DIR if os.getenv(STARBUG_DATA_DIR) else "DEFAULT_DIR",
        data_name))
    _generate_psfs()

    _miri_ap_corr: str = JWST_MIRI_APCORR_0010_FITS_URL

    # noinspection SpellCheckingInspection
    _nircam_ap_corr: str = JWST_NIRCAM_APCORR_0004_FITS_URL

    # noinspection SpellCheckingInspection
    printf("Downloading APPCORR CRDS files. NB: "
           "\x1b[1mTHESE MAY NOT BE THE LATEST!\x1b[0m\n")
    printf("-> %s\n" % _miri_ap_corr)
    printf("-> %s\n" % _nircam_ap_corr)

    # noinspection SpellCheckingInspection
    wget(_miri_ap_corr, "%s/apcorr_miri.fits" % data_name)

    # noinspection SpellCheckingInspection
    wget(_nircam_ap_corr, "%s/apcorr_nircam.fits" % data_name)

    # noinspection SpellCheckingInspection
    printf("Downloading ABVEGA offsets.\n")

    # noinspection SpellCheckingInspection
    wget(JWST_MIRI_ABVEGA_OFFSET_URL,
         "%s/abvegaoffset_miri.asdf" % data_name)

    # noinspection SpellCheckingInspection
    wget(JWST_NIRCAM_ABVEGA_OFFSET_URL,
         "%s/abvegaoffset_nircam.asdf" % data_name)
    puts("Downloading The Junior Colour Encyclopedia of Space\n")

# noinspection SpellCheckingInspection
def _generate_psfs() -> None:
    """
    Generate the psf files inside a given directory

    utilises the star bug data patj to generate the directory to generate info
    :return:
    """
    dname: str = StarbugBase.get_data_path()
    if os.getenv(WEBBPSF_PATH_ENV_VAR):
        dname = os.path.expandvars(dname)
        if not os.path.exists(dname):
            os.makedirs(dname)

        printf("Generating PSFs --> %s\n"%dname)

        load: Loading = Loading(145, msg="initialising")
        load.show()

        # type hitns
        filter_string: str
        filter_data: FilterStruct

        for filter_string, filter_data in STAR_BUG_FILTERS.items():
            if filter_data.instr == NIRCAM:
                if filter_data.length == DetectorLengths.SHORT:
                    detectors: List[Optional[str]] = NIRCAM_SHORT_DETECTORS
                else:
                    detectors = NIRCAM_LONG_DETECTORS
            else:
                detectors = [None]

            det: str
            for det in detectors:
                load.msg = "%6s %5s" % (filter_string, det)
                load.show()
                psf: fits.PrimaryHDU | None = generate_psf(
                    filter_string, det, None)
                if psf:
                    psf.writeto(
                        "%s/%s%s.fits" % (
                            dname, filter_string, "" if det is None else det),
                        overwrite=True)
                load()
                load.show()

    else:
        p_error(
            "WARNING: Cannot generate PSFs, no environment variable "
            "'WEBBPSF_PATH', please see "
            "https://webbpsf.readthedocs.io/en/latest/installation.html\n")


# noinspection SpellCheckingInspection
def generate_psf(
    filter_string: str,
    detector: Optional[str] = None,
    fov_pixels: Optional[int] = None) -> fits.PrimaryHDU | None:
    # noinspection SpellCheckingInspection
    """
    Generate a single PSF for JWST

    :param filter_string: the filter string from the default set of filters (
                          e.g. F444W)
    :type filter_string: str
    :param detector: Instrument detector module e.g. NRCA1
    :type detector: str
    :param fov_pixels: size of PSF
    :type fov_pixels: int
    :return: the generated psfs
    :rtype fits.PrimaryHDU
    """

    # define types
    psf: Optional[fits.PrimaryHDU] = None
    model: Optional[stpsf.JWInstrument] = None

    # ensure fov pixels is greater than 0
    if fov_pixels is not None and fov_pixels <= 0:
        fov_pixels = None

    if filter_string in list(STAR_BUG_FILTERS.keys()):
        the_filter = STAR_BUG_FILTERS.get(filter_string)
        assert the_filter is not None
        if detector is None:
            if (the_filter.instr == NIRCAM
                    and the_filter.length == DetectorLengths.SHORT):
                detector = "NRCA1"
            elif (the_filter.instr == NIRCAM
                    and the_filter.length == DetectorLengths.LONG):
                detector = "NRCA5"
            elif the_filter.instr == STAR_BUG_MIRI:
                detector = "MIRIM"
            else:
                detector = "MIRIM"

        # need to use getattr as these are not found by the IDE automatically.
        mode: stpsf.JWInstrument
        if the_filter.instr == NIRCAM:
            model = getattr(stpsf, "NIRCam")()
        elif the_filter.instr == STAR_BUG_MIRI:
            model = getattr(stpsf, "MIRI")()

        if model:
            model.filter = filter_string
            if detector:
                model.detector = detector
            try:
                # this actually works, as lower stream code will check if
                # fox_pixels is set to None and utilise sensible defaults.
                # so basically bad docing in dependency causes this issue.
                # noinspection PyTypeChecker
                image_hdu: fits.ImageHDU | Any = (
                    model.calc_psf(
                        fov_pixels=fov_pixels,
                        nlambda=the_filter.nlambda)["DET_SAMP"])
                psf: fits.PrimaryHDU = (
                    fits.PrimaryHDU(
                        data=image_hdu.data, header=image_hdu.header))
            except (KeyError, AttributeError, ValueError) as e:
                p_error("\x1b[2KSomething went wrong with: %s %s with "
                        "error %s\n" % (filter_string, detector, str(e)))
            except RuntimeError as e:
                import traceback
                traceback.format_exc()
                p_error(f"during detector on {the_filter.instr} with wave"
                        f" length {filter_string}. something failed. {str(e)}")
        else:
            p_error(
                "Unable to determine instrument from filter_string '%s'\n" %
                filter_string)
    else:
        p_error("Unable to locate '%s' in JWST filter list\n" % filter_string)
    return psf