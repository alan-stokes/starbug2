"""
Miscellaneous functions...
"""

import os, stat, numpy as np
from typing import List, Optional, TextIO, Dict, Tuple, Any

from starbug2.constants import (
    JWST_MIRI_APCORR_0010_FITS_URL, JWST_NIRCAM_APCORR_0004_FITS_URL,
    JWST_MIRI_ABVEGA_OFFSET_URL, JWST_NIRCAM_ABVEGA_OFFSET_URL, NIRCAM,
    SHORT, WEBBPSF_PATH_ENV_VAR, LONG, FITS_EXTENSION, FILE_NAME, FILTER,
    OBS, VISIT, DETECTOR, EXPOSURE, STARBUG_DATA_DIR)
from starbug2.constants import STAR_BUG_MIRI
from starbug2.filters import STAR_BUG_FILTERS, FilterStruct
from astropy.io import fits
from astropy.table import Table

from starbug2.matching.generic_match import GenericMatch
from starbug2.starbug import StarbugBase
from starbug2.utils import (
    printf, wget, puts, Loading, p_error, split_file_name)

# A clear, Type Alias for the deep data nested structure Format:
# Dict[KeyType, ValueType], str mapping being FILTER, OBS, VISIT, DETECTOR
ExposureMapping = (
    Dict[Optional[int], Dict[Optional[int], Dict[Optional[int],
    Dict[Optional[int], List[fits.HDUList]]]]])


##########################
# One time run functions #
##########################
def init_starbug() -> None:
    """
    Initialise Starbug..
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
    generate_psfs()

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
def generate_psfs() -> None:
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
                if filter_data.length == SHORT:
                    detectors: List[Optional[str]] = [
                        "NRCA1","NRCA2","NRCA3","NRCA4","NRCB1","NRCB2",
                        "NRCB3","NRCB4"]
                else:
                    detectors = ["NRCA5","NRCB5"]
            else:
                detectors = [None]

            det: str
            for det in detectors:
                load.msg = "%6s %5s" % (filter_string, det)
                load.show()
                psf: fits.PrimaryHDU = generate_psf(filter_string, det, None)
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
        fov_pixels: Optional[int] = None) -> fits.PrimaryHDU:
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

    # ABS again, why are we importing here?
    import stpsf

    # define types
    psf: Optional[fits.PrimaryHDU] = None
    model: Optional[stpsf.JWInstrument] = None

    # ensure fov pixels is greater than 0
    if fov_pixels is not None and fov_pixels <= 0:
        fov_pixels = None

    if filter_string in list(STAR_BUG_FILTERS.keys()):
        the_filter = STAR_BUG_FILTERS.get(filter_string)
        if detector is None:
            if the_filter.instr == NIRCAM and the_filter.length == SHORT:
                detector = "NRCA1"
            elif the_filter.instr == NIRCAM and the_filter.length == LONG:
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
                    model.calc_psf(fov_pixels=fov_pixels)["DET_SAMP"])
                psf: fits.PrimaryHDU = (
                    fits.PrimaryHDU(
                        data=image_hdu.data, header=image_hdu.header))
            except (KeyError, AttributeError, ValueError) as e:
                p_error("\x1b[2KSomething went wrong with: %s %s with "
                        "error %s\n" % (filter_string, detector, str(e)))
            except RuntimeError as e:
                p_error(f"during detector on {the_filter.instr} with wave"
                        f" length {filter_string}. something failed. {str(e)}")
        else:
            p_error(
                "Unable to determine instrument from filter_string '%s'\n" %
                filter_string)
    else:
        p_error("Unable to locate '%s' in JWST filter list\n" % filter_string)
    return psf


def generate_runscript(
        f_names: List[str], args: str = "starbug2 ") -> None:
    """
    generate the run script

    :param f_names: file names
    :type f_names: list of str.
    :param args: arguments
    :type args: str
    :return: None
    """
    runfile: str = "./run.sh"
    fits_files: List[fits.HDUList] = []

    fp: TextIO = open(runfile, "w")
    fp.write("#!/bin/bash\n")
    fp.write("CMDS=\"-vf\"\n")
    for f_name in f_names:
        if os.path.exists(f_name):
            d_name: str
            name: str
            ext: str
            d_name, name, ext = split_file_name(f_name)
            if ext == FITS_EXTENSION:
                fits_file: fits.HDUList = fits.open(f_name)
                fits_file[0].header[FILE_NAME] = f_name
                fits_files.append(fits_file)
            else:
                p_error(
                    "file %s must be type '.fits' not '%s'\n" % (name, ext))
        else:
            p_error("file \x1b[1;31m%s\x1b[0m not found\n" % f_name)

    sorted_exposures: ExposureMapping = sort_exposures(fits_files)

    # print exps.
    for _, obs in sorted_exposures.items():
        for _, visits in obs.items():
            for _, destinations in visits.items():
                for _, exps in destinations.items():
                    s = f"{args}${{CMDS}} -n{len(exps)} "
                    for exp in exps:
                        s += "%s " % exp[0].header[FILE_NAME]
                    fp.write("%s\n" % s)
    fp.close()
    os.chmod(
        runfile, stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR | stat.S_IRGRP |
                 stat.S_IROTH)
    printf("->%s\n" % runfile)


# ABS DOES THIS METHOD NEED TO EXIST?
def calc_instrumental_zero_point(
        psf_table: Table,
        ap_table: Table,
        filter_string: Optional[str] = None) -> Tuple[np.array, np.array]:
    """
    calculates the zero points.

    :param psf_table: the psf table
    :type psf_table: astropy.Table
    :param ap_table: the ap table
    :type ap_table: astropy.Table
    :param filter_string: the filter string
    :type filter_string: str
    :return: tuple of mean zero point and its standard deviation
    :rtype: (np.array, np.array)
    """
    if (filter_string is None
            and not (filter_string := psf_table.meta.get(FILTER))):
        p_error("Unable to determine filter, set with '--set FILTER=F000W'.\n")
        return None, None
    printf("Calculating instrumental zero point %s.\n" % filter_string)

    matcher: GenericMatch = GenericMatch(
        threshold=0.1, col_names=["RA", "DEC", filter_string])
    matched: Table = matcher([psf_table, ap_table], join_type="and")
    dist: np.array = np.array(
        (matched["%s_2" % filter_string]
         - matched["%s_1" % filter_string]).value)
    instr_zp: np.array = np.nanmedian(dist)
    zp_std: np.array = np.nanstd(dist)
    printf("-> zp=%.3f +/- %.2g\n" % (float(instr_zp), float(zp_std)))
    return instr_zp, zp_std


def sort_exposures(catalogues: List[fits.HDUList]) -> ExposureMapping:
    """
     Given a list of catalogue files, this will return the fitsHDULists as a
     series of nested dictionaries sorted by:
    >   BAND
    >   OBSERVATION ID
    >   VISIT ID
    >   DETECTOR                -- These two have been switched
    >   DITHER (EXPOSURE)       -- These two have been switched

    :param catalogues: the catalogues to sort exposures of.
    :type catalogues: list of fits.HDUList
    :return: a dictionary of sorted catalogues
    :rtype: ExposureMapping
    """
    out: ExposureMapping = {}
    for cat in catalogues:
        info = exp_info(cat)

        if info[FILTER] not in out.keys():
            out[info[FILTER]] = {}

        if info[OBS] not in out[info[FILTER]].keys():
            out[info[FILTER]][info[OBS]] = {}

        if info[VISIT] not in out[info[FILTER]][info[OBS]].keys():
            out[info[FILTER]][info[OBS]][info[VISIT]] = {}

        if (info[DETECTOR] not in
            out[info[FILTER]][info[OBS]][info[VISIT]].keys()):
            out[info[FILTER]][
                info[OBS]][info[VISIT]][info[DETECTOR]] = []
        out[info[FILTER]][
            info[OBS]][info[VISIT]][info[DETECTOR]].append(cat)
    return out


def parse_mask(string, table) -> np.ndarray:
    """
    Parse a commandline mask string to be passed into a matching routine
    Example: --mask=F444W!=nan

    :param string: Raw mask sting to be parsed
    :type string: str
    :param table: Table to work on
    :type table: astropy.table.Table
    :return: Boolean mask array to index into a table or array
    :rtype: np.ndarray
    """
    mask: Optional[np.ndarray] = None

    col_name: str
    for col_name in table.colnames:
        string: str = string.replace(col_name, "table[\"%s\"]" % col_name)

    try:
        mask = eval(string)
        if not isinstance(mask, np.ndarray):
            raise Exception
    except NameError as e:
        p_error("Unable to create mask: %s\n" % repr(e))
    except Exception as e:
        p_error(repr(e))

    return mask


def exp_info(hdu_list) -> Dict[str, int]:
    """
    Get the exposure information about a hdu list
    :param hdu_list: HDUList or ImageHDU or BinTableHDU
    :return: dictionary of relevant information
    (filter, obs, visit exposure, detector)
    :rtype dict(str, Optional[int])
    """
    info: Dict[str, int] = {
        FILTER : None,
        OBS : 0,
        VISIT : 0,
        EXPOSURE : 0,
        DETECTOR : None
    }

    if type(hdu_list) in (fits.ImageHDU, fits.BinTableHDU):
        hdu_list: fits.HDUList = fits.HDUList(hdu_list)

    for hdu in hdu_list:
        for key in info:
            if key in hdu.header:
                info[key] = hdu.header[key]
    return info