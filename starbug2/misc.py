"""
Miscellaneous functions...
"""

import os, stat, numpy as np
from starbug2.constants import (
    JWST_MIRI_APCORR_0010_FITS_URL, JWST_NIRCAM_APCORR_0004_FITS_URL,
    JWST_MIRI_ABVEGA_OFFSET_URL, JWST_NIRCAM_ABVEGA_OFFSET_URL, NIRCAM,
    SHORT, WEBBPSF_PATH_ENV_VAR, LONG, FITS_EXTENSION, FILE_NAME, FILTER)
from starbug2.constants import MIRI as STAR_BUG_MIRI
from starbug2.filters import filters
from astropy.io import fits

from starbug2.matching.generic_match import GenericMatch
from starbug2.starbug import StarbugBase
from starbug2.utils import (
    printf, wget, puts, Loading, p_error, split_file_name)


##########################
# One time run functions #
##########################
def init_starbug():
    """
    Initialise Starbug..
        - generate PSFs
        - download crds files
    INPUT:
        data_name : data directory
    """
    printf("Initialising StarbugII\n")

    data_name = StarbugBase.get_data_path()

    # noinspection SpellCheckingInspection
    printf("-> using %s=%s\n" % (
        "STARBUG_DATDIR" if os.getenv("STARBUG_DATDIR") else "DEFAULT_DIR",
        data_name))
    generate_psfs()

    _miri_ap_corr = JWST_MIRI_APCORR_0010_FITS_URL

    # noinspection SpellCheckingInspection
    _nircam_ap_corr = JWST_NIRCAM_APCORR_0004_FITS_URL

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
def generate_psfs():
    """
    Generate the psf files inside a given directory

    utilises the star bug data patj to generate the directory to generate info
    :return:
    """
    dname = StarbugBase.get_data_path()
    if os.getenv(WEBBPSF_PATH_ENV_VAR):
        dname = os.path.expandvars(dname)
        if not os.path.exists(dname):
            os.makedirs(dname)

        printf("Generating PSFs --> %s\n"%dname)

        load = Loading(145, msg="initialising")
        load.show()
        for fltr, _f in filters.items():
            if _f.instr == NIRCAM:
                if _f.length == SHORT:
                    detectors = [
                        "NRCA1","NRCA2","NRCA3","NRCA4","NRCB1","NRCB2",
                        "NRCB3","NRCB4"]
                else:
                    detectors = ["NRCA5","NRCB5"]
            else:
                detectors = [None]

            for det in detectors:
                load.msg = "%6s %5s" % (fltr, det)
                load.show()
                psf = generate_psf(fltr, det, None)
                if psf: 
                    psf.writeto(
                        "%s/%s%s.fits" % (
                            dname, fltr, "" if det is None else det),
                        overwrite=True)
                load()
                load.show()

    else:
        p_error(
            "WARNING: Cannot generate PSFs, no environment variable "
            "'WEBBPSF_PATH', please see "
            "https://webbpsf.readthedocs.io/en/latest/installation.html\n")


# noinspection SpellCheckingInspection
def generate_psf(filter_string, detector=None, fov_pixels=None):
    # noinspection SpellCheckingInspection
    """
    Generate a single PSF for JWST

    :param filter_string: the filter string from the default set of filters (
                          e.g. F444W)
    :type filter_string: str
    :param detector: Instrument detector module e.g. NRCA1
    :type detector: str
    :param fov_pixels: size of PSF
    :return: the generated psfs
    :rtype fits.HDUlist
    """

    # ABS again, why are we importing here?
    from webbpsf.webbpsf_core import NIRCam, MIRI
    psf = None
    model = None
    if fov_pixels is not None and fov_pixels <= 0:
        fov_pixels = None

    if filter_string in list(filters.keys()):
        the_filter = filters.get(filter_string)
        if detector is None:
            if the_filter.instr == NIRCAM and the_filter.length == SHORT:
                detector = "NRCA1"
            elif the_filter.instr == NIRCAM and the_filter.length == LONG:
                detector="NRCA5"
            elif the_filter.instr == STAR_BUG_MIRI:
                detector = "MIRIM"

        if the_filter.instr == NIRCAM:
            model = NIRCam()
        elif the_filter.instr == STAR_BUG_MIRI:
            model = MIRI()

        if model:
            model.filter = filter_string
            if detector:
                model.detector = detector
            try:
                # this actually works, as lower stream code will check if
                # fox_pixels is set to None and utilise sensible defaults.
                # so basically bad docing in dependency causes this issue.
                # noinspection PyTypeChecker
                psf = model.calc_psf(fov_pixels=fov_pixels)["DET_SAMP"]
                psf = fits.PrimaryHDU(data=psf.data, header=psf.header)
            except (KeyError, AttributeError, ValueError) as e:
                p_error("\x1b[2KSomething went wrong with: %s %s with "
                        "error %s\n" % (filter_string, detector, str(e)))
        else:
            p_error(
                "Unable to determine instrument from filter_string '%s'\n" %
                filter_string)
    else:
        p_error("Unable to locate '%s' in JWST filter list\n" % filter_string)
    return psf


def generate_runscript(f_names, args="starbug2 "):
    """
    generate the run script

    :param f_names: file names
    :type f_names: list of str.
    :param args: arguments
    :type args: str
    :return: None
    """
    runfile = "./run.sh"
    fits_files = []

    fp = open(runfile, "w")
    fp.write("#!/bin/bash\n")
    fp.write("CMDS=\"-vf\"\n")
    for f_name in f_names:
        if os.path.exists(f_name):
            d_name, name, ext = split_file_name(f_name)
            if ext == FITS_EXTENSION:
                fits_file = fits.open(f_name)
                fits_file[0].header[FILE_NAME] = f_name
                fits_files.append(fits_file)
            else:
                p_error(
                    "file %s must be type '.fits' not '%s'\n" % (name, ext))
        else:
            p_error("file \x1b[1;31m%s\x1b[0m not found\n" % f_name)

    sorted_exposures = sort_exposures(fits_files)

    for band, obs in sorted_exposures.items():
        for ob, visits in obs.items():
            for visit, destinations in visits.items():
                for destination, exps in destinations.items():
                    s = f"{args}${{CMDS}} -n{len(exps)} "
                    for exp in exps:
                        s += "%s " % exp[0].header[FILE_NAME]
                    fp.write("%s\n" % s)
    fp.close()
    os.chmod(
        runfile, stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR | stat.S_IRGRP |
                 stat.S_IROTH)
    printf("->%s\n" % runfile)



def calc_instrumental_zero_point(psf_table, ap_table, filter_string=None):
    """
    calculates the zero points.

    :param psf_table: the psf table
    :param ap_table: the ap table
    :param filter_string: the filter string
    :return: tuple of mean zero point and its standard deviation
    :rtype: (float, float)
    """
    if (filter_string is None
            and not (filter_string := psf_table.meta.get(FILTER))):
        p_error("Unable to determine filter, set with '--set FILTER=F000W'.\n")
        return None
    printf("Calculating instrumental zero point %s.\n" % filter_string)

    m = GenericMatch(threshold=0.1, col_names=["RA", "DEC", filter_string])
    matched = m([psf_table, ap_table], join_type="and")
    dist = np.array(
        (matched["%s_2" % filter_string]
         - matched["%s_1" % filter_string]).value)
    instr_zp = np.nanmedian(dist)
    zp_std = np.nanstd(dist)
    printf("-> zp=%.3f +/- %.2g\n" % (float(instr_zp), float(zp_std)))
    return instr_zp, zp_std



def sort_exposures(catalogues):
    """
     Given a list of catalogue files, this will return the fitsHDULists as a
     series of nested dictionaries sorted by:
    >   BAND
    >   OBSERVATION ID
    >   VISIT ID
    >   DETECTOR                -- These two have been switched
    >   DITHER (EXPOSURE)       -- These two have been switched

    :param catalogues: the catalogues to sort exposures of.
    :return: a dictionary of sorted catalogues
    """
    out = {}
    for cat in catalogues:
        info = exp_info(cat)

        if info[_FILTER] not in out.keys():
            out[info[_FILTER]] = {}

        if info[_OBS] not in out[info[_FILTER]].keys():
            out[info[_FILTER]][info[_OBS]] = {}

        if info[_VISIT] not in out[info[_FILTER]][info[_OBS]].keys():
            out[info[_FILTER]][info[_OBS]][info[_VISIT]] = {}

        if (info[_DETECTOR] not in
            out[info[_FILTER]][info[_OBS]][info[_VISIT]].keys()):
            out[info[_FILTER]][
                info[_OBS]][info[_VISIT]][info[_DETECTOR]] = []
        out[info[_FILTER]][
            info[_OBS]][info[_VISIT]][info[_DETECTOR]].append(cat)
    return out


def parse_mask(string, table):
    """
    Parse an commandline mask string to be passed into a matching routine
    Example: --mask=F444W!=nan

    :param string: Raw mask sting to be parsed
    :type string: str
    :param table: Table to work on
    :type table: astropy.table.Table
    :return: Boolean mask array to index into a table or array
    :rtype: np.ndarray
    """
    mask=None

    for col_name in table.colnames: string=string.replace(
        col_name, "table[\"%s\"]" % col_name)

    try:
        mask = eval(string)
        if not isinstance(mask, np.ndarray):
            raise Exception
    except NameError as e:
        p_error("Unable to create mask: %s\n" % repr(e))
    except Exception as e:
        p_error(repr(e))

    return mask


def exp_info(hdu_list):
    """
    Get the exposure information about a hdu list
    INPUT:  HDUList or ImageHDU or BinTableHDU
    RETURN: dictionary of relevant information:
            >   EXPOSURE, DETECTOR, FILTER
    """
    info={  _FILTER : None,
            _OBS : 0,
            _VISIT : 0,
            _EXPOSURE : 0,
            _DETECTOR : None
            }

    if type(hdu_list) in (fits.ImageHDU, fits.BinTableHDU):
        hdu_list=fits.HDUList(hdu_list)

    for hdu in hdu_list:
        for key in info:
            if key in hdu.header:
                info[key] = hdu.header[key]
    return info