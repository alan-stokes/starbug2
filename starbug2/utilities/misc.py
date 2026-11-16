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
import os
import stat
import numpy as np

from starbug2.constants import (
    FileExtensions, FILE_NAME, HeaderTags, ImageHeaderTags)
from astropy.io import fits
from starbug2.utilities.utils import printf, p_error, split_file_name

# A clear, Type Alias for the deep data nested structure Format:
# Dict[KeyType, ValueType], str mapping being FILTER, OBS, VISIT, DETECTOR
ExposureMapping = (
    dict[
        int | None, dict[
            int | None, dict[
                int | None, dict[
                    int | None, list[fits.HDUList]]]]])


def generate_runscript(
        f_names: list[str], args: str = "starbug2 ") -> None:
    """
    generate the run script

    :param f_names: file names
    :type f_names: list of str.
    :param args: arguments
    :type args: str
    :return: None
    """
    runfile: str = "./run.sh"
    fits_files: list[fits.HDUList] = []

    for f_name in f_names:
        if os.path.exists(f_name):
            _, name, ext = split_file_name(f_name)
            if ext == FileExtensions.FITS:
                fits_file: fits.HDUList = fits.open(f_name)
                fits_file[0].header[FILE_NAME] = f_name
                fits_files.append(fits_file)
            else:
                p_error(f"file {name} must be type '.fits' not '{ext}'\n")
        else:
            p_error(f"file \x1b[1;31m{f_name}\x1b[0m not found\n")

    sorted_exposures: ExposureMapping = sort_exposures(fits_files)

    with open(runfile, "w") as fp:
        fp.write("#!/bin/bash\n")
        fp.write('CMDS="-vf"\n')

        # print exps.
        for obs in sorted_exposures.values():
            for visits in obs.values():
                for destinations in visits.values():
                    for exps in destinations.values():
                        s = f"{args}${{CMDS}} -n{len(exps)} "
                        for exp in exps:
                            s += f"{exp[0].header[FILE_NAME]} "
                        fp.write(f"{s}\n")

    os.chmod(
        runfile, stat.S_IXUSR | stat.S_IWUSR | stat.S_IRUSR | stat.S_IRGRP |
        stat.S_IROTH)
    printf(f"->{runfile}\n")


def sort_exposures(catalogues: list[fits.HDUList]) -> ExposureMapping:
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

        if info[HeaderTags.FILTER] not in out:
            out[info[HeaderTags.FILTER]] = {}

        if info[HeaderTags.OBS] not in out[info[HeaderTags.FILTER]]:
            out[info[HeaderTags.FILTER]][info[HeaderTags.OBS]] = {}

        if (info[HeaderTags.VISIT] not in
                out[info[HeaderTags.FILTER]][info[HeaderTags.OBS]]):
            out[info[HeaderTags.FILTER]][
                info[HeaderTags.OBS]][info[HeaderTags.VISIT]] = {}

        if (info[ImageHeaderTags.DETECTOR] not in
            out[info[HeaderTags.FILTER]][info[HeaderTags.OBS]][
                info[HeaderTags.VISIT]]):
            out[info[HeaderTags.FILTER]][
                info[HeaderTags.OBS]][info[HeaderTags.VISIT]][
                info[ImageHeaderTags.DETECTOR]] = []
        out[info[HeaderTags.FILTER]][
            info[HeaderTags.OBS]][info[HeaderTags.VISIT]][
            info[ImageHeaderTags.DETECTOR]].append(cat)
    return out


def parse_mask(string, table) -> np.ndarray | None:
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
    mask: np.ndarray | None = None

    col_name: str
    for col_name in table.colnames:
        string: str = string.replace(col_name, f"table[\"{col_name}\"]")

    try:
        mask = eval(string)
        if not isinstance(mask, np.ndarray):
            raise TypeError
    except NameError as e:
        p_error(f"Unable to create mask: {e!r}\n")
    except TypeError as e:
        p_error(repr(e))

    return mask


def exp_info(hdu_list) -> dict[str, int | None]:
    """
    Get the exposure information about a hdu list
    :param hdu_list: HDUList or ImageHDU or BinTableHDU
    :return: dictionary of relevant information
    (HeaderTags.FILTER, obs, visit exposure, detector)
    :rtype dict(str, Optional[int])
    """
    info: dict[str, int | None] = {
        HeaderTags.FILTER: None,
        HeaderTags.OBS: 0,
        HeaderTags.VISIT: 0,
        HeaderTags.EXPOSURE: 0,
        ImageHeaderTags.DETECTOR: None
    }

    if type(hdu_list) in (fits.ImageHDU, fits.BinTableHDU):
        hdu_list: fits.HDUList = fits.HDUList(hdu_list)

    for hdu in hdu_list:
        for key in info:
            if key in hdu.header:
                info[key] = hdu.header[key]
    return info
