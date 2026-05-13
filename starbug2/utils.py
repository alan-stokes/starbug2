import time
import os, sys, numpy as np
from importlib import metadata
from astropy.table import Table, hstack, Column, MaskedColumn, vstack
from astropy.io import fits
from astropy.wcs import WCS
import starbug2
import requests
from importlib.metadata import PackageNotFoundError

from starbug2.constants import (
    CAT_NUM, DEFAULT_COLOUR, RA, DEC, TMP_OUT, FLUX, TMP_FITS,
    FITS_EXTENSION, FILTER, N_MIS_MATCHES, EXIT_SUCCESS, EXIT_FAIL,
    REST_SUCCESS_CODE)

# different print methods (why are we not using loggers?)
printf = sys.stdout.write
p_error = sys.stderr.write
puts = lambda s:printf("%s\n"%s)
s_bold = lambda s: "\x1b[1m%s\x1b[0m" % s
warn = lambda s:p_error("%s%s" % (s_bold("Warning: "), s))

def append_chars(s, n, c):
    """
    append n characters to s.
    :param s: the base string
    :param n: the number of times to add the character
    :param c: the characters to add.
    :return: the adjusted string.
    """
    for _ in range(n):
        s += c
    return s

def repeat_print(n, c):
    """
    prints out a repeated string.
    NOTE: this seems unused.

    :param n: the number of times to repeat
    :param c: the string to repeat.
    :return: None
    """
    printf(append_chars("", n, c))

def split_file_name(path):
    """
    breaks apart a path into folder, filename and extension.
    :param path: the path to split
    :return: (folder, file name, extension)
    :rtype: tuple of str, str, str
    """
    folder, file = os.path.split(path)
    file_name, ext = os.path.splitext(file)
    if not folder:
        folder = '.'
    return folder, file_name, ext

class Loading(object):
    # how long the bar is
    bar = 40

    # current length
    n = 0

    # no idea
    length = 1

    # loading bar message
    msg=""

    def __init__(self, length, msg="", res=1):
        self.set_len(length)
        self.msg = msg
        self.start_time = time.time()
        self.res = int(res)

    def set_len(self, length):
        self.length = abs(length)

    def __call__(self):
        self.n += 1
        return self.n <= self.length

    def show(self):
        dec= self.n / self.length
        ## only show once per self.res loads
        if (dec == 1) or (not self.n % self.res):
            out = "%s|" % self.msg
            for i in range(self.bar + 0):
                out += ('=' if (i < (self.bar*dec)) else ' ')
            out += "|%.0f%%" % (100 * dec)
            
            if self.n:
                etc = (
                    (time.time() - self.start_time) *
                    (self.length - self.n) / self.n)
                n_hrs = etc // 3600
                n_minutes = (etc - (n_hrs * 3600)) // 60
                n_secs = (etc - (n_hrs * 3600) - (n_minutes * 60))
                stime = ""
                if n_hrs:
                    stime += "%dh" % int(n_hrs)
                if n_minutes:
                    stime += "%dm" % int(n_minutes)

                stime += "%ds" % int(n_secs)
                out += " ETC:%s" % stime

            printf("\x1b[2K%s\r" % out)
            sys.stdout.flush()
        if dec == 1:
            printf("\n")

def combine_tables(base, tab):
    """
    Is this the same as vstack?
    """
    if not base:
        return tab
    else:
        return vstack([base,tab])

def export_region(
        tab, colour=DEFAULT_COLOUR, scale_radius=1, region_radius=3, x_col=RA,
        y_col=DEC, wcs=1, f_name=TMP_OUT):
    """
    A handy function to convert the detections in a DS9 region file

    :param tab: Source list table with some kind of positional columns
    :type tab: Table
    :param colour: Region Colour
    :type colour: str
    :param scale_radius: Scale region radius with flux ? true/false
    :type scale_radius: int
    :param region_radius: Otherwise, use this region radius in pixels
    :type region_radius: int
    :param x_col: X column name to use
    :type x_col: str
    :param y_col: Y column name to use
    :type y_col: str
    :param wcs: Boolean if the xycols use WCS system
    :type wcs: int.
    :param f_name: Filename to output to
    :type f_name: str
    :return:
    """

    if x_col not in tab.colnames:
        x_cols= list(filter(lambda s: 'x' == s[0], tab.colnames))
        if x_cols:
            x_col = x_cols[0]
            printf("Using '%s' as x position column\n" % s_bold(x_col))
            wcs = 0

    if y_col not in tab.colnames:
        y_cols = list(filter(lambda s: 'y' == s[0], tab.colnames))
        if y_cols:
            y_col = y_cols[0]
            printf("Using '%s' as y position column\n" % s_bold(y_col))
            wcs = 0

    if "flux" in tab.colnames and scale_radius: 
        r = (-40.0 / np.log10(tab[FLUX]))
        r[r < region_radius] = region_radius
        r[np.isnan(r)] = region_radius
    else:
        r = np.ones(len(tab)) * region_radius

    prefix = "fk5;" if wcs else ""

    with open(f_name, 'w') as fp:
        fp.write("global color=%s width=2\n" % colour)
        if tab:
            for src, ri in zip(tab, r[r>0]):
                fp.write("%scircle %f %f %fi\n" % (
                    prefix, src[x_col], src[y_col], ri))
        else:
            p_error("unable to open %f\n" % f_name)

def parse_unit(raw):
    """
    Take a value with the ability to be cast into several units and parse it
    i.e. 123p -> 123 'pixels'

    Recognised units are:
    p : pixels
    s : arcsec
    m : arcmin
    d : degree

    :param raw: Raw input string to operate on
    :type raw: str
    :return: Numerical value of unit
    :rtype float
    """

    recognised = {
        'p': starbug2.PIX,
        's': starbug2.ARCSEC,
        'm': starbug2.ARCMIN,
        'd': starbug2.DEG}
    value = None
    unit = None
    if raw:
        try: 
            value = float(raw)
            unit = None
        except ValueError:
            try:
                value = float(raw[:-1])
                unit = recognised.get(raw[-1])
            except ValueError:
                p_error("unable to parse '%s'\n" % raw)
    return value, unit

def tab2array(tab, col_names=None):
    """
    Returns the contents of the table as a normal 2D numpy array
    NB: this is different from Table.asarray(), which returns an array of
    numpy.voids

    if col_names not None, return the subset of the table corresponding to
    this list

    :param tab: Table to operate on
    :type tab: astropy.Table
    :param col_names: Column names in table to include in the array
    :type col_names: list of str
    :return: Array from the table
    :rtype: numpy.ndarray
    """
    if not col_names:
        col_names=tab.col_names
    else:
        col_names=remove_duplicates(col_names)
    return np.array(tab[col_names].as_array().tolist())

def collapse_header(header):
    """
    Convert a dictionary to a Header.
    Parameters in PARAMFILES have keys longer than 8 chars
    which can cause issues in the fits format. This function turns
    those to comment cards.

    :param header: Header or dictionary to convert to collapse header
    :type header: dict, fits.Header
    :return: Collapsed Header
    :rtype fits.Header
    """
    out = fits.Header()
    for key,value in header.items():
        if len(key) > 8:
            out["comment"] = ":".join([key, str(value)])
        else: out[key] = value
    return out


def export_table(table, f_name=None, header=None):
    """
    Export table with correct dtypes

    :param table: Table to export.
    :type table: astropy.Table
    :param f_name: Filename to export to.
    :type f_name: str
    :param header: Optional header file to include in fits table
    :type header: dict, fits.Header
    :return: None
    """
    dtypes = []
    if CAT_NUM not in table.col_names:
        table = reindex(table)
    for name in table.colnames:
        if name == CAT_NUM:
            dtypes.append(str)
        elif name == "flag":
            dtypes.append(np.uint16)
        else:
            dtypes.append(table[name].dtype)
    table = fill_nan(Table(table, dtype=dtypes))

    if not f_name:
        f_name = TMP_FITS
    fits.BinTableHDU(data=table, header=header).writeto(
        f_name, overwrite=True, output_verify="fix")

def import_table(f_name, verbose=0):
    """
    Slight tweak to `astropy.table.Table.read`. This makes sure that the
    proper column dtypes are maintained

    :param f_name: Path to binary fits table file
    :type f_name: str
    :param verbose: Display verbose information
    :type verbose: boolean
    :return: Loading table
    :rtype: atrophy.Table
    """

    tab = None
    if os.path.exists(f_name):
        if os.path.splitext(f_name)[1] == FITS_EXTENSION:
            tab = fill_nan(Table.read(f_name, format="fits"))
            if not tab.meta.get(FILTER):
                if filter_string := find_filter(tab):
                    tab.meta[FILTER] = filter_string
            if verbose:
                printf("-> loaded %s (%s:%d)\n" % (
                    f_name, tab.meta.get(FILTER), len(tab)))
        else:
            p_error("Table must fits format\n")
    else:
        p_error("Unable to locate \"%s\"\n" % f_name)
    return tab

def fill_nan(table):
    """
    Fill empty values in table with nans. This is useful for tables that
    have columns that don't support nans (e.g. starbug flag). These will be
    set to zero instead

    :param table: table to operate on
    :type table: atrophy.table
    :return: Input table with masked vales filled in as nan
    :rtype: atrophy.table
    """
    for i, name in enumerate(table.col_names):
        match table.dtype[i].kind:
            case 'f': fill_val=np.nan
            case 'i' | 'u': fill_val=0
            case _: fill_val=np.nan
        if type(table[name]) == MaskedColumn:
            table[name]=table[name].filled(fill_val)
    return table

def find_col_names(tab, basename):
    """
    Find substring (basename) within the table colnames. Searches for
    substring at the beginning of the word I.E search for "flux" in
    ("flux_out","flux_err","dflux") returns as ("flux_out","flux_err")

    :param tab: Table to operate on
    :type tab: atrophy.table
    :param basename: String basename to search
    :type basename: str
    :return: List of all matching column names
    :rtype: list of str
    """
    return [
        col_name for col_name in tab.col_names
            if col_name[:len(basename)] == basename]

def combine_file_names(f_names, n_mismatch=N_MIS_MATCHES):
    """
    when matching catalogues, combines the file names into an appropriate
    combination of all the inputs.

    :param f_names: list of file names
    :type f_names: list of str
    :param n_mismatch: The number of mismatched characters it will allow
    :type n_mismatch: int
    :return: Combined filenames
    :rtype: str
    """

    trys = 0
    f_name = ""
    d_name, _, ext = split_file_name(f_names[0])
    f_names = [split_file_name(name)[1] for name in f_names]
    
    for i in range(len(f_names[0])):
        chars = [name[i] for name in f_names if len(name) > i]
        if len(set(chars)) == 1:
            f_name += chars[0]
        else: 
            f_name += "(%s)" % "".join(sorted(set(chars)))
            trys += 1
        if trys > n_mismatch:
            return None
    while ")(" in f_name:
        f_name = f_name.replace(")(","")
    return "%s/%s%s" % (d_name, f_name, ext)


def h_cascade(tables, col_names=None):
    """
    Similar use as hstack Except rather than adding a full new column,
    the inserted value is placed into the leftmost empty column

    :param tables: Table to h_cascade.
    :type tables: list of atrophy.Table
    :param col_names: List of column names to include in the stacking.
        If colnames=None, use all possible columns
    :type col_names: list of str
    :return: Single combined table
    :rtype: atrophy.Table
    """
    tab = fill_nan(hstack(tables))

    if not col_names:
        col_names = tables[0].colnames
    for name in col_names:
        cols = find_col_names(tab, name)
        if not cols:
            continue
        move = 1
        while move:
            move = 0
            for n in range(len(cols) - 1, 0, -1):
                ##everything that has a value
                curr_mask = np.invert(np.isnan(tab[cols[n]]))

                ##everything empty in left neighbouring column
                left_mask = np.isnan(tab[cols[n - 1]])

                ##cur has value and left is empty
                mask = np.logical_and(curr_mask, left_mask)

                tab[cols[n]] = MaskedColumn(tab[cols[n]])
                tab[cols[n - 1]][mask] = tab[cols[n]][mask]
                tab[cols[n]][mask] = tab[cols[n]].info.mask_val
                if sum(mask):
                    move = 1
        cols = find_col_names(tab, name)
        if cols:
            tab.rename_columns(
                cols, ["%s_%d" % (name, i + 1) for i in range(len(cols))])

    for name in tab.colnames:
        col = tab[name]

        # Use getattr to safely check for n_bad without crashing
        n_bad = getattr(col.info, 'n_bad', 0)

        if n_bad == len(col):
            tab.remove_column(name)
    return tab

def ext_names(hdu_list):
    """
    Return list of HDU extension names

    :param hdu_list: fits hdu_list to operate on
    :type hdu_list: HDUList
    :return: List of extension names
    :rtype: list of str
    """
    return list(ext.name for ext in hdu_list)

def flux2mag(flux, flux_err=None, zp=1):
    """
    Convert flux to magnitude in an arbitrary system

    :param flux: List of source flux values
    :type flux: list of floats or float
    :param flux_err: List of known flux uncertainties
    :type flux_err: list of floats or float
    :param zp: Zero point flux value
    :type zp: float
    :return: tuple of (Source magnitudes, Magnitude errors )
    :rtype: tuple (float, float)
    """

    ## sort any type issues in FLUX
    if type(flux) != np.array:
        flux = np.array(flux)
    if not flux.shape:
        flux = np.array([flux])

    # sort type issues in FLUXERR
    if flux_err is None:
        flux_err = np.zeros(len(flux))
    if type(flux_err) != np.array:
        flux_err = np.array(flux_err)
    if not flux_err.shape:
        flux_err = np.array([flux_err])

    mag = np.full(len(flux), np.nan)
    mag_err = np.full(len(flux), np.nan)

    mask_flux = (flux > 0)
    mask_f_err = (flux_err >= 0)
    mask= np.logical_and(mask_flux, mask_f_err)

    mag[mask_flux]= -2.5 * np.log10(flux[mask_flux] / zp)
    mag_err[mask] = 2.5 * np.log10(1.0 + (flux_err[mask] / flux[mask]))

    return mag, mag_err


def flux_2_ab_mag(flux, flux_err=None):
    """
    Convert flux to AB magnitudes

    :param flux: Source flux values.
    :type flux: float
    :param flux_err: Source flux error values if known.
    :type flux_err: float
    :return: Magnitude in AB system
    :rtype: float
    """
    return flux2mag(flux, flux_err, zp=3631.0)


def wget(address, f_name=None):
    """
    A really simple "implementation" of wget.

    :param address: URL to download
    :type address: str
    :param f_name: Filename to save output to
    :type f_name: str
    :return: 0 on success, 1 on failure
    :rtype int
    """
    r = requests.get(address)
    if r.status_code == REST_SUCCESS_CODE:
        f_name = f_name if f_name else os.path.basename(address)
        with open(f_name, "wb") as fp:
            for chunk in r.iter_content(chunk_size=128):
                fp.write(chunk)
        return EXIT_SUCCESS
    else:
        p_error("Unable to download \"%s\"\n" % address)
        return EXIT_FAIL


def reindex(table):
    """
    Add indexes into a table

    :param table: the table to reindex
    :type table: atrophy.Table
    :return: the reindex-ed table
    :rtype: atrophy.Table
    """

    if CAT_NUM in table.col_names:
        table.remove_column(CAT_NUM)
    column = Column(
        ["CN%d" % i for i in range(len(table))], name=CAT_NUM)
    table.add_column(column, index=0)
    return table

def colour_index(table, keys):
    """
    Allow table indexing with A-B

    :param table: table to colour index.
    :type table: atrophy.Table
    :param keys: column names to index.
    :return: A table which has only columns defined in keys.
    :rtype: atrophy.Table
    """

    out = Table()
    for key in keys:
        if key in table.col_names:
            out.add_column(table[key])
        elif '-' in key:
            a, b = key.split('-')
            out.add_column(table[a] - table[b], name=key)
    return out

def get_mj_ysr2jy_scale_factor(ext):
    """
    Find the unit scale factor to convert an image from MJy/sr to Jy
    Header file must contain the keyword "PIXAR_SR"

    :param ext: Fits extension with header file
    :type ext: PrimaryHDU, ImageHDU, BinaryTableHDU
    :return: Value of scaling factor from the header
    :rtype float
    """
    scale_factor = 1
    if ext.header.get("BUNIT") == "MJy/sr":
        if "PIXAR_SR" in ext.header:
            scale_factor = 1e6 * float(ext.header["PIXAR_SR"])
    return scale_factor

def find_filter(table):
    """
    Attempt to identify filter for a table from the metadata or column names

    :param table: Table to work on.
    :type table: astropy.table.Table
    :return: Identified filter value, otherwise None.
    :rtype str
    """
    if not (filter_string := table.meta.get(FILTER)):
        lst = (set(table.colnames) & set(starbug2.filters.keys()))
        if lst:
            filter_string = lst.pop()
            return filter_string
    return None

def get_version():
    """
    Try to determine the installed starbug version on the system

    :return: the StarBugII version string
    :rtype str
    """

    try:
        version = metadata.version("starbug2")
    except (AttributeError, TypeError, PackageNotFoundError):
        ## GitHub pytest work around for now
        version = "UNKNOWN"
    return version

def remove_duplicates(seq):
    """
     Take a sequence and rm its duplicates while preserving the order
    of the input

    :param seq: Input list to work on
    :type seq: list of <T>
    :return: A copy of the list with the duplicate elements removed
    :rtype list of <T>
    """
    seen = set()
    return [x for x in seq if not (x in seen or seen.add(x))]

def crop_hdu(hdu, x_limit=None, y_limit=None):
    """
    Crop an image with multiple extensions. Retaining the extensions

    :param hdu: A multi frame fits HDUList
    :type hdu: fits.HDUList
    :param x_limit: Pixel X bounds to crop image between
    :type x_limit: list of ?????
    :param y_limit: Pixel Y bounds to crop image
    :type y_limit: list of ????
    :return: The full HDUList that has been spatially cropped
    :rtype fits.HDUList
    """
    if x_limit is None or y_limit is None: return None

    for ext in hdu:
        if type(ext) not in (fits.PrimaryHDU, fits.ImageHDU):
            continue
        if not ext.header["NAXIS"]:
            continue
        
        ctype = ext.header.get("CTYPE")
        ext.header["CTYPE"] = "%s-SIP" % ctype

        w = WCS(ext.header, relax=False)
        ext.data = ext.data[x_limit[0]:x_limit[1], y_limit[0]:y_limit[1]]
        ext.header.update(
            w[x_limit[0]:x_limit[1], y_limit[0]:y_limit[1]].to_header())
    return hdu


if __name__ == "__main__":
    print(parse_unit(""))
    print(parse_unit("10p"))
    print(parse_unit("10 p"))
    print(parse_unit("10 D"))
    print(parse_unit("10"))
    print(parse_unit("p10"))


