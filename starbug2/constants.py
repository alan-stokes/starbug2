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

# noinspection SpellCheckingInspection
from typing import List, Final

STARBUG_DATA_DIR: Final[str] = "STARBUG_DATDIR"
WEBBPSF_PATH_ENV_VAR: Final[str] = "WEBBPSF_PATH"
STAR_BUG_PARAMS: Final[str] = "STARBUGII PARAMETERS"
STAR_BUG_TEST_DAT_ENV: Final[str] = "STARBUG_TEST_DIR"

# default value for full width half max when nothing sets it
DEFAULT_FULL_WIDTH_HALF_MAX = 2.0
DEFAULT_PSF_FILE_NAME = "psf.fits"

# url to docs
URL_DOCS: Final[str] = (
    "https://raw.githubusercontent.com/conornally/starbug2/"
    "refs/heads/main/docs/source/_static/images/starbug.png")
READ_THE_DOCS_URL: Final[str] = "https://starbug2.readthedocs.io/en/latest/"

# fit urls
JWST_MIRI_APCORR_0010_FITS_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_miri_apcorr_0010.fits"
)
JWST_NIRCAM_APCORR_0004_FITS_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_nircam_apcorr_0004.fits"
)

# abvega offset urls
JWST_MIRI_ABVEGA_OFFSET_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_miri_abvegaoffset_0001.asdf"
)
JWST_NIRCAM_ABVEGA_OFFSET_URL: Final[str] = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_nircam_abvegaoffset_0002.asdf"
)


# problematic paths
PLOT_MAIN_TABLE_PATH: Final[str] = (
    "/home/conor/sci/proj/ngc6822/overview/dat/ngc6822.fits")
MASK_MAIN_TABLE_PATH: Final[str] = (
    "/home/conor/sci/proj/ngc6822/paper1/dat/ngc6822.fits")

# paths to temp files.
TMP_OUT: Final[str] = "/tmp/out.reg"
TMP_FITS: Final[str] = "/tmp/starbug.fits"

# the fits file extension
FITS_EXTENSION: Final[str] = ".fits"
FILE_NAME: Final[str] = "FILENAME"

# HDU extension names
DQ: Final[str] = "DQ"
AREA: Final[str] = "AREA"
WHT: Final[str] = "WHT"
ERR: Final[str] = "ERR"

# file types
AP_FILE: Final[str] = "AP_FILE"
BGD_FILE: Final[str] = "BGD_FILE"
PSF_FILE: Final[str] = "PSF_FILE"

# colours
DEFAULT_COLOUR: Final[str] = "green"

## SOURCE FLAGS
SRC_GOOD: Final[int] = 0
SRC_BAD: Final[int] = 0x01
SRC_JMP: Final[int] = 0x02
##source frame mean >5% different from median
SRC_VAR: Final[int] = 0x04
##psf fit with fixed centroid
SRC_FIX: Final[int] = 0x08
##source unknown (this isnt used anywhere!)
SRC_UKN: Final[int] = 0x10

##DQ FLAGS
DQ_DO_NOT_USE: Final[int] = 0x01
DQ_SATURATED: Final[int] = 0x02
DQ_JUMP_DET: Final[int] = 0x04

# e name common names
SCI: Final[str] = "SCI"
BGD: Final[str] = "BGD"
RES: Final[str] = "RES"

# test states
EXIT_SUCCESS: Final[int] = 0
EXIT_FAIL: Final[int] = 1
EXIT_EARLY: Final[int] = 2
EXIT_MIXED: Final[int] = 3

# rest success
REST_SUCCESS_CODE: Final[int] = 200

# tag used table col names
CAT_NUM: Final[str] = "Catalogue_Number"
RA: Final[str] = "RA"
DEC: Final[str] = "DEC"
FLUX: Final[str] = "flux"
E_FLUX: Final[str] = "eflux"
FLUX_2: Final[str] = "flux_2"
X_CENTROID: Final[str] = "xcentroid"
Y_CENTROID: Final[str] = "ycentroid"
X_PEAK: Final[str] = "x_peak"
Y_PEAK: Final[str] = "y_peak"
EE_FRACTION: Final[str] = "eefraction"
RADIUS: Final[str] = "radius"
AP_CORR: Final[str] = "apcorr"
STD_FLUX: Final[str] = "stdflux"
NUM: Final[str] = "NUM"
FLAG: Final[str] = "flag"
FLUX_DET: Final[str] = "flux_det"
FLUX_FIT: Final[str] = "flux_fit"
FLUX_ERR: Final[str] = "flux_err"
OUT_FLUX: Final[str] = "outflux"
X_0: Final[str] = "x_0"
Y_0: Final[str] = "y_0"
X_DET: Final[str] = "x_det"
Y_DET: Final[str] = "y_det"
ID: Final[str] = "id"
MAG: Final[str] = "mag"
STATUS: Final[str] = "status"
REC: Final[str] = "rec"
PARAM: Final[str] = "PARAM"
X_INIT: Final[str] = "x_init"
Y_INIT: Final[str] = "y_init"
XY_DEV: Final[str] = "xydev"
XY_DEV_: Final[str] = "_xydev"
ERR_LOWER: Final[str] = "err"
OFF: Final[str] = "off"
X_FIT: Final[str] = "x_fit"
Y_FIT: Final[str] = "y_fit"
Q_FIT: Final[str] = "qfit"
PUPIL: Final[str] = "pupil"
SKY: Final[str] = "sky"
SMOOTHNESS: Final[str] = "smoothness"

# Q table col names
SUM_ERR_0: Final[str] = "aperture_sum_err_0"
SUM_0: Final[str] = "aperture_sum_0"
SUM_1: Final[str] = "aperture_sum_1"

## DEFAULT MATCHING COLS
MATCH_COLS: List[str] = [RA, DEC, FLAG, FLUX, E_FLUX, NUM]

# tag for header
FILTER_LOWER: Final[str] = "filter"
FILTER: Final[str] = "FILTER"
EXT: Final[str] = "XTENSION"
IMAGE: Final[str] = "IMAGE"
BIN_TABLE: Final[str] = "BINTABLE"
OUTPUT: Final[str] = "OUTPUT"
STAR_BUG: Final[str] = "STARBUG"
CALIBRATION_LV: Final[str] = "CALIBLEVEL"
NAXIS: Final[str] = "NAXIS"
NAXIS1: Final[str] = "NAXIS1"
NAXIS2: Final[str] = "NAXIS2"
C_TYPE: Final[str] = "CTYPE"

# tags for image header
DETECTOR: Final[str] = "DETECTOR"
TELESCOPE: Final[str] = "TELESCOP"
INSTRUMENT: Final[str] = "INSTRUME"
BUN_IT: Final[str] = "BUNIT"
PIXAR_A2: Final[str] = "PIXAR_A2"
PIXAR_SR: Final[str] = "PIXAR_SR"
JWST: Final[str] = "JWST"

# tag used for param file.
VERBOSE_TAG: Final[str] = "VERBOSE"

# mode labels.
DETECTION: Final[str] = "DETECTION"
BACKGROUND: Final[str] = "BACKGROUND"
APP_HOT: Final[str] = "APPHOT"
PSFP_HOT: Final[str] = "PSFPHOT"
MATCH_OUTPUTS: Final[str] = "MATCHOUTPUTS"
CLEAR: Final[str] = "CLEAR"

#info tags / keys for catalogue fields.
OBS: Final[str] = "OBSERVTN"
VISIT: Final[str] = "VISIT"
EXPOSURE: Final[str] = "EXPOSURE"


## HASHDEFS
STAR_BUG_MIRI: Final[int] = 1
NIRCAM: Final[int] = 2
NIRCAM_STRING: Final[str] = "NIRCAM"

NULL: Final[int] = 0
LONG: Final[int] = 1
SHORT: Final[int] = 2

# enum unit
PIX: Final[int] = 0
ARCSEC: Final[int] = 1
ARCMIN: Final[int] = 2
DEG: Final[int] = 3


# how many characters we will allow by default.
N_MIS_MATCHES: Final[int] = 10

# text based logo (using raw string to bypass escape characters)
LOGO: Final[str] = r"""
                   *          *  __  *  __   - * --   - 
 STARBUGII               *      / ___ /    \   --  -      - 
 ---------                 *___---.    .___/  -   --   -
 JWST photometry in       ./=== \  \.     \      * 
 complex crowded fields   | (O)  |  |     |           *
                           \._._/ ./    _(\)   *   
 conor.nally@ed.ac.uk     /   ~--\ ----~   \      *
 alan.stokes@stfc.ac.uk  ---      ___       ---      
 > %s
"""

# dictionary of help strings for specific modes (
# DETECTION, BACKGROUND, APPHOT, PSFPHOT, MATCHOUTPUTS).
HELP_STRINGS = {
    DETECTION :
        """
            Source Detection
            ----------------
        
            This routine locates point sources in an image. The input is a
            FITS image and the output is a FITS table, containing a list of
            point source locations, their geometric properties and 
            flux/magnitude measurements as calculated by aperture photometry. 
            The output file will have the suffix "-ap", note this is the same 
            as the output for the aperture photometry routine.
        
            To run this routine, use the core command:
        
                $~ starbug2 -D image.fits
        
            Alter the parameter file options under "DETECTION" to tune the 
            performance of starbug2. Two of the key parameters are:
            
                - SIGSKY : Set the background level of the image
                - SIGSRC : Set the detection threshold of the sources
        
            Full documentation is at https://starbug2.readthedocs.io
        """,
    BACKGROUND:
        """
            Diffuse Background Estimations
            ------------------------------
        
            This routine estimates the "dusty" emissions in an image, given
            a source list. It is used to subtract from the image, thus removing
            the flux contribution on a source brightness from the dusty 
            environment.
            
            The routine requires a list of sources to be generated (by source 
            detection) or loaded with [-d sourcelist.fits] and requires a FITS
            image to work on. The routine will ouput a FITS image, with the 
            same dimensions and spatial coverage as the input image, with the 
            suffix "-bgd". This background image can be used in the photometry
             later.
        
            To run the routine, use the core command:
        
                $~ starbug2 -B -d sourcelist.fits image.fits
        
            Alter the parameter file options under "BACKGROUND ESTIMATION" 
            to tune the performance of starbug2. Two key parameters are:
        
                - BGD_R    : Set a fixed aperture mask radius around each 
                             source
                - BOX_SIZE : Set the estimation resolution (larger will be 
                             more blurred)
        
            Full documentation is at https://starbug2.readthedocs.io
        """,
    APP_HOT:
        """
            Aperture Photometry
            -------------------
        
            This routine conducts aperture photometry on an image given a list
            of sources. It requires a FITS image to run on and a FITS table 
            source list with either RA/DEC columns, or x/y_centroid or x/y_0
            columns. The routine outputs a table with the suffix "-ap". Note
            this filename is the same as the source detection routine because 
            aperture photometry is automatically run at the end of the source 
            detection step. The output table contains 2flux/magnitude 
            information on every source
        
            To run this routine, use the core command:
        
                $~ starbug2 -A -d sourcelist.fits image.fits
        
            Alter the parameter file options under "APERTURE PHOTOMETRY" to 
            tune the performance of starbug2. Three key parameters are:
        
                - APPHOT_R : Set the aperture radius for photometry (in pixels)
                - SKY_RIN  : Set the inner sky annulus radius (in pixels)
                - SKY_ROUT : Set the outer sky annulus radius (in pixels)
        
            Full documentation is at https://starbug2.readthedocs.io
        """,
    PSFP_HOT:
        """
            PSF Photometry
            --------------
        
            This routine conducts PSF fitting photometry on an image given
            a list of sources. Its requires a FITS image to run on and a FITS
            table sourcelist with either RA/DEC columns, or x/y_centroid or 
            x/y_0 columns. The routine outputs a table with the suffix "-psf". 
            The output table contains 2flux/magnitude information on every 
            source.
            
            To run this routine, use the core command:
                
                $~ starbug2 -P -d sourcelist image.fits 
        
            Alter the parameter file options under "PHOTOMETRY" to tune the
            performance of starbug2. Two key parameters are:
        
                - FORCE_POS    : Hold the cetroid positions of source fixed 
                                 (forced photometry)
                - GEN_RESIDUAL : Generate a residual image from all the fit 
                                 source
        
            Full documentation is at https://starbug2.readthedocs.io
        """,
    MATCH_OUTPUTS:
        """
            Match Outputs
            -------------
        
            This option is set if the user wishes to combine all the output
            catalogues from starbug together. It would be used in the case 
            that a routine is being ran on a list of images (either in series 
            or parallel) and the final catalogues should all be combined into 
            a single source list. It outputs two files, one with the suffix 
            "full" and another with "match". The first is all columns from all 
            table preserved into a single large catalogue, the second averages 
            all the similar columns into a reduced table.
        
            To run this routine, use the core code:
                
                $~ starbug2 -DM image1.fits image2.fits image3.fits ...
        
            Alter the parameter file options under "CATALOGUE MATCHING" to 
            tune the performance of starbug2. Two key parameters are:
        
                - MATCH_THRESH : Set the separation threshold (arcsec) to match 
                                 two sources
                - NEXP_THRESH  : Set the minimum number of catalogues a source 
                                 must be present in
        """,
}