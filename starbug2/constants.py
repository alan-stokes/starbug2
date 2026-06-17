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
from enum import Enum

# the filter id which we've had to adjsut the bin size to allow it to
# initilise without errors.
PROBLEMATIC_FILTER_ID = "F150W2"
PROBLEMATIC_FILTER_WARNING = (
    "Caution needed with F150W2 photometric accuracy - please check out "
    "carefully. More Info can be found in ("
    "https://github.com/alan-stokes/starbug2/issues/2)")

STARBUG_DATA_DIR: Final[str] = "STARBUG_DATDIR"
WEBBPSF_PATH_ENV_VAR: Final[str] = "WEBBPSF_PATH"
STAR_BUG_PARAMS: Final[str] = "STARBUGII PARAMETERS"
STAR_BUG_TEST_DAT_ENV: Final[str] = "STARBUG_TEST_DIR"

# default values
DEFAULT_FULL_WIDTH_HALF_MAX = 2.0
DEFAULT_PSF_FILE_NAME = "psf.fits"
DEFAULT_COLOUR: Final[str] = "green"
# how many characters we will allow by default.
N_MIS_MATCHES: Final[int] = 10

# rest success
REST_SUCCESS_CODE: Final[int] = 200

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

## SOURCE FLAGS
class SourceFlags(int, Enum):
    SRC_GOOD = 0
    SRC_BAD = 0x01
    SRC_JMP = 0x02
    ##source frame mean >5% different from median
    SRC_VAR = 0x04
    ##psf fit with fixed centroid
    SRC_FIX = 0x08
    ##source unknown (this isnt used anywhere!)
    SRC_UKN = 0x10

##DQ FLAGS
class DQFlags(int, Enum):
    DQ_DO_NOT_USE = 0x01
    DQ_SATURATED = 0x02
    DQ_JUMP_DET = 0x04

# e name common names
SCI: Final[str] = "SCI"
BGD: Final[str] = "BGD"
RES: Final[str] = "RES"

# test states
class ExitStates(int, Enum):
    EXIT_SUCCESS = 0
    EXIT_FAIL = 1
    EXIT_EARLY = 2
    EXIT_MIXED = 3

# table column enum to be used to amtch table col names
class TableColumn(str, Enum):
    """Table column names used across the pipeline."""

    CAT_NUM = "Catalogue_Number"
    RA = "RA"
    DEC = "DEC"
    FLUX = "flux"
    E_FLUX = "eflux"
    FLUX_2 = "flux_2"
    X_CENTROID = "x_centroid"
    Y_CENTROID = "y_centroid"
    X_PEAK = "x_peak"
    Y_PEAK = "y_peak"
    EE_FRACTION = "eefraction"
    RADIUS = "radius"
    AP_CORR = "apcorr"
    STD_FLUX = "stdflux"
    NUM = "NUM"
    FLAG = "flag"
    FLUX_DET = "flux_det"
    FLUX_FIT = "flux_fit"
    FLUX_ERR = "flux_err"
    OUT_FLUX = "outflux"
    X_0 = "x_0"
    Y_0 = "y_0"
    X_DET = "x_det"
    Y_DET = "y_det"
    ID = "id"
    MAG = "mag"
    MAG_UPPER = "MAG"
    ERROR_MAG = "eMAG"
    STATUS = "status"
    REC = "rec"
    PARAM = "PARAM"
    X_INIT = "x_init"
    Y_INIT = "y_init"
    XY_DEV = "xydev"
    XY_DEV_ = "_xydev"
    ERR_LOWER = "err"
    OFF = "off"
    X_FIT = "x_fit"
    Y_FIT = "y_fit"
    Q_FIT = "qfit"
    PUPIL = "pupil"
    SKY = "sky"
    SMOOTHNESS = "smoothness"
    SHARPNESS = "sharpness"
    ROUNDNESS1 = "roundness1"
    ROUNDNESS2 = "roundness2"
    RA_1 = "RA_1"
    RA_2 = "RA_2"

    # needed as the table system doenst seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system doenst seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)

## DEFAULT MATCHING COLS
MATCH_COLS: List[str] = [
    TableColumn.RA, TableColumn.DEC, TableColumn.FLAG, TableColumn.FLUX,
    TableColumn.E_FLUX, TableColumn.NUM]

# Q table col names
class QTableColNames(str, Enum):
    SUM_ERR_0 = "aperture_sum_err_0"
    SUM_0 = "aperture_sum_0"
    SUM_1 = "aperture_sum_1"

    # needed as the table system doenst seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system doenst seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)

# tag for header
class HeaderTags(str, Enum):
    FILTER_LOWER = "filter"
    FILTER = "FILTER"
    EXT = "XTENSION"
    IMAGE = "IMAGE"
    BIN_TABLE = "BINTABLE"
    OUTPUT = "OUTPUT"
    STAR_BUG = "STARBUG"
    CALIBRATION_LV = "CALIBLEVEL"
    NAXIS = "NAXIS"
    NAXIS1 = "NAXIS1"
    NAXIS2 = "NAXIS2"
    C_TYPE = "CTYPE"
    OBS = "OBSERVTN"
    VISIT = "VISIT"
    EXPOSURE = "EXPOSURE"

    # needed as the table system doenst seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system doenst seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)

# tags for image header
class ImageHeaderTags(str, Enum):
    DETECTOR = "DETECTOR"
    TELESCOPE = "TELESCOP"
    INSTRUMENT = "INSTRUME"
    BUN_IT = "BUNIT"
    PIXAR_A2 = "PIXAR_A2"
    PIXAR_SR = "PIXAR_SR"
    JWST = "JWST"
    FILTER = "FILTER"

    # needed as the table system doenst seem to handle enums properly
    def __str__(self) -> str:
        return self.value

    # needed as the table system doenst seem to handle enums properly
    def __format__(self, format_spec: str) -> str:
        return self.value.__format__(format_spec)

# tag used for param file.
VERBOSE_TAG: Final[str] = "VERBOSE"

# mode labels.
class Modes(str, Enum):
    DETECTION = "DETECTION"
    BACKGROUND = "BACKGROUND"
    APP_HOT = "APPHOT"
    PSFP_HOT = "PSFPHOT"
    MATCH_OUTPUTS = "MATCHOUTPUTS"
    CLEAR = "CLEAR"


## HASHDEFS
STAR_BUG_MIRI: Final[int] = 1
NIRCAM: Final[int] = 2
NIRCAM_STRING: Final[str] = "NIRCAM"

class DetectorLengths(int, Enum):
    NULL = 0
    LONG = 1
    SHORT = 2

# enum unit
class Units(int, Enum):
    PIX = 0
    ARCSEC = 1
    ARCMIN = 2
    DEG = 3

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
    Modes.DETECTION :
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
    Modes.BACKGROUND:
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
    Modes.APP_HOT:
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
    Modes.PSFP_HOT:
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
    Modes.MATCH_OUTPUTS:
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

# Named Constant Template for Default Parameter Files
DEFAULT_PARAM_TEMPLATE: Final[str] = """## STARBUG CONFIG FILE
# Generated with starbug2-v{version_str}
PARAM       =  STARBUGII PARAMETERS     // COMMENT

## GENERIC
// (0:false 1:true)
VERBOSE     = {VERBOSE}

// Directory or filename to output to 
OUTPUT      = {OUTPUT}

// If using a non standard HDU name, name it here (str or int)
HDUNAME     = {HDUNAME}

// Set a custom filter for the image
FILTER      = {FILTER}

## DETECTION 
// Custom FWHM for image (-1 to use WEBBPSF)
FWHM        = {FWHM}

// Number of sigma above the median to clip out as background
SIGSKY      = {SIGSKY}

// Source value minimum N sigma above background
SIGSRC      = {SIGSRC}

// Run background2D step (usually finds more sources but takes time)
DOBGD2D     = {DOBGD2D}

// Run convolution step (usually finds more sources)
DOCONVL     = {DOCONVL}

// Run source cleaning after detection (removes likely contaminants)
CLEANSRC    = {CLEANSRC}

// Lower limit of source sharpness (0 is not sharp)
SHARP_LO    = {SHARP_LO}

// Upper limit of source sharpness (1 is sharp)
SHARP_HI    = {SHARP_HI}

// Limit of source roundness1 (|roundness|>>0 is less round)
ROUND1_HI   = {ROUND1_HI} 

// Limit of source roundness2 (|roundness|>>0 is less round)
ROUND2_HI   = {ROUND2_HI}

// Lower limit on source smoothness (0 is not smooth)
SMOOTH_LO   = {SMOOTH_LO}

// Upper limit on source smoothness (1 is smooth)
SMOOTH_HI   = {SMOOTH_HI}

// Radius (pix) of ricker wavelet 
RICKER_R    = {RICKER_R}

## APERTURE PHOTOMETRY
// Radius in number of pixels
APPHOT_R    = {APPHOT_R}

// Fraction encircled energy (mutually exclusive with APPHOT_R)
ENCENERGY   = {ENCENERGY} 

// Sky annulus inner radius
SKY_RIN     = {SKY_RIN} 

// Sky annulus outer radius
SKY_ROUT    = {SKY_ROUT}  

// Aperture correction file. See full manual for details
APCORR_FILE = {APCORR_FILE}

## BACKGROUND ESTIMATION
// Aperture masking fixed radius (if zero, starbug will scale radii)
BGD_R       = {BGD_R} 

// Aperture mask radius profile scaling factor
PROF_SCALE  = {PROF_SCALE}

// Aperture mask radius profile slope
PROF_SLOPE  = {PROF_SLOPE} 

// Background estimation kernel size (pix)
BOX_SIZE    = {BOX_SIZE}

// Output region file to check the aperture mask radii
BGD_CHECKFILE = {BGD_CHECKFILE}

## PHOTOMETRY
// Detection file to use instead of detecting
AP_FILE     = {AP_FILE}

// Background estimation file
BGD_FILE    = {BGD_FILE}

// Non default PSF file
PSF_FILE    = {PSF_FILE}

// When loading an AP_FILE, do you want to use WCS or xy values (if available)
USE_WCS     = {USE_WCS}

// Zero point (mag) to add to the magnitude columns 
ZP_MAG      = {ZP_MAG} 

// Minimum distance for grouping (pixels) between two sources
CRIT_SEP    = {CRIT_SEP}

// Force centroid position (1) or allow psf fitting to fit position too (0)
FORCE_POS   = {FORCE_POS}

// Maximum deviation from initial guess centroid position
MAX_XYDEV   = {MAX_XYDEV}

// Set fit size of psf (>0) or -1 to take PSF file dimensions
PSF_SIZE    = {PSF_SIZE}

// Generate a residual image
GEN_RESIDUAL = {GEN_RESIDUAL}

## SOURCE STATS
// Run crowding metric calculation (execution time scales N^2)
CALC_CROWD  = {CALC_CROWD}

## CATALOGUE MATCHING
// Matching separation threshold in units arcsec
MATCH_THRESH = {MATCH_THRESH}

// EXTRA columns to include in output matched table i.e sharpness
MATCH_COLS   = {MATCH_COLS}

// Keep sources that appear in NUM >= NEXP_THRESH (if -1 keep everything)
NEXP_THRESH  = {NEXP_THRESH}

// Bridge --band matching NIRCam and MIRI catalogues by ensuring NIRCam 
// catalogue has a match in BRIDGE_COL
BRIDGE_COL   = {BRIDGE_COL}

## ARTIFICIAL STAR TESTS
// Number of artificial star tests
NTESTS      = {NTESTS}

// Number of stars per artificial test
NSTARS      = {NSTARS}

// Number of pixels to crop around artificial star
SUBIMAGE    = {SUBIMAGE}

// Bright limit of test magnitude
MAX_MAG     = {MAX_MAG}

// Faint limit of test magnitude
MIN_MAG     = {MIN_MAG}

// Output AST result as image with this filename
PLOTAST     = {PLOTAST}

## MISC EXTRAS
// DS9 region colour
REGION_COL  = {REGION_COL}

// Scale region to flux if possible
REGION_SCAL = {REGION_SCAL}

// Region radius default
REGION_RAD  = {REGION_RAD}

// X column name to use for region
REGION_XCOL = {REGION_XCOL}

// Y column name to use for region
REGION_YCOL = {REGION_YCOL}

// If X/Y column names correspond to WCS values
REGION_WCS  = {REGION_WCS}

// detector name used within psf generation
DET_NAME = {DET_NAME}

// region table file name for generating regions
REGION_TAB = {REGION_TAB}

## ADDITIONAL UNSETTABLE / DERIVED PARAMETERS
// Custom analytical parameter grouping identifier string
PARAM_TAG   = {PARAM}

## PIPELINE SWITCHES (BYPASS COMMAND LINE FLAGS)
// Run aperture photometry stage (0:false 1:true)
RUN_APPHOT  = {RUN_APPHOT}

// Run background estimation stage (0:false 1:true)
RUN_BGD_EST = {RUN_BGD_EST}

// Run star detection stage (0:false 1:true)
RUN_DETECT  = {RUN_DETECT}

// Run source geometry analysis stage (0:false 1:true)
RUN_GEOM    = {RUN_GEOM}

// Run catalogue matching stage (0:false 1:true)
RUN_MATCH   = {RUN_MATCH}

// Run PSF photometry routine (0:false 1:true)
RUN_PSFPHOT = {RUN_PSFPHOT}

// Run background subtraction stage (0:false 1:true)
RUN_BGDSUB  = {RUN_BGDSUB}

## SYSTEM & EXECUTION CONTROL
// Number of processor cores to use for calculation
NCORES      = {NCORES}

// Find files automatically (0:false 1:true)
FIND_FILE   = {FIND_FILE}

## EXTRA RUN GENERATION COMMAND SWITCHES
// Execute JWST data initialization steps (0:false 1:true)
INIT_JWST   = {INIT_JWST}

// Trigger PSF generation logic (0:false 1:true)
GEN_PSF     = {GEN_PSF}

// Trigger automation run generation scripts (0:false 1:true)
GEN_RUN     = {GEN_RUN}

// Filename target string to output generated region file
GEN_REGION  = {GEN_REGION}

## ADVANCED ARTIFICIAL STAR TEST CONTROLS
// Recover previous artificial star test state (0:false 1:true)
AST_RECOVER = {AST_RECOVER}

// Save frequency of progress during artificial star tests
AST_AUTOSAVE = {AST_AUTOSAVE}

// Disable background logic during artificial star tests (0:false 1:true)
AST_NO_BGD  = {AST_NO_BGD}

// Disable PSF photometry during artificial star tests (0:false 1:true)
AST_NO_PSF  = {AST_NO_PSF}

## CATALOGUE MATCHING MODE SWITCHES
// Process matched catalogue across multiple bands (0:false 1:true)
MATCH_BAND  = {MATCH_BAND}

// Run matching in cascade execution sequence (0:false 1:true)
MATCH_CASCADE = {MATCH_CASCADE}

// Use dither offsets during matching calculations (0:false 1:true)
MATCH_DITHER = {MATCH_DITHER}

// Require exact row criteria matches across tables (0:false 1:true)
MATCH_EXACT = {MATCH_EXACT}

// Force a full evaluation run across matching pipelines (0:false 1:true)
MATCH_FULL  = {MATCH_FULL}

// Use generic operating specifications for matching (0:false 1:true)
MATCH_GENERIC = {MATCH_GENERIC}

// Error column label name target string for evaluation
MATCH_ERR_COL = {MATCH_ERR_COL}

// Expression filter string for processing match masks
MATCH_MASK_EVAL = {MATCH_MASK_EVAL}

## DIAGNOSTIC PLOTTING MODULE SWITCHES
// Enable interactive test mode plotting panels (0:false 1:true)
PLOT_TEST   = {PLOT_TEST}

// Dark frame visualization profile layout (0:false 1:true)
PLOT_DARK   = {PLOT_DARK}

// Visual parameter validation target inspection parameter key
PLOT_INSPECT = {PLOT_INSPECT}

// Target design stylesheet pattern configuration profile name
PLOT_STYLE  = {PLOT_STYLE}
"""