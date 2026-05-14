# noinspection SpellCheckingInspection

STARBUG_DATA_DIR = "STARBUG_DATDIR"
WEBBPSF_PATH_ENV_VAR = "WEBBPSF_PATH"

# url to docs
URL_DOCS = (
    "https://raw.githubusercontent.com/conornally/starbug2/"
    "refs/heads/main/docs/source/_static/images/starbug.png")
READ_THE_DOCS_URL = "https://starbug2.readthedocs.io/en/latest/"

# fit urls
JWST_MIRI_APCORR_0010_FITS_URL = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_miri_apcorr_0010.fits"
)
JWST_NIRCAM_APCORR_0004_FITS_URL = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_nircam_apcorr_0004.fits"
)

# abvega offset urls
JWST_MIRI_ABVEGA_OFFSET_URL = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_miri_abvegaoffset_0001.asdf"
)
JWST_NIRCAM_ABVEGA_OFFSET_URL = (
    "https://jwst-crds.stsci.edu/unchecked_get/references/jwst/"
    "jwst_nircam_abvegaoffset_0002.asdf"
)


# problematic paths
PLOT_MAIN_TABLE_PATH = (
    "/home/conor/sci/proj/ngc6822/overview/dat/ngc6822.fits")
MASK_MAIN_TABLE_PATH = (
    "/home/conor/sci/proj/ngc6822/paper1/dat/ngc6822.fits")

# paths to temp files.
TMP_OUT = "/tmp/out.reg"
TMP_FITS = "/tmp/starbug.fits"

# the fits file extension
FITS_EXTENSION = ".fits"
FILE_NAME = "FILENAME"

# HDU extension names
DQ = "DQ"
AREA = "AREA"
WHT = "WHT"
ERR = "ERR"

# file types
AP_FILE = "AP_FILE"
BGD_FILE = "BGD_FILE"

# init parameters
DET_NAME = "DET_NAME"
PSF_SIZE = "PSF_SIZE"
REGION_COL = "REGION_COL"
REGION_SCAL = "REGION_SCAL"
REGION_RAD = "REGION_RAD"
REGION_X_COL = "REGION_XCOL"
REGION_Y_COL = "REGION_YCOL"
REGION_WCS = "REGION_WCS"

# colours
DEFAULT_COLOUR = "green"

## SOURCE FLAGS
SRC_GOOD = 0
SRC_BAD = 0x01
SRC_JMP = 0x02
SRC_VAR = 0x04 ##source frame mean >5% different from median
SRC_FIX = 0x08 ##psf fit with fixed centroid
SRC_UKN = 0x10 ##source unknown

##DQ FLAGS
DQ_DO_NOT_USE = 0x01
DQ_SATURATED = 0x02
DQ_JUMP_DET = 0x04


# some binary values.
VERBOSE = 0x01
KILLPROC = 0x02
STOPPROC = 0x04
SHOWHELP = 0x08

DODETECT = 0x100
DOBGDEST = 0x200
DOPHOTOM = 0x400
FINDFILE = 0x800

DOARTIFL = 0x1000
DOMATCH = 0x2000
DOAPPHOT= 0x4000
DOBGDSUB = 0x8000
DOGEOM = 0x10000

GENRATPSF = 0x100000
GENRATRUN = 0x200000
GENRATREG = 0x400000
INITSB = 0x800000
UPDATEPRM = 0x1000000
DODEBUG = 0x2000000
CALCINSTZP = 0x4000000
APPLYZP = 0x8000000

# option names
HDU_NAME = "HDUNAME"

# e name common names
SCI = "SCI"
BGD = "BGD"
RES = "RES"

# test states
EXIT_SUCCESS = 0
EXIT_FAIL = 1
EXIT_EARLY = 2
EXIT_MIXED = 3

# rest success
REST_SUCCESS_CODE = 200

# tag used table col names
CAT_NUM = "Catalogue_Number"
RA = "RA"
DEC = "DEC"
FLUX = "flux"
E_FLUX = "eflux"
X_CENTROID = "xcentroid"
Y_CENTROID = "ycentroid"
X_PEAK = "x_peak"
Y_PEAK = "y_peak"
EE_FRACTION = "eefraction"
RADIUS = "radius"
AP_CORR = "apcorr"

## DEFAULT MATCHING COLS
match_cols = [RA, DEC, "flag", FLUX, "eflux", "NUM"]

# tag for header
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

# tags for image header
DETECTOR = "DETECTOR"
TELESCOPE = "TELESCOP"
INSTRUMENT = "INSTRUME"
BUN_IT = "BUNIT"
PIXAR_A2 = "PIXAR_A2"
PIXAR_SR = "PIXAR_SR"
JWST ="JWST"

# tag used for param file.
PARAM_FILE_TAG = "PARAMFILE"
REGION_TAB = "REGION_TAB"
VERBOSE_TAG = "VERBOSE"

# mode labels.
DETECTION = "DETECTION"
BACKGROUND = "BACKGROUND"
APP_HOT = "APPHOT"
PSFP_HOT = "PSFPHOT"
MATCH_OUTPUTS = "MATCHOUTPUTS"

# options
N_CORES = "NCORES"
FWHM = "FWHM"
USE_WCS = "USE_WCS"
CRIT_SEP = "CRIT_SEP"
FORCE_POS = "FORCE_POS"
MAX_XY_DEV = "MAX_XYDEV"
CALC_CROWD = "CALC_CROWD"
APCORR_FILE = "APCORR_FILE"
APPHOT_R = "APPHOT_R"
ENCENERGY = "ENCENERGY"
SKY_RIN = "SKY_RIN"
SKY_ROUT = "SKY_ROUT"
SIGSKY = "SIGSKY"
ZP_MAG = "ZP_MAG"
CLEANSRC = "CLEANSRC"
QUIETMODE = "QUIETMODE"
BOX_SIZE = "BOX_SIZE"
BGD_R = "BGD_R"
PROF_SCALE = "PROF_SCALE"
PROF_SLOPE = "PROF_SLOPE"
BGD_CHECKFILE = "BGD_CHECKFILE"
PSF_FILE = "PSF_FILE"
GEN_RESIDUAL = "GEN_RESIDUAL"

# match options
MATCH_THRESH = "MATCH_THRESH"

# match params
NEXP_THRESH = "NEXP_THRESH"
ZP_MAG = "ZP_MAG"


## HASHDEFS
MIRI = 1
NIRCAM = 2
NIRCAM_STRING = "NIRCAM"

NULL = 0
LONG = 1
SHORT = 2

# enum unit
PIX = 0
ARCSEC = 1
ARCMIN = 2
DEG = 3


# how many characters we will allow by default.
N_MIS_MATCHES = 10

# text based logo (using raw string to bypass escape characters)
LOGO = r"""
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