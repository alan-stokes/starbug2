import os
from parse import parse
from starbug2.utils import printf,p_error,get_version

# noinspection SpellCheckingInspection
default = """## STARBUG CONFIG FILE
# Generated with starbug2-v%s
PARAM       =  STARBUGII PARAMETERS     // COMMENT

## GENERIC
// (0:false 1:true)
VERBOSE     = 0

// Directory or filename to output to 
OUTPUT      = 

// If using a non standard HDU name, name it here (str or int)
HDUNAME     = SCI

// Set a custom filter for the image
FILTER      = 

## DETECTION 
// Custom FWHM for image (-1 to use WEBBPSF)
FWHM        = -1

// Number of sigma above the median to clip out as background
SIGSKY      = 2.0

// Source value minimum N sigma above background
SIGSRC      = 5.0

// Run background2D step (usually finds more sources but takes time)
DOBGD2D     = 1

// Run convolution step (usually finds more sources)
DOCONVL     = 1

// Run source cleaning after detection (removes likely contaminants)
CLEANSRC    = 1

// Lower limit of source sharpness (0 is not sharp)
SHARP_LO    = 0.4

// Upper limit of source sharpness (1 is sharp)
SHARP_HI    = 0.9

// Limit of source roundness1 (|roundness|>>0 is less round)
ROUND1_HI   = 1.0 

// Limit of source roundness2 (|roundness|>>0 is less round)
ROUND2_HI   = 1.0

// Lower limit on source smoothness (0 is not smooth)
SMOOTH_LO   = 0.0

// Upper limit on source smoothness (1 is smooth)
SMOOTH_HI   = 1.0

// Radius (pix) of ricker wavelet 
RICKER_R    = 1.0

## APERTURE PHOTOMETRY
// Radius in number of pixels
APPHOT_R    = 1.5

// Fraction encircled energy (mutually exclusive with APPHOT_R)
ENCENERGY   = -1 

// Sky annulus inner radius
SKY_RIN     = 3.0 

// Sky annulus outer radius
SKY_ROUT    = 4.5  

// Aperture correction file. See full manual for details
APCORR_FILE = 

## BACKGROUND ESTIMATION
// Aperture masking fixed radius (if zero, starbug will scale radii)
BGD_R       = 0 

// Aperture mask radius profile scaling factor
PROF_SCALE  = 1.0

// Aperture mask radius profile slope
PROF_SLOPE  = 0.5 

// Background estimation kernel size (pix)
BOX_SIZE    = 2

// Output region file to check the aperture mask radii
BGD_CHECKFILE = 

## PHOTOMETRY
// Detection file to use instead of detecting
AP_FILE     = 

// Background estimation file
BGD_FILE    = 

// Non default PSF file
PSF_FILE    = 

// When loading an AP_FILE, do you want to use WCS or xy values (if available)
USE_WCS     = 1

// Zero point (mag) to add to the magnitude columns 
ZP_MAG      = 8.9 

// Minimum distance for grouping (pixels) between two sources
CRIT_SEP    = 1.0

// Force centroid position (1) or allow psf fitting to fit position too (0)
FORCE_POS   = 0

// If allowed to fit position, max separation (arcsec) from source list 
centroid
DPOS_THRESH = -1

// Maximum deviation from initial guess centroid position
MAX_XYDEV   = 3.0

// Set fit size of psf (>0) or -1 to take PSF file dimensions
PSF_SIZE    = -1

// Generate a residual image
GEN_RESIDUAL = 0

## SOURCE STATS
// Run crowding metric calculation (execution time scales N^2)
CALC_CROWD  = 1

## CATALOGUE MATCHING
// Matching separation threshold in units arcsec
MATCH_THRESH = 0.1

// EXTRA columns to include in output matched table i.e sharpness
MATCH_COLS   = 

// Keep sources that appear in NUM >= NEXP_THRESH (if -1 keep everything)
NEXP_THRESH  = -1

// Remove sources with SN ratio < SN_THRESH before matching 
(default -1 to not apply this cut)
SN_THRESH    = -1

// Bridge --band matching NIRCam and MIRI catalogues by ensuring NIRCam 
catalogue has a match in BRIDGE_COL
BRIDGE_COL   = 

## ARTIFICIAL STAR TESTS
// Number of artificial star tests
NTESTS      = 100

// Number of stars per artificial test
NSTARS      = 10

// Number of pixels to crop around artificial star
SUBIMAGE    = 500

// Bright limit of test magnitude
MAX_MAG     = 18.0

// Faint limit of test magnitude
MIN_MAG     = 28.0

// Output AST result as image with this filename
PLOTAST     = 

## MISC EXTRAS
// DS9 region colour
REGION_COL  = green

// Scale region to flux if possible
REGION_SCAL = 1

// Region radius default
REGION_RAD  = 3

// X column name to use for region
REGION_XCOL = RA

// Y column name to use for region
REGION_YCOL = DEC

// If X/Y column names correspond to WCS values
REGION_WCS  = 1
""" % get_version()

def parse_param(line):
    """
    Parse a parameter line
    """
    param={}
    if line and line[0] not in "# \t\n":
        if "//" in line:
            key, value, _ = parse("{}={}//{}", line)
        else:
            key, value = parse("{}={}",line)
        key = key.strip().rstrip()
        value = value.strip().rstrip()
        try:
            if '.' in value:
                value = float(value)
            else:
                value = int(value)
        except (ValueError, AttributeError, TypeError):
            pass

        ## Special case values
        if key in ("OUTPUT", "AP_FILE", "BGD_FILE", "PSF_FILE"):
            value = os.path.expandvars(value)
        param[key] = value
    return param



def load_default_params():
    config = {}
    for line in default.split('\n'):
        config.update(parse_param(line))
    return config

def load_params(f_name):
    """
    Convert a parameter file into a dictionary of options

    :param f_name: path/to/file.param
    :type f_name: str
    :return: dictionary of options
    :rtype: dict of string, string
    """
    config = {}
    if f_name is None:
        config = load_default_params()
    elif os.path.exists(f_name):
        with open(f_name, "r") as fp:
            for line in fp.readlines():
                config.update(parse_param(line))
    else:
        p_error("config file \"%s\" does not exist\n" % f_name)
    return config

def local_param():
    """
    reads a local param file.
    :return: None
    """
    with open("starbug.param", "w") as fp:
        fp.write(default)

def update_param_file(f_name):
    """
    When the local parameter file is from an older version, add or remove the
    new or obsolete keys

    :param f_name: local file to update
    :type f_name: str
    :return: None
    """
    default_param = load_default_params()
    current_param = load_params(f_name)

    if os.path.exists(f_name):
        printf("Updating \"%s\"\n" % f_name)
        fpi = open(f_name, 'r')
        fpo = open("/tmp/starbug.param",'w')

        add_keys = set(default_param.keys()) - set(current_param.keys())
        del_keys = set(current_param.keys()) - set(default_param.keys())
        if add_keys:
            printf("-> adding: %s  \n"%(', '.join(add_keys)))
        if del_keys:
            printf("-> removing: %s\n"%(', '.join(del_keys)))
        
        if not len(add_keys|del_keys): 
            printf("-> No updates needed\n")
            return 

        for inline in default.split("\n"):
            if inline and inline[0] not in "# \t\n":

                key, value, comment = parse("{}={}//{}", inline)
                key = key.strip().rstrip()

                if key not in add_keys:
                    value = current_param[key]
                outline = (
                    "%-24s" % ("%-12s" % key + "= " +str(value)) +
                    " //" + comment)
            else: outline = inline

            fpo.write("%s\n" % outline)
        fpi.close()
        fpo.close()
        os.system("mv /tmp/starbug.param %s" % f_name)
    else:
        p_error("local parameter file '%s' does not exist\n" % f_name)

