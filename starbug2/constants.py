STARBUG_DATA_DIR = "STARBUG_DATDIR"

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

# test states
EXIT_SUCCESS = 0
EXIT_FAIL = 1
EXIT_EARLY = 2
EXIT_MIXED = 3

# tag used for param file.
PARAM_FILE_TAG = "PARAMFILE"
REGION_TAB = "REGION_TAB"

# mode labels.
DETECTION = "DETECTION"
BACKGROUND = "BACKGROUND"
APPHOT = "APPHOT"
PSFPHOT = "PSFPHOT"
MATCHOUTPUTS = "MATCHOUTPUTS"

# ?????
OUTPUT = "OUTPUT"
VERBOSE_TAG = "VERBOSE"
NCORES = "NCORES"

# text based logo (using raw string to bypass escape characters)
LOGO = r"""
                   *          *  __  *  __   - * --   - 
 STARBUGII               *      / ___ /    \   --  -      - 
 ---------                 *___---.    .___/  -   --   -
 JWST photometry in       ./=== \  \.     \      * 
 complex crowded fields   | (O)  |  |     |           *
                           \._._/ ./    _(\)   *   
 conor.nally@ed.ac.uk     /   ~--\ ----~   \      *
                        ---      ___       ---      
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
            point source locations, their geometric properties and flux/magnitude
            measurements as calculated by aperture photometry. The output file
            will have the suffix "-ap", note this is the same as the output
            for the aperture photometry routine.
        
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
            the flux contribution on a source brightness from the dusty environment.
            
            The routine requires a list of sources to be generated (by source 
            detection) or loaded with [-d sourcelist.fits] and requires a FITS
            image to work on. The routine will ouput a FITS image, with the same
            dimensions and spatial coverage as the input image, with the suffix
            "-bgd". This background image can be used in the photometry later.
        
            To run the routine, use the core command:
        
                $~ starbug2 -B -d sourcelist.fits image.fits
        
            Alter the parameter file options under "BACKGROUND ESTIMATION" 
            to tune the performance of starbug2. Two key parameters are:
        
                - BGD_R    : Set a fixed aperture mask radius around each source
                - BOX_SIZE : Set the estimation resolution (larger will be more blurred)
        
            Full documentation is at https://starbug2.readthedocs.io
        """,
    APPHOT:
        """
            Aperture Photometry
            -------------------
        
            This routine conducts aperture photometry on an image given a list
            of sources. It requires a FITS image to run on and a FITS table source
            list with either RA/DEC columns, or x/y_centroid or x/y_0 columns. 
            The routine outputs a table with the suffix "-ap". Note this filename
            is the same as the source detection routine because aperture photometry 
            is automatically run at the end of the source detection step. The output
            table contains 2flux/magnitude information on every source
        
            To run this routine, use the core command:
        
                $~ starbug2 -A -d sourcelist.fits image.fits
        
            Alter the parameter file options under "APERTURE PHOTOMETRY" to tune
            the performance of starbug2. Three key parameters are:
        
                - APPHOT_R : Set the aperture radius for photometry (in pixels)
                - SKY_RIN  : Set the inner sky annulus radius (in pixels)
                - SKY_ROUT : Set the outer sky annulus radius (in pixels)
        
            Full documentation is at https://starbug2.readthedocs.io
        """,
    PSFPHOT:
        """
            PSF Photometry
            --------------
        
            This routine conducts PSF fitting photometry on an image given
            a list of sources. Its requires a FITS image to run on and a FITS table
            sourcelist with either RA/DEC columns, or x/y_centroid or x/y_0 columns.
            The routine outputs a table with the suffix "-psf". The output table 
            contains 2flux/magnitude information on every source
            
            To run this routine, use the core command:
                
                $~ starbug2 -P -d sourcelist image.fits 
        
            Alter the parameter file options under "PHOTOMETRY" to tune the
            performance of starbug2. Two key parameters are:
        
                - FORCE_POS    : Hold the cetroid positions of source fixed (forced photometry)
                - GEN_RESIDUAL : Generate a residual image from all the fit source
        
            Full documentation is at https://starbug2.readthedocs.io
        """,
    MATCHOUTPUTS:
        """
            Match Outputs
            -------------
        
            This option is set if the user wishes to combine all the output catalogues
            from starbug together. It would be used in the case that a routine is
            being ran on a list of images (either in series or parallel) and the
            final catalogues should all be combined into a single source list.
            It outputs two files, one with the suffix "full" and another with "match".
            The first is all columns from all table preserved into a single large
            catalogue, the second averages all the similar columns into a reduced
            table.
        
            To run this routine, use the core code:
                
                $~ starbug2 -DM image1.fits image2.fits image3.fits ...
        
            Alter the parameter file options under "CATALOGUE MATCHING" to tune the
            performance of starbug2. Two key parameters are:
        
                - MATCH_THRESH : Set the separation threshold (arcsec) to match two sources
                - NEXP_THRESH  : Set the minimum number of catalogues a source must be present in
        """,
}