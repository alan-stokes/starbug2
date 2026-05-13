import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter("ignore", category=AstropyWarning)
warnings.simplefilter("ignore", category=RuntimeWarning) ## bit dodge that

motd = "https://starbug2.readthedocs.io/en/latest/"

from os import getenv

_ = getenv("STARBUG_DATDIR")
DATDIR = _ if _ else "%s/.local/share/starbug"%(getenv("HOME"))
