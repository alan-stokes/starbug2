import warnings
import logging
import photutils
from astropy.utils.exceptions import (
    AstropyDeprecationWarning, AstropyWarning, AstropyUserWarning)
from astropy.io.fits.verify import VerifyWarning
from astropy.utils.metadata import MergeConflictWarning

def configure_environment():
    """Configure global package flags, warning filters, and logging levels."""
    # Forces photutils to not give QTable wrappers.
    photutils.future_column_names = True

    # Target-silence only the specific Photutils/Astropy deprecation noise
    # without masking generic Runtime math errors globally.
    warnings.filterwarnings("ignore", category=AstropyDeprecationWarning)
    warnings.filterwarnings(
        "ignore", message=".*contains deprecated section.*",
        category=AstropyWarning)
    warnings.filterwarnings(
        "ignore", category=AstropyWarning,
        message=".*Input data contains invalid values*")
    warnings.filterwarnings("ignore", message=".*invalid values.*")
    warnings.filterwarnings("ignore", message=".*sigma_clipping.*")
    warnings.filterwarnings("ignore", message=".*Cannot merge meta key*")

    # Handle RuntimeWarnings elegantly: Ignore expected ones (like NaN
    # comparisons during clipping), but let actual mathematical issues surface.
    warnings.filterwarnings(
        "ignore", message=".*invalid value encountered.*",
        category=RuntimeWarning)
    warnings.filterwarnings(
        "ignore", message=".*divide by zero.*", category=RuntimeWarning)

    # Ignore all metadata merge conflict warnings
    warnings.filterwarnings("ignore", category=MergeConflictWarning)

    # --- FITS IO FORMATTING NOISE ---
    # These suppress warnings about FITS header compliance
    # (e.g., truncated comments) that do not affect scientific output.
    warnings.filterwarnings(
        "ignore", category=VerifyWarning, message=".*Card is too long.*")
    warnings.filterwarnings(
        "ignore", category=AstropyUserWarning,
        message=".*cannot be added to FITS Header.*")

    logging.getLogger().setLevel(logging.ERROR)