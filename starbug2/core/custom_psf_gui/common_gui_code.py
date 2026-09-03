import copy
import os
import sys
from pathlib import Path
from typing import Tuple

import numpy as np
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication
from astropy.table import Table

from starbug2.constants import ExitStates, STAR_BUG_TEST_DAT_ENV
from starbug2.core.star_bug_config import StarBugMainConfig
from starbug_main import StarbugBase

# size of grid for ech star from centre
STAR_IMAGE_SIZE: int = 25

def detect_stars(
    config: StarBugMainConfig) -> (
        Tuple[np.ndarray | None, Table | None, ExitStates]):
    """
    runs basic starbug to generate the image and detections.
    :param config: the main config which will contain detection params.
    :type config: StarBugMainConfig
    :return: a tuple containing the main image, and the detections, and
             exit state. If not successful, main image and detections will
             be None.
    :rtype: Tuple[np.ndarray | None, Table | None, ExitStates]
    """
    config_copy = update_config(config)
    return run_starbug_for_detection(config_copy)


def run_starbug_for_detection(config: StarBugMainConfig) -> (
    Tuple[np.ndarray | None, Table | None, ExitStates]):
    """
    runs basic starbug to generate the image and detections.
    :param config: the config
    :type config: StarBugMainConfig
    :return: a tuple containing the main image, and the detections, and
             exit state. If not successful, main image and detections will
             be None.
    :rtype: Tuple[np.ndarray | None, Table | None, ExitStates]
    """
    # create new base and execute
    starbug_base: StarbugBase = StarbugBase(
        config.fits_images[0], config, ap_file=None,
        bkg_file=None)
    result: ExitStates = starbug_base.run_starbug(config)

    # if not successful, pass upwards.
    if result != ExitStates.EXIT_SUCCESS:
        return None, None, result
    else:
        return (starbug_base.main_image().data, starbug_base.detections,
                result)


def run_starbug_for_image_and_ap_file(config: StarBugMainConfig) -> (
        Tuple[np.ndarray | None, Table | None, ExitStates]):
    """
    reads in the ap file and image.
    :param config: the main config file
    :type config: StarBugMainConfig
    :return: a tuple containing the main image, and the detections, and
             exit state. If not successful, main image and detections will
             be None.
    :rtype: Tuple[np.ndarray | None, Table | None, ExitStates]
    """
    starbug_base: StarbugBase = StarbugBase(
        config.fits_images[0], config, ap_file=config.ap_file,
        bkg_file=config.background_file)
    return (starbug_base.main_image().data, starbug_base.detections,
            ExitStates.EXIT_SUCCESS)


def update_config(config: StarBugMainConfig) -> StarBugMainConfig:
    """
    updates a config to just do star detection and aperture photometry.
    :param config: the main config.
    :type config: StarBugMainConfig
    :return: a modified config.
    :rtype: StarBugMainConfig
    """
    config_copy: StarBugMainConfig = copy.deepcopy(config)
    config_copy.unfreeze()
    config_copy.do_star_detection = True
    config_copy.do_aperture_photometry = True
    config_copy.do_photometry_routine = False
    config_copy.do_bgd_estimate = False
    config_copy.clean_sources = False

    # turn off any other functionality.
    config_copy.execute_jwst_initialisation = False
    config_copy.update_param = False
    config_copy.generate_local_param_file = False
    config_copy.generate_region = False
    config_copy.do_artificial_star_test = False
    config_copy.generate_psf = False
    config_copy.do_custom_psf_gui = False
    config_copy.freeze()
    return config_copy


def register_desktop_entry(
    icon_path: str, app_id: str = "starbug2") -> None:
    """Creates a temporary .desktop entry so Linux window managers map
    the taskbar icon correctly."""
    desktop_dir = Path.home() / ".local" / "share" / "applications"
    desktop_dir.mkdir(parents=True, exist_ok=True)

    desktop_file = desktop_dir / f"{app_id}.desktop"

    # Simple launcher definition linking the App ID to the image path
    content = f"""[Desktop Entry]
    Type=Application
    Name=starbug2 PSF generator
    Exec=python3
    Icon={os.path.abspath(icon_path)}
    Terminal=false
    StartupWMClass={app_id}
    """
    try:
        desktop_file.write_text(content)
    except Exception as e:
        print(f"Warning: Could not write desktop entry: {e}")

def create_gui_instance_with_icon() -> Tuple[QApplication, QIcon]:
    """
    creates the application and icon and hands back for ui choice.
    :return: the application and the icon.
    :rtype: Tuple[QApplication, QIcon]
    """
    test_path: str | None = os.environ.get(STAR_BUG_TEST_DAT_ENV)
    assert test_path is not None
    test_path: str = os.path.dirname(os.path.dirname(
        os.path.dirname(test_path)))
    image_path: str = os.path.join(
        test_path, "docs", "source", "_static", "images")
    icon_path: str = os.path.join(image_path, "starbug.png")

    # this allows debugging to work as well.
    gui_app = QApplication.instance() or QApplication(sys.argv)
    app_id = "starbug2-psf-generator"

    if os.path.exists(icon_path):
        register_desktop_entry(icon_path, app_id)

    gui_app.setApplicationName("starbug2 PSF generator")
    gui_app.setApplicationDisplayName("starbug2 PSF generator")
    gui_app.setDesktopFileName(app_id)

    # try getting the icon to appear.
    app_icon: QIcon = QIcon()
    if os.path.exists(icon_path):
        app_icon = QIcon(icon_path)
        if not app_icon.isNull():
            gui_app.setWindowIcon(app_icon)
        else:
            print(
                f"Warning: Icon file found at {icon_path} but image "
                f"payload is invalid."
            )
    else:
        print(f"Warning: Icon file not found at {icon_path}")
    return gui_app, app_icon