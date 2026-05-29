#!/usr/bin/env python
# coding: utf-8

"""
SED fitting script using sedfitter (Richardson+ 2024 YSO models)
Handles zipped model directories: extract - fit - recompress - cleanup
"""

import os
import shutil
import tarfile
import logging
from datetime import datetime

import numpy as np
from astropy import units as u
from sedfitter import fit, write_parameters, write_parameter_ranges
from sedfitter.extinction import Extinction

# -----------------------------
# User Configuration
# -----------------------------

#phot_file_to_fit = "/Users/olivia.jones/Science/N79/results/YSOmodels/N79_YSO_Candidate_Catalog.csv"

phot_file_to_fit = "/Users/olivia.jones/Science/N79/results/YSOmodels/N79_YSO_Candidate_Catalog_shortList.csv"
extinction_file = "/Users/olivia.jones/Python/MyPython/YSOfit/kmh94.par"

#r24_modeldir = "/Users/olivia.jones/Science/N79/results/YSOmodels/"
#output_dir = "/Users/olivia.jones/Science/N79/results/YSOmodels/"

r24_modeldir = "/Volumes/T7/YSOmodels_20251117/"
output_dir = "/Volumes/T7/N79/results/YSOmodels/N6/"

os.makedirs(output_dir, exist_ok=True)

observed_filters = [
    "F115W", "F187N", "F200W", "F277W", "F335M", "F444W",
    "F770W", "F1000W", "F1500W", "F2100W"
]

distance_range = [45, 55] * u.kpc

#av_range=[0.0, 10.0]      # Stars usually have lower extinction
av_range=[0.0, 20.0]      # For a YSO first run



# -----------------------------
# YSO Model Geometry Options
# -----------------------------
geometry_list = [
    's---s-i', 'sp--s-i', 'sp--h-i',
    's---smi', 'sp--smi', 'sp--hmi',
    's-p-smi', 's-p-hmi', 's-pbsmi', 's-pbhmi',
    's-u-smi', 's-u-hmi', 's-ubsmi', 's-ubhmi',
    'spu-smi', 'spu-hmi', 'spubsmi', 'spubhmi'
]



#geometry_list = ['s---s-i']  # This is a simple star model. Add more geometries if needed
#geometry_list = ['spubhmi']  # This is the most populous geometry


# -----------------------------
# Validate Inputs
# -----------------------------

if not os.path.exists(phot_file_to_fit):
    raise FileNotFoundError(f"Photometry file not found: {phot_file_to_fit}")

if not os.path.exists(extinction_file):
    raise FileNotFoundError(f"Extinction file not found: {extinction_file}")

# -----------------------------
# Configure Logging
# -----------------------------

run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
log_file = os.path.join(output_dir, f"sedfitter_run_{run_id}.log")

logging.basicConfig(
    filename=log_file,
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
logging.getLogger().addHandler(console_handler)

logging.info("Starting SED fitting process")

# -----------------------------
# JWST Filter Apertures (arcsec)
# -----------------------------
filter_apertures = {
    "F070W": ["NIRCam", 0.029], "F090W": ["NIRCam", 0.033], "F115W": ["NIRCam", 0.040],
    "F140M": ["NIRCam", 0.048], "F150W": ["NIRCam", 0.050], "F150W2": ["NIRCam", 0.045],
    "F182M": ["NIRCam", 0.062], "F187N": ["NIRCam", 0.064], "F200W": ["NIRCam", 0.066],
    "F210M": ["NIRCam", 0.071], "F212N": ["NIRCam", 0.072], "F250M": ["NIRCam", 0.085],
    "F277W": ["NIRCam", 0.092], "F300M": ["NIRCam", 0.100], "F322W2": ["NIRCam", 0.096],
    "F335M": ["NIRCam", 0.111], "F356W": ["NIRCam", 0.116], "F360M": ["NIRCam", 0.120],
    "F410M": ["NIRCam", 0.137], "F430M": ["NIRCam", 0.144], "F444W": ["NIRCam", 0.145],
    "F460M": ["NIRCam", 0.157], "F480M": ["NIRCam", 0.164], "F560W": ["MIRI", 0.207],
    "F770W": ["MIRI", 0.269], "F1000W": ["MIRI", 0.328], "F1130W": ["MIRI", 0.375],
    "F1280W": ["MIRI", 0.420], "F1500W": ["MIRI", 0.488], "F1800W": ["MIRI", 0.591],
    "F2100W": ["MIRI", 0.674], "F2550W": ["MIRI", 0.803]
}

# -----------------------------
# Build filter and aperture lists dynamically
# -----------------------------
filters = []
apertures = []

for f in observed_filters:
    if f in filter_apertures:
        instrument, aperture_size = filter_apertures[f]
        filters.append(f"JWST/{instrument}.{f}")
        apertures.append(aperture_size)
    else:
        logging.warning(f"Filter {f} not found in JWST aperture dictionary. Skipping.")

if not filters:
    raise ValueError("No valid JWST filters found. Check observed_filters list.")

apertures = np.array(apertures) * u.arcsec

logging.info(f"Selected Filters: {filters}")
logging.info(f"Apertures (arcsec): {apertures}")

# -----------------------------
# Load Extinction Law
# -----------------------------

extinction = Extinction.from_file(
    extinction_file,
    columns=[0, 3],
    wav_unit=u.micron,
    chi_unit=u.cm**2 / u.g
)

# -----------------------------
# Downselection format
# -----------------------------
# By default, functions in `sedfitter` will work with the best-fitting model (i.e. `('N', 1)`).
# Keep only models within  $\chi^2-\chi^2_{\rm best}$-per-data-point $<$ 3.
select_format = ('F', 3)

# # -----------------------------
# # Helper: Compress model directory
# # -----------------------------
#
# def compress_model(model_path, tar_path):
#     """Compress a model directory into a tar.gz archive."""
#     try:
#         with tarfile.open(tar_path, "w:gz") as tar_ref:
#             tar_ref.add(model_path, arcname=os.path.basename(model_path))
#         logging.info(f"Compressed {model_path} to {tar_path}")
#     except Exception as e:
#         logging.error(f"Compression failed for {model_path}: {e}")


# -----------------------------
# Loop over geometries
# -----------------------------

for geometry in geometry_list:
    logging.info(f"Processing geometry: {geometry}")

    tar_path = os.path.join(r24_modeldir, f"{geometry}.tar.gz")
    extract_root = os.path.join(r24_modeldir, "temp_extract")
    model_path = os.path.join(extract_root, "r+24_models-1.2", geometry)

    # Extract tar.gz
    if os.path.exists(tar_path):
        logging.info(f"Extracting {tar_path} to {extract_root}")
        os.makedirs(extract_root, exist_ok=True)

        with tarfile.open(tar_path, "r:gz") as tar_ref:
            tar_ref.extractall(extract_root, filter="data")
    else:
        logging.error(f"Tar file not found: {tar_path}. Skipping.")
        continue

    # Validate extracted model
    if not os.path.exists(os.path.join(model_path, "models.conf")):
        logging.error(f"Missing models.conf in {model_path}. Skipping.")
        shutil.rmtree(extract_root, ignore_errors=True)
        continue

    # Run SED fitting
    output_file = os.path.join(output_dir, f"fitinfo_{geometry}_{run_id}.fits")

    try:
        start_time = datetime.now()

        fit(
            phot_file_to_fit,
            filters,
            apertures,
            model_path,
            output_file,
            extinction_law=extinction,
            distance_range=distance_range,
            av_range=[0.0, 20.0]     #av_range
        )

        # # R24 models include parameters that vary with aperture (e.g., sphere masses). Apertures are log-spaced.
        # # To pick the right aperture:
        # all_apertures = np.logspace(2,6,20) * u.AU
        #
        # # This finds the aperture closest to 7200 AU (corresponding to a 0.145 arcseconds aperture at 50 kpc).
        # mid_distance_ap = 7200 * u.AU
        # ap_num = np.argmin(abs(mid_distance_ap - all_apertures))



        # Write downselected parameters
        logging.info("Writing parameters...")
        write_parameters(
            output_file,
            os.path.join(output_dir, f"params_{geometry}_{run_id}.txt"),
            select_format=select_format
        )

        logging.info("Writing parameter ranges...")
        write_parameter_ranges(
            output_file,
            os.path.join(output_dir, f"param_ranges_{geometry}_{run_id}.txt"),
            select_format=select_format
        )

        # #  Plot SEDs for selected models
        # logging.info("Generating plots...")
        # plot(
        #     output_file,
        #     os.path.join(output_dir, f"plots_{geometry}_{run_id}"),
        #     plot_max=10,
        #     select_format=select_format
        # )

        elapsed = (datetime.now() - start_time).total_seconds()
        logging.info(f"Completed geometry: {geometry} in {elapsed:.2f}s")

        # Re-compress and cleanup
        #compress_model(model_path, tar_path)
        shutil.rmtree(extract_root, ignore_errors=True)

    except Exception as e:
        logging.error(f"Error processing geometry {geometry}: {e}")
        shutil.rmtree(extract_root, ignore_errors=True)

logging.info("SED fitting process completed.")
print(f"SED fitting completed for all geometries. Log: {log_file}")
