#!/usr/bin/env python3
"""
Combine per-source summaries from multiple model grids into one CSV.
Each source appears only once, based on the lowest chi2 across all grids.
Includes a column indicating which model grid the best fit came from.

Example:
    python sedfit_best_model_all_grids.py \
        --input_dir /path/to/modelgrid_Fit_summaries \
        --output_file /path/to/combined_summary.csv
"""


import pandas as pd
import glob
import os
import argparse

def combine_summaries(input_dir, output_file):
    """
    Combine all summary CSVs into one DataFrame, keeping best chisq per source.

    Parameters
    ----------
    input_dir : str
        Directory containing per-grid summary CSV files.
    output_file : str
        Path to save the combined summary CSV.
    """
    # Find all summary files
    all_files = glob.glob(os.path.join(input_dir, "ysofit_per_source_summary_*.csv"))
    if not all_files:
        raise FileNotFoundError("No summary files found in input directory.")

    dfs = []
    for f in all_files:
        print(f"Reading: {f}")
        df = pd.read_csv(f)
        dfs.append(df)

    # Combine all columns (outer join)
    combined_df = pd.concat(dfs, axis=0, ignore_index=True, sort=False)

    # Sort by chi² and keep best per source
    combined_df = combined_df.sort_values(by="chi2", ascending=True)
    best_df = combined_df.groupby("source_name", as_index=False).first()

    # Ensure model_grid column is present
    if "model_grid" not in best_df.columns:
        best_df["model_grid"] = ""

    # Sort by model_grid before saving
    best_df = best_df.sort_values(by="model_grid")

    # Save combined file
    best_df.to_csv(output_file, index=False)
    print(f"Combined summary saved to {output_file}")
    print(f"Sources combined: {len(best_df)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Combine model grid summaries into one CSV (best chi2 per source).")
    parser.add_argument("--input_dir", required=True, help="Directory containing summary CSV files")
    parser.add_argument("--output_file", help="Path to save combined summary CSV (default: input_dir/combined_bestModel_fit.csv)")
    args = parser.parse_args()

    # Set default output file if not provided
    output_file = args.output_file if args.output_file else os.path.join(args.input_dir, "combined_bestModel_fit.csv")

    combine_summaries(args.input_dir, output_file)
