#!/usr/bin/env python3
"""
Process YSO model parameter files from sedfitter, compute per-source summaries,
and save formatted CSV outputs for each model grid.

Example:
    python sedfit_grid_summary.py \
        --input_dir /path/to/input \
        [--output_dir /path/to/output]
If --output_dir is not provided, results are saved in input_dir/modelgrid_summaries.
"""

import pandas as pd
import numpy as np
import os
import glob
import argparse

# -----------------------------
def parse_sed_file_to_dataframe(filename):
    """
    Parse a sedfitter parameter file into a pandas DataFrame.

    Expected file format:
    - First three lines: metadata (ignored except for column headers on line 2)
    - Subsequent blocks: source header followed by model parameter lines

    Each source block contains:
    - A header line: <source_name> <n_data> <n_mod>
    - n_mod lines of model parameters corresponding to column headers

    Parameters
    ----------
    filename : str
        Path to the sedfitter params file.

    Returns
    -------
    df : pandas.DataFrame
        DataFrame where each row represents a model for a source.
        Includes columns: source_name, n_data, n_mod, and all model parameters.
    column_headers : list of str
        List of parameter names extracted from the file header.
    """
    rows = []
    with open(filename, 'r') as file:
        # Read metadata
        _ = file.readline().strip()  # meta1
        raw_header = file.readline().strip()  # raw header line
        _ = file.readline().strip()  # meta3

        # Known multi-word phrases to normalize
        phrases = [
            "source luminosity",
            "line-of-sight masses",
            "sphere masses",
            "inner radius",
            "outer radius",
            "spectral index",
            "disk minimum q",
            "line-of-sight mass-weighted temperatures",
            "line-of-sight photon-weighted temperatures",
            "sphere mass-weighted temperatures"
        ]

        # Replace spaces in known phrases with underscores
        for phrase in phrases:
            raw_header = raw_header.replace(phrase, phrase.replace(" ", "_"))

        # Split into tokens
        column_headers = raw_header.split()

        while True:
            header = file.readline()
            if not header:
                break
            parts = header.split()
            if len(parts) < 3:
                continue

            source_name, n_data, n_mod = parts[0], int(parts[1]), int(parts[2])

            # Read model lines for this source
            for _ in range(n_mod):
                model_line = file.readline()
                if not model_line:
                    break
                model_values = model_line.split()

                # Adjust header length to match values
                if len(model_values) < len(column_headers):
                    column_headers = column_headers[:len(model_values)]
                elif len(model_values) > len(column_headers):
                    extra_cols = [f"extra_{i}" for i in range(len(model_values) - len(column_headers))]
                    column_headers = column_headers + extra_cols

                row = {'source_name': source_name, 'n_data': n_data, 'n_mod': n_mod}
                for col, val in zip(column_headers, model_values):
                    if col == "model_name":
                        row[col] = val  # Always keep as string
                    else:
                        try:
                            row[col] = float(val)
                        except ValueError:
                            row[col] = np.nan
                rows.append(row)

    return pd.DataFrame(rows), column_headers

# -----------------------------
def compute_per_source_summary(df, model_grid):
    """
    Compute per-source summary statistics:
    - Median and MAD for all numeric columns
    - Best-fit model based on minimum chi²
    """
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    exclude_cols = ['n_data', 'n_mod']  # keep chi2 for best-fit selection
    numeric_cols = [col for col in numeric_cols if col not in exclude_cols]

    best_models = df.loc[df.groupby('source_name')['chi2'].idxmin()].copy()

    agg_funcs = {}
    for col in numeric_cols:
        agg_funcs[col + '_median'] = (col, 'median')
        agg_funcs[col + '_mad'] = (col, lambda x: np.median(np.abs(x - np.median(x))))

    agg_df = df.groupby('source_name').agg(**agg_funcs).reset_index()

    final_df = agg_df.merge(best_models[['source_name', 'n_data', 'n_mod', 'chi2', 'model_name']], on='source_name')
    final_df['model_grid'] = model_grid

    cols_order = ['source_name', 'model_grid', 'n_data', 'n_mod', 'chi2', 'model_name'] + \
                 [c for c in final_df.columns if c not in ['source_name', 'model_grid', 'n_data', 'n_mod', 'chi2', 'model_name']]
    return final_df[cols_order]

# -----------------------------
def main(input_dir, output_dir):
    """Process all params*.txt files in input_dir and save summaries to output_dir."""
    os.makedirs(output_dir, exist_ok=True)

    for filename in glob.glob(os.path.join(input_dir, "params*.txt")):
        print(f"Processing: {filename}")
        model_grid = os.path.basename(filename).split('_')[1]

        df, _ = parse_sed_file_to_dataframe(filename)
        summary_df = compute_per_source_summary(df, model_grid)

        # Format numeric columns for clean output
        for col in summary_df.select_dtypes(include=[np.number]).columns:
            #summary_df[col] = summary_df[col].apply(lambda x: f"{x:.6g}" if pd.notnull(x) else "")
            summary_df[col] = summary_df[col].apply(lambda x: f"{x:.6g}" if pd.notnull(x) else "NaN")

        out_csv = os.path.join(output_dir, f"ysofit_per_source_summary_{model_grid}.csv")
        summary_df.to_csv(out_csv, index=False)
        print(f"Saved summary to {out_csv}")

# -----------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute per-source summaries from sedfitter params files.")
    parser.add_argument("--input_dir", required=True, help="Directory containing params*.txt files")
    parser.add_argument("--output_dir", help="Directory to save summary CSV files (default: input_dir/modelgrid_summaries)")
    args = parser.parse_args()

    output_dir = args.output_dir if args.output_dir else os.path.join(args.input_dir, "modelgrid_summaries")
    main(args.input_dir, output_dir)