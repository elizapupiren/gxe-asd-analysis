#!/usr/bin/env python3
"""
Pre-processing script for GxE ASD analysis.
Performs normalization and cleaning of raw expression data.
"""

import pandas as pd
import numpy as np
import argparse
import sys

def load_expression_data(filepath):
    """Load raw expression matrix from CSV file.
    
    Args:
        filepath (str): Path to CSV file with genes as rows, samples as columns.
    
    Returns:
        pandas.DataFrame: Expression matrix.
    """
    # TODO: implement loading
    print(f"Loading data from {filepath}")
    df = pd.read_csv(filepath, index_col=0)
    return df


def handle_missing_values(df, method='drop'):
    """Handle missing values in expression matrix.
    
    Args:
        df (pandas.DataFrame): Expression matrix.
        method (str): Either 'drop' to remove genes with missing values,
                     or 'mean' to impute with gene-wise mean.
    
    Returns:
        pandas.DataFrame: Processed matrix.
    """
    # TODO: implement missing value handling
    if method == 'drop':
        df = df.dropna(axis=0)
    elif method == 'mean':
        df = df.fillna(df.mean(axis=0))
    else:
        raise ValueError(f"Unknown method: {method}")
    return df


def log2_normalize(df):
    """Apply log2 transformation to expression values.
    
    Assumes data is in linear scale (e.g., counts). Adds a pseudocount of 1
    to avoid log(0).
    
    Args:
        df (pandas.DataFrame): Expression matrix.
    
    Returns:
        pandas.DataFrame: Log2-normalized matrix.
    """
    # TODO: ensure values are non-negative
    df_log = np.log2(df + 1)
    return df_log


def main():
    parser = argparse.ArgumentParser(
        description='Normalize and clean raw expression data.'
    )
    parser.add_argument('input', help='Path to raw expression CSV')
    parser.add_argument('output', help='Path for processed output CSV')
    parser.add_argument('--missing', default='drop',
                        help="Missing value handling: 'drop' or 'mean'")
    args = parser.parse_args()
    
    # Load raw data
    df = load_expression_data(args.input)
    
    # Handle missing values
    df = handle_missing_values(df, method=args.missing)
    
    # Log2 normalization
    df = log2_normalize(df)
    
    # Save processed data
    df.to_csv(args.output)
    print(f"Processed data saved to {args.output}")


if __name__ == '__main__':
    main()