#!/usr/bin/env python3
"""
Visualization script for GxE ASD analysis.
Generates a volcano plot from interaction analysis results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

def load_analysis_results(filepath):
    """Load results from two-way ANOVA analysis.
    
    Args:
        filepath (str): Path to CSV with columns: gene, p_value.
    
    Returns:
        pandas.DataFrame: Results DataFrame.
    """
    df = pd.read_csv(filepath)
    return df


def compute_log_p_values(p_values):
    """Compute -log10(p-value) for volcano plot.
    
    Args:
        p_values (pd.Series): Raw p-values.
    
    Returns:
        pd.Series: -log10(p-value).
    """
    # Avoid log(0) by replacing zero with smallest positive float
    p_values_clipped = p_values.clip(lower=np.finfo(float).eps)
    return -np.log10(p_values_clipped)


def create_volcano_plot(results_df, output_path, p_threshold=0.05, logfc_threshold=None):
    """Generate volcano plot highlighting significant interaction genes.
    
    Args:
        results_df (pd.DataFrame): Must contain 'gene', 'p_value' columns.
        output_path (str): Path to save plot image.
        p_threshold (float): Significance threshold for p-value.
        logfc_threshold (float): Optional fold-change threshold (if fold-change column exists).
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Compute -log10(p)
    results_df['neg_log10_p'] = compute_log_p_values(results_df['p_value'])
    
    # Determine significance
    # If fold-change column exists, use it; otherwise just use p-value
    if 'log2_fold_change' in results_df.columns:
        x = results_df['log2_fold_change']
        xlabel = 'log2(Fold Change)'
        if logfc_threshold is not None:
            sig = (results_df['p_value'] < p_threshold) & (np.abs(x) > logfc_threshold)
        else:
            sig = results_df['p_value'] < p_threshold
    else:
        # No fold-change data, just plot -log10(p) vs index (or random jitter)
        x = np.random.randn(len(results_df)) * 0.1  # dummy x for visualization
        xlabel = 'Random jitter (no fold-change data)'
        sig = results_df['p_value'] < p_threshold
    
    # Scatter plot
    ax.scatter(x[~sig], results_df.loc[~sig, 'neg_log10_p'],
               alpha=0.5, c='gray', label='Not significant', s=20)
    ax.scatter(x[sig], results_df.loc[sig, 'neg_log10_p'],
               alpha=0.7, c='red', label='Significant', s=30)
    
    # Add significance thresholds
    ax.axhline(y=-np.log10(p_threshold), color='black', linestyle='--', linewidth=1,
               label=f'p = {p_threshold}')
    if logfc_threshold is not None and 'log2_fold_change' in results_df.columns:
        ax.axvline(x=logfc_threshold, color='blue', linestyle=':', linewidth=1)
        ax.axvline(x=-logfc_threshold, color='blue', linestyle=':', linewidth=1)
    
    # Labels and title
    ax.set_xlabel(xlabel)
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('Volcano Plot: Genotype-by-Environment Interaction')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Volcano plot saved to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate volcano plot from interaction analysis results.'
    )
    parser.add_argument('results', help='Path to analysis results CSV')
    parser.add_argument('output', help='Path for output plot image (PNG)')
    parser.add_argument('--p_threshold', type=float, default=0.05,
                        help='Significance threshold for p-value (default: 0.05)')
    parser.add_argument('--logfc_threshold', type=float, default=None,
                        help='Fold-change threshold for significance (optional)')
    args = parser.parse_args()
    
    # Load results
    results = load_analysis_results(args.results)
    
    # Create volcano plot
    create_volcano_plot(results, args.output,
                        p_threshold=args.p_threshold,
                        logfc_threshold=args.logfc_threshold)


if __name__ == '__main__':
    main()