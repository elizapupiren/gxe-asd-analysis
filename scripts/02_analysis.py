#!/usr/bin/env python3
"""
Statistical analysis script for GxE ASD analysis.
Performs two-way ANOVA to test for genotype-by-environment interaction effects.
"""

import pandas as pd
import numpy as np
import argparse
import statsmodels.api as sm
from statsmodels.formula.api import ols
import warnings
warnings.filterwarnings('ignore')

def load_processed_data(filepath):
    """Load pre-processed expression matrix.
    
    Args:
        filepath (str): Path to CSV file with genes as rows, samples as columns.
    
    Returns:
        pandas.DataFrame: Expression matrix.
    """
    df = pd.read_csv(filepath, index_col=0)
    return df


def load_sample_metadata(filepath):
    """Load sample metadata including genotype and environment.
    
    Args:
        filepath (str): Path to CSV with columns: sample, genotype, environment.
    
    Returns:
        pandas.DataFrame: Metadata DataFrame.
    """
    meta = pd.read_csv(filepath)
    return meta


def perform_two_way_anova(expression_vector, meta):
    """Perform two-way ANOVA for a single gene.
    
    Args:
        expression_vector (pandas.Series): Expression values across samples.
        meta (pandas.DataFrame): Metadata with 'genotype' and 'environment' columns.
    
    Returns:
        float: p-value for genotype:environment interaction term.
    """
    # Combine expression with metadata
    df = pd.DataFrame({
        'expression': expression_vector.values,
        'genotype': meta['genotype'].values,
        'environment': meta['environment'].values
    })
    
    # Fit two-way ANOVA model with interaction
    model = ols('expression ~ genotype + environment + genotype:environment', data=df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    
    # Extract p-value for interaction term
    if 'genotype:environment' in anova_table.index:
        p_value = anova_table.loc['genotype:environment', 'PR(>F)']
    else:
        p_value = np.nan
    
    return p_value


def analyze_all_genes(expression_df, meta_df):
    """Run two-way ANOVA for each gene.
    
    Args:
        expression_df (pandas.DataFrame): Expression matrix (genes x samples).
        meta_df (pandas.DataFrame): Metadata for each sample.
    
    Returns:
        pandas.DataFrame: Results with gene IDs and p-values.
    """
    results = []
    for gene in expression_df.index:
        p_val = perform_two_way_anova(expression_df.loc[gene], meta_df)
        results.append({'gene': gene, 'p_value': p_val})
    
    results_df = pd.DataFrame(results)
    return results_df


def main():
    parser = argparse.ArgumentParser(
        description='Perform two-way ANOVA for genotype-by-environment interaction.'
    )
    parser.add_argument('expression', help='Path to processed expression CSV')
    parser.add_argument('metadata', help='Path to sample metadata CSV')
    parser.add_argument('output', help='Path for output results CSV')
    args = parser.parse_args()
    
    # Load data
    expr = load_processed_data(args.expression)
    meta = load_sample_metadata(args.metadata)
    
    # Ensure sample order matches
    # (Assuming column names of expr match meta['sample'])
    # TODO: implement sample alignment
    
    # Perform analysis
    results = analyze_all_genes(expr, meta)
    
    # Save results
    results.to_csv(args.output, index=False)
    print(f"Analysis results saved to {args.output}")


if __name__ == '__main__':
    main()