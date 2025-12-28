# GxE ASD Analysis

A reproducible bioinformatics workflow to identify genes whose expression is affected by an environmental factor differently in an autism spectrum disorder (ASD) genetic model versus a control.

## Objective

This project aims to perform a gene-environment interaction analysis using hypothetical gene expression data from brain tissue of two groups of mice: wild-type (control) and a genetic model for ASD (e.g., *Cntnap2* knockout). Within each group, some mice were exposed to a specific environmental pollutant, while others were not. The goal is to find genes showing a significant "genotype-by-environment" interaction.

## Project Structure

- `data/`: Contains raw and processed data files (hypothetical)
- `scripts/`: Analysis scripts for pre-processing, statistical analysis, and visualization
- `results/`: Output files, plots, and summary tables

## Workflow

1. Data pre-processing: Normalization and cleaning of raw expression data
2. Statistical analysis: Two-way ANOVA to test for genotype-by-environment interaction effects
3. Visualization: Volcano plot highlighting significant interaction genes
