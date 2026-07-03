#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import argparse
import sys

def load_data(filepath, test_mode=False):
    """Load the phenotype data from a tab-separated file."""
    try:
        if test_mode:
            # Load only first 1000 rows in test mode
            df = pd.read_csv(filepath, sep='\t', nrows=10000)
            print(f"Test mode: Loaded first 10000 rows from {filepath}")
        else:
            df = pd.read_csv(filepath, sep='\t')
        return df
    except Exception as e:
        print(f"Error loading {filepath}: {e}")
        sys.exit(1)

def apply_boolean_filter(df, filter_col, file_name):
    """Apply boolean filter to dataframe, keeping only True values."""
    if filter_col not in df.columns:
        print(f"Filter column '{filter_col}' not found in {file_name}")
        print("Available columns:", list(df.columns))
        sys.exit(1)
    
    # Count original rows
    original_count = len(df)
    
    # Apply filter - keep rows where filter_col is True
    # Handle various boolean representations
    df_filtered = df[df[filter_col].astype(str).str.lower().isin(['true', '1', 'yes', 't'])]
    
    filtered_count = len(df_filtered)
    print(f"Applied boolean filter on '{filter_col}' in {file_name}: {original_count} -> {filtered_count} rows")
    
    return df_filtered

def find_shared_phenotypes(df1, df2, pheno_col1=None, pheno_col2=None):
    """Find phenotypes that exist in both datasets."""
    # Use first column if no phenotype column specified
    col1 = pheno_col1 if pheno_col1 else df1.columns[0]
    col2 = pheno_col2 if pheno_col2 else df2.columns[0]
    
    shared = set(df1[col1]) & set(df2[col2])
    return sorted(list(shared)), col1, col2

def create_scatter_plot(df1, df2, column, output_file=None, pheno_col1=None, pheno_col2=None):
    """Create a scatter plot comparing the chosen column between two datasets."""
    # Find shared phenotypes
    shared_phenos, pheno_col1, pheno_col2 = find_shared_phenotypes(df1, df2, pheno_col1, pheno_col2)
    
    if not shared_phenos:
        print("No shared phenotypes found between the two files.")
        return
    
    print(f"Found {len(shared_phenos)} shared phenotypes")
    print(f"Using phenotype columns: '{pheno_col1}' (file 1) and '{pheno_col2}' (file 2)")
    
    # Create merged dataframe with just the columns we need
    df1_subset = df1[[pheno_col1, column]].rename(columns={column: f"{column}_1"})
    df2_subset = df2[[pheno_col2, column]].rename(columns={column: f"{column}_2"})
    
    # Merge on phenotype columns
    merged = df1_subset.merge(
        df2_subset, 
        left_on=pheno_col1, 
        right_on=pheno_col2, 
        how='inner'
    )
    
    # Drop rows with any NA values in the target columns
    merged_clean = merged.dropna(subset=[f"{column}_1", f"{column}_2"])
    
    if merged_clean.empty:
        print(f"No valid (non-NA) values found for column '{column}' in shared phenotypes.")
        return
    
    print(f"Plotting {len(merged_clean)} valid data points")
    
    # Extract values directly from merged dataframe
    values1 = merged_clean[f"{column}_1"].astype(float)
    values2 = merged_clean[f"{column}_2"].astype(float)

    # Create the scatter plot
    plt.figure(figsize=(10, 8))
    plt.scatter(values1, values2, alpha=0.6, s=50)
    # Add diagonal line (y=x) for reference
    min_val = min(min(values1), min(values2))
    max_val = max(max(values1), max(values2))
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5, label='y=x')
    
    # Calculate correlation
    correlation = np.corrcoef(values1, values2)[0, 1]
    
    # Labels and title
    plt.xlabel(f'{column} (File 1)')
    plt.ylabel(f'{column} (File 2)')
    plt.title(f'Scatter Plot: {column} Comparison\n(r = {correlation:.3f}, n = {len(values1)})')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    # Make plot square
    plt.axis('equal')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()
    
    # Print some statistics
    print(f"\nStatistics for {column}:")
    print(f"Correlation coefficient: {correlation:.4f}")
    print(f"File 1 - Mean: {np.mean(values1):.4f}, Std: {np.std(values1):.4f}")
    print(f"File 2 - Mean: {np.mean(values2):.4f}, Std: {np.std(values2):.4f}")
    
    return values1, values2

def main():
    parser = argparse.ArgumentParser(description='Create scatter plot comparing a column between two phenotype files')
    parser.add_argument('file1', help='First phenotype file')
    parser.add_argument('file2', help='Second phenotype file')
    parser.add_argument('column', help='Column to compare (e.g., H2, SE, INT, etc.)')
    parser.add_argument('-o', '--output', help='Output file for the plot (optional)')
    parser.add_argument('--list-columns', action='store_true', help='List available columns and exit')
    parser.add_argument('--pheno-col1', help='Phenotype column name for file 1 (default: first column)')
    parser.add_argument('--pheno-col2', help='Phenotype column name for file 2 (default: first column)')
    parser.add_argument('--filter-col1', help='Boolean filter column for file 1 (keep only True values)')
    parser.add_argument('--filter-col2', help='Boolean filter column for file 2 (keep only True values)')
    parser.add_argument('--test', action='store_true', help='Test mode: use only first 1000 rows from each file')
    
    args = parser.parse_args()
    
    # Load the data
    df1 = load_data(args.file1, test_mode=args.test)
    df2 = load_data(args.file2, test_mode=args.test)
    
    # List columns if requested
    if args.list_columns:
        print("Available columns in file 1:", list(df1.columns))
        print("Available columns in file 2:", list(df2.columns))
        return
    
    # Apply boolean filters if specified
    if args.filter_col1:
        df1 = apply_boolean_filter(df1, args.filter_col1, args.file1)
    
    if args.filter_col2:
        df2 = apply_boolean_filter(df2, args.filter_col2, args.file2)
    
    # Check phenotype columns if specified
    if args.pheno_col1 and args.pheno_col1 not in df1.columns:
        print(f"Phenotype column '{args.pheno_col1}' not found in {args.file1}")
        print("Available columns:", list(df1.columns))
        sys.exit(1)
    
    if args.pheno_col2 and args.pheno_col2 not in df2.columns:
        print(f"Phenotype column '{args.pheno_col2}' not found in {args.file2}")
        print("Available columns:", list(df2.columns))
        sys.exit(1)
    
    # Check if the comparison column exists in both files
    if args.column not in df1.columns:
        print(f"Column '{args.column}' not found in {args.file1}")
        print("Available columns:", list(df1.columns))
        sys.exit(1)
    
    if args.column not in df2.columns:
        print(f"Column '{args.column}' not found in {args.file2}")
        print("Available columns:", list(df2.columns))
        sys.exit(1)
    
    # Create the scatter plot
    create_scatter_plot(df1, df2, args.column, args.output, args.pheno_col1, args.pheno_col2)

if __name__ == "__main__":
    main()
