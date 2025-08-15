from typing import Any

#@title **COMPREHENSIVE STATISTICAL TESTING FOR DIFFERENTIAL EXPRESSION**
import os
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, mannwhitneyu
from statsmodels.stats.multitest import multipletests

print("Comprehensive Statistical Testing...")

# Create output directory
os.makedirs('deg_analysis_output', exist_ok=True)

# Load the filtered expression data (assuming from previous steps)
try:
    expression_filtered = pd.read_csv('deg_analysis_output/expression_filtered.csv', index_col=0)
    sample_metadata = pd.read_csv('deg_analysis_output/sample_metadata.csv', index_col=0)
    print(f"Loaded filtered expression data: {expression_filtered.shape}")
except FileNotFoundError:
    print("Creating mock filtered expression data for demonstration...")
    np.random.seed(42)
    n_genes = 15000
    n_samples = 21

    base_expression = np.random.lognormal(mean=5, sigma=2, size=(n_genes, n_samples))
    gene_names = [f"GENE_{i:05d}" for i in range(n_genes)]
    sample_names = [f"GSM603087{i}" if i < 10 else f"GSM60308{i}" for i in range(n_samples)]

    expression_filtered = pd.DataFrame(base_expression, index=gene_names, columns=sample_names)

    sample_metadata = pd.DataFrame({
        'Sample_ID': sample_names,
        'Group': ['SCFA_Treated'] * 11 + ['Control'] * 10
    }).set_index('Sample_ID')

    # Add differential expression to some genes
    n_deg_genes = 500
    deg_indices = np.random.choice(n_genes, n_deg_genes, replace=False)
    for idx in deg_indices:
        if np.random.random() > 0.5:
            expression_filtered.iloc[idx, :11] *= np.random.uniform(2, 4)
        else:
            expression_filtered.iloc[idx, :11] *= np.random.uniform(0.25, 0.5)

# Define sample groups
treated_group = 'SCFA_Treated'
control_group = 'Control'

treated_samples = sample_metadata[sample_metadata['Group'] == treated_group].index.tolist()
control_samples = sample_metadata[sample_metadata['Group'] == control_group].index.tolist()

print(f"Treatment samples: {len(treated_samples)}")
print(f"Control samples: {len(control_samples)}")

# Log2 transformation for statistical analysis
expression_log2 = np.log2(expression_filtered + 1)

# Comprehensive statistical testing
print("\nPerforming comprehensive statistical tests...")
deg_results = []

for gene in expression_log2.index:
    treated_expr = expression_log2.loc[gene, treated_samples].values
    control_expr = expression_log2.loc[gene, control_samples].values

    # Skip genes with insufficient data
    if len(treated_expr) < 3 or len(control_expr) < 3:
        continue

    # Basic statistics (use sample std: ddof=1)
    treated_mean = np.mean(treated_expr)
    control_mean = np.mean(control_expr)
    treated_std = np.std(treated_expr, ddof=1)
    control_std = np.std(control_expr, ddof=1)

    log2_fold_change = treated_mean - control_mean
    # Keep signed fold-change (FC > 1 means up in treated)
    fold_change = 2 ** (log2_fold_change)

    # 1. Student's t-test (Welch)
    try:
        t_stat, t_pval = ttest_ind(treated_expr, control_expr, equal_var=False, nan_policy='omit')
    except Exception:
        t_stat, t_pval = np.nan, 1.0

    # 2. Mann-Whitney U test (non-parametric)
    try:
        u_stat, u_pval = mannwhitneyu(treated_expr, control_expr, alternative='two-sided')
    except Exception:
        u_stat, u_pval = np.nan, 1.0

    # 3. Effect size (Cohen's d) using pooled sample std
    pooled_numerator = ((len(treated_expr) - 1) * (treated_std ** 2) +
                        (len(control_expr) - 1) * (control_std ** 2))
    pooled_denominator = (len(treated_expr) + len(control_expr) - 2)
    pooled_std = np.sqrt(pooled_numerator / pooled_denominator) if pooled_denominator > 0 else 0.0
    cohens_d = (treated_mean - control_mean) / pooled_std if pooled_std > 0 else 0.0

    deg_results.append({
        'Gene': gene,
        'Treated_Mean': treated_mean,
        'Control_Mean': control_mean,
        'Treated_Std': treated_std,
        'Control_Std': control_std,
        'Log2_Fold_Change': log2_fold_change,
        'Fold_Change': fold_change,
        'T_Statistic': t_stat,
        'T_P_Value': t_pval,
        'U_Statistic': u_stat,
        'U_P_Value': u_pval,
        'Cohens_D': cohens_d
    })

# Convert to DataFrame (handle case where no genes were processed)
if len(deg_results) == 0:
    deg_df = pd.DataFrame(columns=[
        'Gene','Treated_Mean','Control_Mean','Treated_Std','Control_Std',
        'Log2_Fold_Change','Fold_Change','T_Statistic','T_P_Value',
        'U_Statistic','U_P_Value','Cohens_D'
    ])
else:
    deg_df = pd.DataFrame(deg_results)

print("Applying multiple hypothesis correction...")

# Safe multiple testing correction (if there are p-values)
if deg_df.shape[0] > 0:
    try:
        deg_df['T_P_Adjusted'] = multipletests(deg_df['T_P_Value'].fillna(1).astype(float), method='fdr_bh')[1]
        deg_df['U_P_Adjusted'] = multipletests(deg_df['U_P_Value'].fillna(1).astype(float), method='fdr_bh')[1]
    except Exception as e:
        print("Warning: multipletests failed:", str(e))
        deg_df['T_P_Adjusted'] = 1.0
        deg_df['U_P_Adjusted'] = 1.0
else:
    deg_df['T_P_Adjusted'] = []
    deg_df['U_P_Adjusted'] = []

# Define thresholds
p_threshold = 0.05
fc_threshold = 1.5
log2fc_threshold = np.log2(fc_threshold)

# Classify significance safely (works even if deg_df empty)
if deg_df.shape[0] > 0:
    deg_df['T_Significant'] = (deg_df['T_P_Adjusted'] < p_threshold) & (deg_df['Log2_Fold_Change'].abs() > log2fc_threshold)
    deg_df['U_Significant'] = (deg_df['U_P_Adjusted'] < p_threshold) & (deg_df['Log2_Fold_Change'].abs() > log2fc_threshold)
    deg_df['Significant'] = deg_df['T_Significant'] | deg_df['U_Significant']
    deg_df['Regulation'] = np.where(deg_df['Log2_Fold_Change'] > log2fc_threshold, 'Upregulated',
                                   np.where(deg_df['Log2_Fold_Change'] < -log2fc_threshold, 'Downregulated', 'Not_Changed'))
    deg_df = deg_df.sort_values('T_P_Adjusted')
else:
    # ensure columns exist even when empty
    for col in ['T_Significant','U_Significant','Significant','Regulation']:
        deg_df[col] = pd.Series(dtype=object)

print(f"Statistical analysis completed for {len(deg_df)} genes")
if 'Significant' in deg_df.columns:
    print(f" - Significant genes (consensus): {int(deg_df['Significant'].sum()) if len(deg_df)>0 else 0}")
    print(f" - Upregulated: {int((deg_df['Regulation'] == 'Upregulated').sum()) if len(deg_df)>0 else 0}")
    print(f" - Downregulated: {int((deg_df['Regulation'] == 'Downregulated').sum()) if len(deg_df)>0 else 0}")

# Save results
deg_df.to_csv('deg_analysis_output/comprehensive_deg_results.csv', index=False)
print("Saved: deg_analysis_output/comprehensive_deg_results.csv")

