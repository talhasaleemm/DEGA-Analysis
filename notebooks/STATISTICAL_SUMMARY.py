#@title **COMPREHENSIVE STATISTICAL SUMMARY**


print("\n Comprehensive Statistical Summary...")

# Create comprehensive summary statistics
summary_stats = {
    'Total_Genes_Analyzed': len(deg_df),
    'Significant_Genes': deg_df['Significant'].sum(),
    'Upregulated_Genes': (deg_df['Regulation'] == 'Upregulated').sum(),
    'Downregulated_Genes': (deg_df['Regulation'] == 'Downregulated').sum(),
    'Mean_Fold_Change_Significant': deg_df[deg_df['Significant']]['Fold_Change'].mean(),
    'Max_Fold_Change': deg_df['Fold_Change'].max(),
    'Min_P_Value': deg_df['T_P_Adjusted'][deg_df['T_P_Adjusted'] > 0].min(),
    'Percent_Significant': (deg_df['Significant'].sum() / len(deg_df)) * 100
}

print("Statistical Summary:")
for key, value in summary_stats.items():
    if 'P_Value' in key:
        print(f"  {key.replace('_', ' ')}: {value:.2e}")
    elif 'Percent' in key:
        print(f"  {key.replace('_', ' ')}: {value:.2f}%")
    else:
        print(f"  {key.replace('_', ' ')}: {value:.2f}" if isinstance(value, float) else f"  {key.replace('_', ' ')}: {value}")

# Effect size analysis
print(f"\nEffect Size Analysis:")
significant_genes = deg_df[deg_df['Significant']]
if len(significant_genes) > 0:
    print(f"  Average Cohen's D for significant genes: {significant_genes['Cohens_D'].abs().mean():.3f}")
    print(f"  Large effect size (|d| > 0.8): {(significant_genes['Cohens_D'].abs() > 0.8).sum()} genes")
    print(f"  Medium effect size (0.5 < |d| ≤ 0.8): {((significant_genes['Cohens_D'].abs() > 0.5) & (significant_genes['Cohens_D'].abs() <= 0.8)).sum()} genes")

# Save summary
with open('deg_analysis_output/statistical_summary.txt', 'w') as f:
    f.write("DIFFERENTIAL GENE EXPRESSION ANALYSIS SUMMARY\n")
    f.write("=" * 50 + "\n\n")
    f.write(f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    f.write("STATISTICAL SUMMARY:\n")
    for key, value in summary_stats.items():
        if 'P_Value' in key:
            f.write(f"  {key.replace('_', ' ')}: {value:.2e}\n")
        elif 'Percent' in key:
            f.write(f"  {key.replace('_', ' ')}: {value:.2f}%\n")
        else:
            f.write(f"  {key.replace('_', ' ')}: {value:.2f}\n" if isinstance(value, float) else f"  {key.replace('_', ' ')}: {value}\n")

    f.write(f"\nTOP 10 UPREGULATED GENES:\n")
    top_up = deg_df[(deg_df['Regulation'] == 'Upregulated')].head(10)
    for _, gene in top_up.iterrows():
        f.write(f"  {gene['Gene']}: FC={gene['Fold_Change']:.2f}, p={gene['T_P_Adjusted']:.2e}\n")

    f.write(f"\nTOP 10 DOWNREGULATED GENES:\n")
    top_down = deg_df[(deg_df['Regulation'] == 'Downregulated')].head(10)
    for _, gene in top_down.iterrows():
        f.write(f"  {gene['Gene']}: FC={gene['Fold_Change']:.2f}, p={gene['T_P_Adjusted']:.2e}\n")

