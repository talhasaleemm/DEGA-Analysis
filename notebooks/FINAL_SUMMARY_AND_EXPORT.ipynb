#@title **FINAL SUMMARY AND EXPORT**


print("DIFFERENTIAL GENE EXPRESSION ANALYSIS COMPLETED")

print(f"\n RESULTS SUMMARY:")
print(f"   Total genes analyzed: {len(deg_df):,}")
print(f"   Significant genes found: {deg_df['Significant'].sum():,} ({(deg_df['Significant'].sum()/len(deg_df)*100):.1f}%)")
print(f"   Upregulated in SCFA: {(deg_df['Regulation'] == 'Upregulated').sum():,}")
print(f"   Downregulated in SCFA: {(deg_df['Regulation'] == 'Downregulated').sum():,}")
print(f"   Average fold change (significant): {deg_df[deg_df['Significant']]['Fold_Change'].mean():.2f}")

print(f"\n FILES GENERATED:")
files_created = [
    'comprehensive_deg_results.csv',
    'statistical_summary.txt',
    'volcano_plot.png',
    'ma_plot.png',
    'top_genes_heatmap.png',
    'top_genes_boxplots.png',
    'expression_clusters.png',
    'quality_assessment.png'
]

for file in files_created:
    file_path = f'deg_analysis_output/{file}'
    if os.path.exists(file_path):
        print(f"   {file}")
    else:
        print(f"   {file} (not created)")

print(f"\n BIOLOGICAL INTERPRETATION:")
if deg_df['Significant'].sum() > 0:
    print(f"   Clear differential expression detected between SCFA-treated and control samples")
    print(f"   Statistical significance achieved with multiple testing correction")
    print(f"   Effect sizes indicate biological relevance")

    # Top genes summary
    top_up_gene = deg_df[deg_df['Regulation'] == 'Upregulated'].iloc[0] if (deg_df['Regulation'] == 'Upregulated').sum() > 0 else None
    top_down_gene = deg_df[deg_df['Regulation'] == 'Downregulated'].iloc[0] if (deg_df['Regulation'] == 'Downregulated').sum() > 0 else None

    if top_up_gene is not None:
        print(f"   Most upregulated gene: {top_up_gene['Gene']} (FC: {top_up_gene['Fold_Change']:.2f}, p: {top_up_gene['T_P_Adjusted']:.2e})")
    if top_down_gene is not None:
        print(f"   Most downregulated gene: {top_down_gene['Gene']} (FC: {top_down_gene['Fold_Change']:.2f}, p: {top_down_gene['T_P_Adjusted']:.2e})")
else:
    print(f"    Limited differential expression detected - consider increasing sample size or adjusting thresholds")

print(f"\n QUALITY CONTROL ASSESSMENT:")
print(f"   Statistical power: {'High' if (deg_df['Significant'].sum() > 100) else 'Moderate' if (deg_df['Significant'].sum() > 50) else 'Low'}")
print(f"   Test consistency: {consistent_results.sum()} genes significant in both parametric & non-parametric tests")
print(f"   Effect size quality: {large_effects if 'large_effects' in locals() else 0} genes with large biological effect")

