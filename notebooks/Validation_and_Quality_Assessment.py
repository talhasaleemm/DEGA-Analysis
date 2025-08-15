#@title **Final Validation and Quality Assessment**


print("\nFinal Validation and Quality Assessment...")

# 1. Cross-validation of statistical tests
print("Cross-validating statistical results...")
consistent_results = deg_df['T_Significant'] & deg_df['U_Significant']
print(f"  Genes significant in both parametric and non-parametric tests: {consistent_results.sum()}")

# 2. Power analysis estimation
print("Estimating statistical power...")
effect_sizes = deg_df[deg_df['Significant']]['Cohens_D'].abs()
if len(effect_sizes) > 0:
    large_effects = (effect_sizes > 0.8).sum()
    medium_effects = ((effect_sizes > 0.5) & (effect_sizes <= 0.8)).sum()
    small_effects = ((effect_sizes > 0.2) & (effect_sizes <= 0.5)).sum()

    print(f"  Large effect sizes (d > 0.8): {large_effects}")
    print(f"  Medium effect sizes (0.5 < d ≤ 0.8): {medium_effects}")
    print(f"  Small effect sizes (0.2 < d ≤ 0.5): {small_effects}")

# 3. Quality assessment plot
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Quality Assessment: Statistical Analysis Validation', fontsize=16)

# P-value distribution
axes[0, 0].hist(deg_df['T_P_Value'], bins=50, alpha=0.7, color='skyblue', edgecolor='black')
axes[0, 0].set_xlabel('Raw P-values')
axes[0, 0].set_ylabel('Frequency')
axes[0, 0].set_title('P-value Distribution')
axes[0, 0].grid(True, alpha=0.3)

# Effect size distribution
axes[0, 1].hist(deg_df['Cohens_D'], bins=50, alpha=0.7, color='lightgreen', edgecolor='black')
axes[0, 1].axvline(x=0, color='red', linestyle='--', alpha=0.7)
axes[0, 1].set_xlabel("Cohen's D (Effect Size)")
axes[0, 1].set_ylabel('Frequency')
axes[0, 1].set_title('Effect Size Distribution')
axes[0, 1].grid(True, alpha=0.3)

# Correlation between statistical tests
axes[1, 0].scatter(deg_df['T_P_Adjusted'], deg_df['U_P_Adjusted'], alpha=0.6, s=20)
axes[1, 0].plot([0, 1], [0, 1], 'r--', alpha=0.7)
axes[1, 0].set_xlabel('T-test Adjusted P-value')
axes[1, 0].set_ylabel('Mann-Whitney Adjusted P-value')
axes[1, 0].set_title('Correlation Between Statistical Tests')
axes[1, 0].grid(True, alpha=0.3)

# Expression difference vs significance
sig_genes_plot = deg_df[deg_df['Significant']]
non_sig_genes_plot = deg_df[~deg_df['Significant']]

axes[1, 1].scatter(non_sig_genes_plot['Log2_Fold_Change'], -np.log10(non_sig_genes_plot['T_P_Adjusted']),
                  alpha=0.5, s=20, color='gray', label='Not Significant')
axes[1, 1].scatter(sig_genes_plot['Log2_Fold_Change'], -np.log10(sig_genes_plot['T_P_Adjusted']),
                  alpha=0.8, s=30, color='red', label='Significant')
axes[1, 1].set_xlabel('Log2 Fold Change')
axes[1, 1].set_ylabel('-Log10 Adjusted P-value')
axes[1, 1].set_title('Significance vs Effect Size')
axes[1, 1].legend()
axes[1, 1].grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('deg_analysis_output/quality_assessment.png', dpi=300, bbox_inches='tight')
plt.show()
