#@title **EXPLORATORY DATA ANALYSIS(EDA) PREPARATION**

# SAMPLE CORRELATION HEATMAP

# Compute sample-to-sample correlation matrix
corr_matrix = counts_df.corr(method='pearson')

# Plot the heatmap
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, cmap='coolwarm', vmin=0.5, vmax=1.0, square=True,
            xticklabels=sample_ids, yticklabels=sample_ids)
plt.title("Sample Correlation Heatmap")
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()
plt.savefig('sample_correlation_heatmap.png', dpi=300)
plt.close()
# MEDIAN EXPRESSION PER SAMPLE

# Compute median gene count for each sample
medians = counts_df.median(axis=0)

# Prepare colors for SCFA vs Control
palette = {'SCFA': 'salmon', 'Control': 'skyblue'}
colors = [palette[g] for g in group_labels]

# Bar plot of median expression
plt.figure(figsize=(12, 6))
plt.bar(sample_ids, medians, color=colors)
plt.xlabel("Sample")
plt.ylabel("Median Gene Expression (raw counts)")
plt.xticks(rotation=90)
plt.title("Median Gene Expression per Sample")
# Create custom legend
import matplotlib.patches as mpatches
scfa_patch = mpatches.Patch(color=palette['SCFA'], label='SCFA-treated')
ctrl_patch = mpatches.Patch(color=palette['Control'], label='Control')
plt.legend(handles=[scfa_patch, ctrl_patch])
plt.tight_layout()
plt.savefig('median_expression_barplot.png', dpi=300)
plt.close()
# PCA PLOT

# Log-transform the counts (+1 to avoid log(0))
log_counts = np.log2(counts_df + 1)

# PCA on samples (we transpose so rows = samples)
pca = PCA(n_components=2)
pca_result = pca.fit_transform(log_counts.T)

# Percentage of variance explained by PCs
explained = pca.explained_variance_ratio_ * 100
pc_df = pd.DataFrame({
    'PC1': pca_result[:, 0],
    'PC2': pca_result[:, 1],
    'sample': sample_ids,
    'group': group_labels
})

# Plot PCA
plt.figure(figsize=(8, 6))
sns.scatterplot(x='PC1', y='PC2', hue='group', style='group',
                palette=palette, data=pc_df, s=100)
plt.xlabel(f"PC1 ({explained[0]:.1f}% variance)")
plt.ylabel(f"PC2 ({explained[1]:.1f}% variance)")
plt.title("PCA of Log-transformed Counts")
plt.legend(title='Group')
plt.tight_layout()
plt.savefig('PCA_SCFA_vs_Control.png', dpi=300)
plt.close()
# MEAN-VARIANCE SCATTER PLOT

# Compute gene means and variances
gene_means = counts_df.mean(axis=1)
gene_vars = counts_df.var(axis=1)

# Log2-transform for plotting (add 1 to avoid log(0))
log_means = np.log2(gene_means + 1)
log_vars = np.log2(gene_vars + 1)

plt.figure(figsize=(8, 6))
plt.scatter(log_means, log_vars, alpha=0.5, s=10)
plt.xlabel("Log2(Mean Expression + 1)")
plt.ylabel("Log2(Variance)")
plt.title("Mean–Variance Relationship")
plt.tight_layout()
plt.savefig('mean_variance_scatter.png', dpi=300)
plt.close()

