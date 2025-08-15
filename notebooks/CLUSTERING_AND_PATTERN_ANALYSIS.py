#@title **CLUSTERING AND PATTERN ANALYSIS**


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans

print("\nClustering and Pattern Analysis...")

# Select genes (prefer significant ones)
if 'deg_df' not in globals():
    raise RuntimeError("deg_df not found. Run statistical testing first.")

sig_genes = deg_df[deg_df.get('Significant', False)]['Gene'].tolist()
if len(sig_genes) < 10:
    print(" Using top 100 genes by p-value for clustering")
    sig_genes = deg_df.sort_values('T_P_Adjusted', na_position='last').head(100)['Gene'].tolist()

# Ensure genes exist in expression matrix
sig_genes = [g for g in sig_genes if g in expression_log2.index]
if len(sig_genes) == 0:
    raise RuntimeError("No selected genes found in expression_log2. Cannot cluster.")

print(f"Performing clustering analysis on {len(sig_genes)} genes...")

# Extract expression data for significant genes
samples_order = treated_samples + control_samples
clustering_data = expression_log2.loc[sig_genes, samples_order].astype(float).copy()

# Clean data: replace inf, drop all-NaN rows/cols
clustering_data.replace([np.inf, -np.inf], np.nan, inplace=True)
clustering_data.dropna(axis=0, how='all', inplace=True)  # drop genes with all NaN
clustering_data.dropna(axis=1, how='all', inplace=True)  # drop samples with all NaN

# Drop genes with zero variance (they break scaling & clustering)
row_std = clustering_data.std(axis=1, ddof=1)
nonconstant_genes = row_std[row_std > 0].index.tolist()
clustering_data = clustering_data.loc[nonconstant_genes]

if clustering_data.shape[0] < 2:
    raise RuntimeError("Too few non-constant genes remain for clustering after cleaning.")

# Row-wise z-score (standardize each gene across samples)
clustering_z = clustering_data.sub(clustering_data.mean(axis=1), axis=0).div(
    clustering_data.std(axis=1).replace(0, np.nan), axis=0
).fillna(0)

# K-means clustering to identify expression patterns
n_clusters = 4
n_genes_available = clustering_z.shape[0]
if n_genes_available < n_clusters:
    print(f" fewer genes ({n_genes_available}) than requested clusters ({n_clusters}), reducing clusters.")
    n_clusters = max(1, n_genes_available)

kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
# fit_predict on the matrix genes x features (samples)
cluster_labels = kmeans.fit_predict(clustering_z.values)

# Add cluster information to results
cluster_df = pd.DataFrame({
    'Gene': clustering_z.index.tolist(),
    'Cluster': cluster_labels
})

# Analyze cluster patterns using original log2 expression values
print("\nCluster Analysis:")
for cluster_id in range(n_clusters):
    cluster_genes = cluster_df[cluster_df['Cluster'] == cluster_id]['Gene'].tolist()
    cluster_size = len(cluster_genes)
    if cluster_size == 0:
        print(f"  Cluster {cluster_id}: 0 genes")
        continue

    cluster_data_orig = clustering_data.loc[cluster_genes]  # original (non-z) values
    treated_mean = cluster_data_orig[treated_samples].mean(axis=1).mean()
    control_mean = cluster_data_orig[control_samples].mean(axis=1).mean()

    print(f"  Cluster {cluster_id}: {cluster_size} genes")
    print(f"    Average SCFA expression: {treated_mean:.3f}")
    print(f"    Average Control expression: {control_mean:.3f}")
    print(f"    Pattern: {'Higher in SCFA' if treated_mean > control_mean else 'Higher in Control'}")

# Visualize cluster patterns
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Gene Expression Clusters: Average Patterns', fontsize=16)
axes = axes.flatten()

for cluster_id in range(n_clusters):
    ax = axes[cluster_id]
    cluster_genes = cluster_df[cluster_df['Cluster'] == cluster_id]['Gene'].tolist()
    if not cluster_genes:
        ax.set_title(f'Cluster {cluster_id} (0 genes)')
        ax.axis('off')
        continue

    cluster_data_orig = clustering_data.loc[cluster_genes]

    # Plot individual gene trajectories (thin gray lines)
    for gene in cluster_genes[:50]:  # limit number of thin lines for clarity
        control_vals = cluster_data_orig.loc[gene, control_samples].values
        treated_vals = cluster_data_orig.loc[gene, treated_samples].values
        ax.plot([0, 1], [np.mean(control_vals), np.mean(treated_vals)], color='gray', alpha=0.25, linewidth=0.6)

    # Plot average pattern (thick line)
    avg_control = cluster_data_orig[control_samples].mean(axis=1).mean()
    avg_treated = cluster_data_orig[treated_samples].mean(axis=1).mean()
    ax.plot([0, 1], [avg_control, avg_treated], color='red', linewidth=3, marker='o', markersize=8)

    ax.set_xticks([0, 1])
    ax.set_xticklabels(['Control', 'SCFA'])
    ax.set_ylabel('Average Log2 Expression')
    ax.set_title(f'Cluster {cluster_id} ({len(cluster_genes)} genes)')
    ax.grid(True, alpha=0.3)

# Hide any unused axes (when n_clusters < 4)
for i in range(n_clusters, 4):
    axes[i].axis('off')

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig('deg_analysis_output/expression_clusters.png', dpi=300, bbox_inches='tight')
plt.show()

