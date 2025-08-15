#@title **ADVANCED VISUALIZATIONS**
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

print("\nCreating Advanced Visualizations...")

# Defensive checks
if 'deg_df' not in globals():
    raise RuntimeError("deg_df not found. Run the statistical testing block first.")
if deg_df.shape[0] == 0:
    raise RuntimeError("deg_df is empty. No visualization possible.")

# Ensure thresholds exist (fall back to defaults if not)
p_threshold = globals().get('p_threshold', 0.05)
log2fc_threshold = globals().get('log2fc_threshold', np.log2(1.5))

# --- Volcano Plot ---
plt.figure(figsize=(12, 8))
# Protect against zeros or invalid p-values
tp_adj = deg_df.get('T_P_Adjusted', pd.Series(np.ones(len(deg_df)), index=deg_df.index)).astype(float)

# Replace zeros and negative/NaN with the smallest positive float > 0 for plotting
min_nonzero = tp_adj[tp_adj > 0].min() if (tp_adj > 0).any() else np.finfo(float).tiny
tp_adj_safe = tp_adj.fillna(min_nonzero).replace(0, min_nonzero)
x = deg_df['Log2_Fold_Change'].astype(float)
y = -np.log10(tp_adj_safe)

# Color points based on significance and fold change (defensively)
colors = []
for _, row in deg_df.iterrows():
    sig = bool(row.get('Significant', False))
    log2fc = float(row.get('Log2_Fold_Change', 0.0))
    if sig:
        if log2fc > log2fc_threshold:
            colors.append('red')       # Upregulated
        elif log2fc < -log2fc_threshold:
            colors.append('blue')      # Downregulated
        else:
            colors.append('orange')    # Significant but low FC
    else:
        colors.append('gray')         # Not significant

plt.scatter(x, y, c=colors, alpha=0.6, s=20)

# Add threshold lines
plt.axvline(x=log2fc_threshold, color='black', linestyle='--', alpha=0.5)
plt.axvline(x=-log2fc_threshold, color='black', linestyle='--', alpha=0.5)
plt.axhline(y=-np.log10(p_threshold), color='black', linestyle='--', alpha=0.5)

plt.xlabel('Log2 Fold Change (SCFA vs Control)')
plt.ylabel('-Log10 Adjusted P-value')
plt.title('Volcano Plot: Differential Gene Expression\nSCFA-Treated vs Control')

# Fixed legend parentheses and counts
legend_elements = [
    Patch(facecolor='red', label=f'Upregulated (n={sum([c=="red" for c in colors])})'),
    Patch(facecolor='blue', label=f'Downregulated (n={sum([c=="blue" for c in colors])})'),
    Patch(facecolor='gray', label=f'Not Significant (n={sum([c=="gray" for c in colors])})')
]
plt.legend(handles=legend_elements, loc='upper right')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('deg_analysis_output/volcano_plot.png', dpi=300, bbox_inches='tight')
plt.show()

# --- MA Plot ---
plt.figure(figsize=(12, 8))
A = (deg_df['Treated_Mean'].astype(float) + deg_df['Control_Mean'].astype(float)) / 2
M = deg_df['Log2_Fold_Change'].astype(float)

plt.scatter(A, M, c=colors, alpha=0.6, s=20)
plt.axhline(y=0, color='black', linestyle='-', alpha=0.8)
plt.axhline(y=log2fc_threshold, color='red', linestyle='--', alpha=0.5)
plt.axhline(y=-log2fc_threshold, color='red', linestyle='--', alpha=0.5)

plt.xlabel('Average Expression (Log2)')
plt.ylabel('Log2 Fold Change (SCFA vs Control)')
plt.title('MA Plot: Expression vs Fold Change')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('deg_analysis_output/ma_plot.png', dpi=300, bbox_inches='tight')
plt.show()

# --- Heatmap of Top Differentially Expressed Genes ---
print("Creating heatmap of top differentially expressed genes...")

# Pick top genes (prefer significant ones)
top_genes = deg_df[deg_df.get('Significant', False)].head(50)
if len(top_genes) < 10:
    print("Few significant genes found, using top 50 genes by p-value")
    top_genes = deg_df.sort_values('T_P_Adjusted', na_position='last').head(50)

selected_genes = top_genes['Gene'].tolist()
# Make sure genes exist in expression_log2
available_genes = [g for g in selected_genes if g in expression_log2.index]
if len(available_genes) == 0:
    raise RuntimeError("None of the selected top genes exist in expression_log2. Cannot make heatmap.")

heatmap_data = expression_log2.loc[available_genes, treated_samples + control_samples].astype(float).copy()

# Clean: replace inf with nan, drop rows/cols that are all nan
heatmap_data.replace([np.inf, -np.inf], np.nan, inplace=True)
heatmap_data.dropna(axis=0, how='all', inplace=True)
heatmap_data.dropna(axis=1, how='all', inplace=True)

# Drop genes with zero variance (they break z-scoring & clustering)
row_std = heatmap_data.std(axis=1, ddof=1)
nonconstant_rows = row_std[row_std > 0].index.tolist()
heatmap_data = heatmap_data.loc[nonconstant_rows]

# If still NaNs: fill by row mean
heatmap_data = heatmap_data.apply(lambda r: r.fillna(r.mean()), axis=1)

# If too few genes/samples remain for clustering, fall back to a simple heatmap
do_clustermap = (heatmap_data.shape[0] >= 2 and heatmap_data.shape[1] >= 2)

if do_clustermap:
    # Row z-score (standardize rows) but avoid divide-by-zero (we dropped zero-variance rows above)
    heatmap_z = heatmap_data.sub(heatmap_data.mean(axis=1), axis=0).div(heatmap_data.std(axis=1).replace(0, np.nan), axis=0).fillna(0)

    sample_colors = ['red' if s in treated_samples else 'blue' for s in heatmap_z.columns]
    try:
        # clustermap creates its own figure, don't call plt.figure() before it
        cg = sns.clustermap(heatmap_z,
                            col_colors=sample_colors,
                            cmap='RdBu_r',
                            center=0,
                            z_score=None,          # already z-scored manually
                            figsize=(16, 12),
                            method='average',
                            metric='euclidean',
                            dendrogram_ratio=(.1, .2),
                            cbar_kws={'label': 'Z-score'})
        plt.suptitle(f'Heatmap: Top {heatmap_z.shape[0]} Differentially Expressed Genes\n(Red=SCFA, Blue=Control)', y=1.02, fontsize=14)
        cg.savefig('deg_analysis_output/top_genes_heatmap.png')
        plt.show()
    except Exception as e:
        print("Clustermap failed (falling back to simple heatmap). Error:", str(e))
        do_clustermap = False

if not do_clustermap:
    # Simple heatmap fallback
    plt.figure(figsize=(12, max(4, heatmap_data.shape[0] * 0.25)))
    sns.heatmap(heatmap_data, cmap='RdBu_r', center=0)
    plt.title(f'Heatmap (no clustering): Top {heatmap_data.shape[0]} Genes')
    plt.xlabel('Samples')
    plt.ylabel('Genes')
    plt.tight_layout()
    plt.savefig('deg_analysis_output/top_genes_heatmap_simple.png', dpi=300, bbox_inches='tight')
    plt.show()

# --- Box plots for top genes ---
print("Creating box plots for top differentially expressed genes...")

top_12_genes = deg_df[deg_df.get('Significant', False)].head(12)
if len(top_12_genes) < 12:
    top_12_genes = deg_df.sort_values('T_P_Adjusted', na_position='last').head(12)

fig, axes = plt.subplots(3, 4, figsize=(20, 15))
fig.suptitle('Expression Levels: Top 12 Differentially Expressed Genes', fontsize=16)
axes = axes.flatten()

for i, (_, gene_info) in enumerate(top_12_genes.iterrows()):
    if i >= 12:
        break
    gene = gene_info['Gene']
    if gene not in expression_log2.index:
        axes[i].text(0.5, 0.5, f'{gene}\n(not found)', ha='center')
        axes[i].set_axis_off()
        continue

    treated_values = expression_log2.loc[gene, treated_samples].values
    control_values = expression_log2.loc[gene, control_samples].values

    data = [control_values, treated_values]
    labels = ['Control', 'SCFA']

    bp = axes[i].boxplot(data, labels=labels, patch_artist=True)
    try:
        bp['boxes'][0].set_facecolor('lightblue')
        bp['boxes'][1].set_facecolor('lightcoral')
    except Exception:
        pass

    pval = gene_info.get('T_P_Adjusted', np.nan)
    fc = gene_info.get('Fold_Change', np.nan)
    axes[i].set_title(f'{gene}\nFC: {fc:.2f}, p: {pval:.2e}')
    axes[i].set_ylabel('Log2 Expression')
    axes[i].grid(True, alpha=0.3)

# hide unused subplots
for j in range(i+1, 12):
    axes[j].set_axis_off()

plt.tight_layout()
plt.savefig('deg_analysis_output/top_genes_boxplots.png', dpi=300, bbox_inches='tight')
plt.show()

print("All visualizations created (saved into deg_analysis_output/).")


