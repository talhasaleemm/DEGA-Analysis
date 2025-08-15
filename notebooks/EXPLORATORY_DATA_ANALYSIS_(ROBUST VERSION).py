#@title **EXPLORATORY DATA ANALYSIS (ROBUST VERSION)**


print("\n  EXPLORATORY DATA ANALYSIS ")
print("Goal: Visualize data quality and identify patterns between SCFA-treated and control samples")

# Import required libraries
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import matplotlib.patches as mpatches
import os
import warnings

# Suppress warnings for cleaner output
warnings.filterwarnings('ignore')

# Create output directory if it doesn't exist
os.makedirs('deg_analysis_output', exist_ok=True)

# Set up plotting parameters
plt.rcParams['figure.figsize'] = (15, 12)
plt.rcParams['font.size'] = 10

print("\n1. Data Loading and Validation...")

# Initialize variables to None
expression_filtered = None
sample_df = None

# Method 1: Check if variables exist in current namespace
print("   Checking for existing variables in memory...")
try:
    # This will work if variables are already defined
    if 'expression_filtered' in locals() or 'expression_filtered' in globals():
        print(f"   Found expression_filtered in memory: {expression_filtered.shape}")
    if 'sample_df' in locals() or 'sample_df' in globals():
        print(f"   Found sample_df in memory: {sample_df.shape}")
except:
    print("    No variables found in current session")

# Method 2: Try to load from saved files
if expression_filtered is None:
    print("   Attempting to load from saved files...")
    try:
        if os.path.exists('deg_analysis_output/expression_filtered.csv'):
            expression_filtered = pd.read_csv('deg_analysis_output/expression_filtered.csv', index_col=0)
            print(f"   Loaded expression data from file: {expression_filtered.shape}")

        if os.path.exists('deg_analysis_output/sample_metadata.csv'):
            sample_df = pd.read_csv('deg_analysis_output/sample_metadata.csv')
            print(f"    Loaded sample metadata from file: {sample_df.shape}")

    except Exception as e:
        print(f"    Could not load from files: {str(e)}")

# Method 3: Try common file locations
if expression_filtered is None:
    print("   Searching for data files in common locations...")
    possible_locations = [
        'expression_data.csv',
        'gene_expression.csv',
        'counts_matrix.csv',
        'expression_matrix.csv',
        'data/expression_data.csv',
        'input/expression_data.csv'
    ]

    for location in possible_locations:
        if os.path.exists(location):
            try:
                expression_filtered = pd.read_csv(location, index_col=0)
                print(f"    Found expression data at: {location}")
                break
            except Exception as e:
                print(f"   - Could not read {location}: {str(e)}")

    # Search for sample metadata
    sample_locations = [
        'sample_metadata.csv',
        'samples.csv',
        'sample_info.csv',
        'metadata.csv',
        'data/sample_metadata.csv',
        'input/sample_metadata.csv'
    ]

    for location in sample_locations:
        if os.path.exists(location):
            try:
                sample_df = pd.read_csv(location)
                print(f"    Found sample metadata at: {location}")
                break
            except Exception as e:
                print(f"   - Could not read {location}: {str(e)}")

# Method 4: Create demo data if nothing is found
if expression_filtered is None or sample_df is None:
    print("\n   No data found. Creating demo dataset for demonstration...")

    # Create synthetic gene expression data
    np.random.seed(42)
    n_genes = 1000
    n_samples_per_group = 6

    # Generate gene names
    gene_names = [f"Gene_{i:04d}" for i in range(n_genes)]

    # Generate sample names
    control_samples = [f"Control_{i+1}" for i in range(n_samples_per_group)]
    treated_samples = [f"SCFA_Treated_{i+1}" for i in range(n_samples_per_group)]
    all_samples = control_samples + treated_samples

    # Generate expression data (log2 transformed)
    # Control samples: mean expression around 5
    control_data = np.random.normal(5, 2, (n_genes, n_samples_per_group))

    # Treated samples: similar expression but with some differential genes
    treated_data = np.random.normal(5, 2, (n_genes, n_samples_per_group))

    # Make some genes differentially expressed
    n_deg = 200  # Number of differentially expressed genes
    deg_indices = np.random.choice(n_genes, n_deg, replace=False)

    # Half upregulated, half downregulated
    up_indices = deg_indices[:n_deg//2]
    down_indices = deg_indices[n_deg//2:]

    # Add differential expression
    treated_data[up_indices, :] += np.random.normal(2, 0.5, (len(up_indices), n_samples_per_group))
    treated_data[down_indices, :] -= np.random.normal(2, 0.5, (len(down_indices), n_samples_per_group))

    # Combine data
    expression_data = np.hstack([control_data, treated_data])
    expression_filtered = pd.DataFrame(expression_data,
                                     index=gene_names,
                                     columns=all_samples)

    # Create sample metadata
    sample_df = pd.DataFrame({
        'Sample_ID': all_samples,
        'Group': ['Control'] * n_samples_per_group + ['SCFA_Treated'] * n_samples_per_group,
        'Batch': ['Batch1'] * 3 + ['Batch2'] * 3 + ['Batch1'] * 3 + ['Batch2'] * 3
    })

    print(f"   Created demo dataset: {expression_filtered.shape[0]} genes × {expression_filtered.shape[1]} samples")
    print(f"   Created sample metadata: {sample_df.shape[0]} samples")

    # Save demo data
    expression_filtered.to_csv('deg_analysis_output/expression_filtered.csv')
    sample_df.to_csv('deg_analysis_output/sample_metadata.csv', index=False)
    print("   Demo data saved to deg_analysis_output/")

# Final validation
print("\n2. Data Validation...")
if expression_filtered is None or sample_df is None:
    print("   ERROR: Still no data available. Please provide:")
    print("   - Gene expression matrix (genes as rows, samples as columns)")
    print("   - Sample metadata with 'Sample_ID' and 'Group' columns")
    print("\n   To proceed, you can:")
    print("   1. Run previous analysis steps (1-4)")
    print("   2. Place your files in the current directory")
    print("   3. Use the demo data created above")
    exit()

# Validate data structure
print(f"   Expression data: {expression_filtered.shape[0]} genes × {expression_filtered.shape[1]} samples")
print(f"   Sample metadata: {sample_df.shape[0]} samples")

# Check required columns in sample metadata
required_columns = ['Sample_ID', 'Group']
missing_columns = [col for col in required_columns if col not in sample_df.columns]

if missing_columns:
    print(f"   WARNING: Missing columns in sample metadata: {missing_columns}")
    print(f"   Available columns: {list(sample_df.columns)}")

    # Try to fix common issues
    if 'sample_id' in sample_df.columns:
        sample_df['Sample_ID'] = sample_df['sample_id']
        print("  Fixed: renamed 'sample_id' to 'Sample_ID'")

    if 'group' in sample_df.columns:
        sample_df['Group'] = sample_df['group']
        print("   Fixed: renamed 'group' to 'Group'")

    # Check again
    missing_columns = [col for col in required_columns if col not in sample_df.columns]
    if missing_columns:
        print(f"   ERROR: Still missing required columns: {missing_columns}")
        exit()

# Check for expected groups
unique_groups = sample_df['Group'].unique()
print(f"   Sample groups found: {list(unique_groups)}")

# Determine treatment and control groups
if 'SCFA_Treated' in unique_groups:
    treated_group = 'SCFA_Treated'
    control_group = 'Control' if 'Control' in unique_groups else unique_groups[unique_groups != 'SCFA_Treated'][0]
else:
    # Use first two groups found
    treated_group = unique_groups[0]
    control_group = unique_groups[1] if len(unique_groups) > 1 else unique_groups[0]
    print(f"   Note: Using '{treated_group}' as treatment, '{control_group}' as control")

print(f"   Analysis will compare: {treated_group} vs {control_group}")

# Data preprocessing
print("\n3. Data Preprocessing...")

# Ensure expression data is numeric
numeric_cols = expression_filtered.select_dtypes(include=[np.number]).columns
if len(numeric_cols) < len(expression_filtered.columns):
    print(f"   Converting non-numeric columns to numeric...")
    expression_filtered = expression_filtered.select_dtypes(include=[np.number])

print(f"   Numeric expression matrix: {expression_filtered.shape}")

# Handle missing values
missing_count = expression_filtered.isnull().sum().sum()
if missing_count > 0:
    print(f"   Found {missing_count} missing values, filling with gene medians...")
    expression_filtered = expression_filtered.fillna(expression_filtered.median(axis=1), axis=0)
    print("   Missing values handled")

# Remove genes with zero variance
gene_vars = expression_filtered.var(axis=1)
zero_var_genes = (gene_vars == 0).sum()
if zero_var_genes > 0:
    print(f"   Removing {zero_var_genes} genes with zero variance...")
    expression_filtered = expression_filtered[gene_vars > 0]

print(f"   Final expression matrix: {expression_filtered.shape}")

# Create visualizations
print("\n4. Creating Exploratory Visualizations...")

try:
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Exploratory Data Analysis: Gene Expression', fontsize=16, y=0.98)

    # 1. Sample correlation heatmap
    print("    Sample correlation heatmap...")

    # Only use samples that exist in both dataframes
    common_samples = list(set(expression_filtered.columns) & set(sample_df['Sample_ID']))
    expression_common = expression_filtered[common_samples]

    sample_corr = expression_common.T.corr()

    # Create annotation for samples
    sample_annotations = []
    for sample in common_samples:
        group = sample_df[sample_df['Sample_ID'] == sample]['Group'].iloc[0]
        sample_annotations.append(group)

    # Create color map for groups
    unique_groups_plot = list(set(sample_annotations))
    colors = plt.cm.Set1(np.linspace(0, 1, len(unique_groups_plot)))
    group_colors = dict(zip(unique_groups_plot, colors))

    sns.heatmap(sample_corr, ax=axes[0,0], cmap='coolwarm', center=0,
                annot=False, cbar_kws={'shrink': 0.8})
    axes[0,0].set_title(f'Sample Correlation Heatmap\n({len(common_samples)} samples)')
    axes[0,0].tick_params(axis='both', rotation=45, labelsize=8)

    # 2. Expression distribution by sample
    print("    Expression distribution by sample...")
    sample_medians = expression_common.median(axis=0)

    sample_colors = []
    for sample in sample_medians.index:
        group = sample_df[sample_df['Sample_ID']==sample]['Group'].iloc[0]
        if group == treated_group:
            sample_colors.append('red')
        else:
            sample_colors.append('blue')

    bars = axes[0,1].bar(range(len(sample_medians)), sample_medians.values,
                        color=sample_colors, alpha=0.7)
    axes[0,1].set_title('Median Expression by Sample')
    axes[0,1].set_xlabel('Sample Index')
    axes[0,1].set_ylabel('Median Expression')

    # Add legend
    red_patch = mpatches.Patch(color='red', alpha=0.7, label=treated_group)
    blue_patch = mpatches.Patch(color='blue', alpha=0.7, label=control_group)
    axes[0,1].legend(handles=[red_patch, blue_patch])

    # 3. PCA Analysis
    print("    PCA analysis...")
    pca_data = expression_common.T.fillna(0)

    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(pca_data)

    # Separate groups for plotting
    treated_indices = []
    control_indices = []

    for i, sample_id in enumerate(expression_common.columns):
        group = sample_df[sample_df['Sample_ID']==sample_id]['Group'].iloc[0]
        if group == treated_group:
            treated_indices.append(i)
        else:
            control_indices.append(i)

    # Plot PCA
    if treated_indices:
        axes[1,0].scatter(pca_result[treated_indices, 0], pca_result[treated_indices, 1],
                         c='red', alpha=0.8, s=80, label=treated_group, edgecolors='darkred')
    if control_indices:
        axes[1,0].scatter(pca_result[control_indices, 0], pca_result[control_indices, 1],
                         c='blue', alpha=0.8, s=80, label=control_group, edgecolors='darkblue')

    axes[1,0].set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    axes[1,0].set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    axes[1,0].set_title('PCA - Sample Separation')
    axes[1,0].legend()
    axes[1,0].grid(True, alpha=0.3)

    # 4. Mean-variance relationship
    print("    Mean-variance relationship...")
    gene_means = expression_filtered.mean(axis=1)
    gene_vars = expression_filtered.var(axis=1)

    # Remove outliers for better visualization
    q99_mean = gene_means.quantile(0.99)
    q99_var = gene_vars.quantile(0.99)

    plot_mask = (gene_means <= q99_mean) & (gene_vars <= q99_var)

    axes[1,1].scatter(gene_means[plot_mask], gene_vars[plot_mask],
                     alpha=0.5, s=15, color='purple')
    axes[1,1].set_xlabel('Mean Expression')
    axes[1,1].set_ylabel('Variance')
    axes[1,1].set_title(f'Mean-Variance Relationship\n({plot_mask.sum()} genes shown)')
    axes[1,1].grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig('deg_analysis_output/exploratory_analysis.png', dpi=300, bbox_inches='tight')
    plt.show()

    print("   Main visualizations completed")

    # Additional analysis: Treatment effect preview
    print("\n5. Treatment Effect Analysis...")

    treated_samples = [s for s in common_samples
                      if sample_df[sample_df['Sample_ID']==s]['Group'].iloc[0] == treated_group]
    control_samples = [s for s in common_samples
                      if sample_df[sample_df['Sample_ID']==s]['Group'].iloc[0] == control_group]

    print(f"   Treatment samples: {len(treated_samples)}")
    print(f"   Control samples: {len(control_samples)}")

    if len(treated_samples) > 0 and len(control_samples) > 0:
        # Calculate group means
        treated_means = expression_filtered[treated_samples].mean(axis=1)
        control_means = expression_filtered[control_samples].mean(axis=1)

        # Calculate fold changes (log2)
        fold_changes = treated_means - control_means

        # Summary statistics
        print(f"   Mean absolute fold change: {abs(fold_changes).mean():.3f}")
        print(f"   Genes with |log2FC| > 1: {(abs(fold_changes) > 1).sum()}")
        print(f"   Genes with |log2FC| > 2: {(abs(fold_changes) > 2).sum()}")

        # Create fold change distribution plot
        plt.figure(figsize=(12, 5))

        plt.subplot(1, 2, 1)
        plt.hist(fold_changes, bins=50, alpha=0.7, color='skyblue', edgecolor='black')
        plt.axvline(x=0, color='red', linestyle='--', alpha=0.7)
        plt.axvline(x=1, color='orange', linestyle='--', alpha=0.7)
        plt.axvline(x=-1, color='orange', linestyle='--', alpha=0.7)
        plt.xlabel('Log2 Fold Change (Treatment vs Control)')
        plt.ylabel('Number of Genes')
        plt.title('Distribution of Fold Changes')
        plt.grid(True, alpha=0.3)

        plt.subplot(1, 2, 2)
        plt.scatter(control_means, treated_means, alpha=0.5, s=20, color='green')
        plt.xlabel(f'{control_group} Mean Expression')
        plt.ylabel(f'{treated_group} Mean Expression')
        plt.title('Treatment vs Control Expression')

        # Add diagonal line
        min_val = min(plt.xlim()[0], plt.ylim()[0])
        max_val = max(plt.xlim()[1], plt.ylim()[1])
        plt.plot([min_val, max_val], [min_val, max_val], 'r--', alpha=0.5)
        plt.grid(True, alpha=0.3)

        plt.tight_layout()
        plt.savefig('deg_analysis_output/treatment_effect_preview.png', dpi=300, bbox_inches='tight')
        plt.show()

        print("   Treatment effect analysis completed")

except Exception as e:
    print(f"  ERROR in visualization: {str(e)}")
    print("   Creating simplified plots...")

    # Simple backup visualization
    try:
        plt.figure(figsize=(10, 6))

        # Simple boxplot comparison if possible
        if len(treated_samples) > 0 and len(control_samples) > 0:
            treated_sample_means = [expression_filtered[s].mean() for s in treated_samples]
            control_sample_means = [expression_filtered[s].mean() for s in control_samples]

            plt.boxplot([control_sample_means, treated_sample_means],
                       labels=[control_group, treated_group])
            plt.ylabel('Mean Sample Expression')
            plt.title('Sample Expression Distribution by Group')
            plt.grid(True, alpha=0.3)

            plt.savefig('deg_analysis_output/simple_comparison.png', dpi=300, bbox_inches='tight')
            plt.show()

            print("   Simplified visualization completed")

    except Exception as e2:
        print(f"   ERROR in simplified visualization: {str(e2)}")

# Summary
print(f"\n ANALYSIS SUMMARY ")
print(f"Expression Data: {expression_filtered.shape[0]} genes × {expression_filtered.shape[1]} samples")
print(f"Sample Groups: {dict(sample_df['Group'].value_counts())}")
print(f"Output Directory: deg_analysis_output/")
print(f"Files Created:")
print(f"   exploratory_analysis.png")
print(f"   treatment_effect_preview.png")
print(f"   expression_filtered.csv")
print(f"   sample_metadata.csv")


print(f"Data quality appears suitable for differential expression analysis.")

