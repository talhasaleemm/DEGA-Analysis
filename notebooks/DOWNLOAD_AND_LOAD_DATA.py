#@title  **DOWNLOAD AND LOAD DATA**
# Install and import required packages
# pandas, numpy, matplotlib, seaborn, scikit-learn are usually pre-installed
!pip install GEOparse

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import sys

sns.set(style='whitegrid')  # set seaborn style for plots



# Construct the URL for the gene-level count matrix (from GEO supplemental files).
# GSE200309 is under series GSE200nnn on the NCBI FTP server.
counts_url = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE200nnn/GSE200309/"
    "suppl/GSE200309_TxImport.GeneLevel.counts.GEO.txt.gz"
)
try:
    counts_df = pd.read_csv(counts_url, sep='\t', compression='gzip', index_col=0)
except Exception as e:
    sys.exit(f"Error downloading or reading count data: {e}")



# Verify we have the expected number of samples
expected_samples = 21
actual_samples = counts_df.shape[1]
assert actual_samples == expected_samples, f"Expected {expected_samples} samples, found {actual_samples}"

# Check for missing values and drop genes with any NA
if counts_df.isnull().values.any():
    print("Warning: Missing values found; dropping incomplete genes.")
    counts_df.dropna(axis=0, inplace=True)

# FILTER LOW-EXPRESSION GENES

# Remove genes with very low total counts (sum across all samples ≤ 10).
# This follows DESeq2 recommendations for pre-filtering:contentReference[oaicite:4]{index=4}.
gene_counts_sum = counts_df.sum(axis=1)
counts_df = counts_df[gene_counts_sum > 10]
assert counts_df.shape[0] > 0, "No genes left after filtering low counts."

print(f"Data loaded with {counts_df.shape[0]} genes and {counts_df.shape[1]} samples after filtering.")

