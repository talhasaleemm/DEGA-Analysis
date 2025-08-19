---
# notebooks/ — Notebooks and recommended execution order

This folder contains Colab-ready notebooks extracted from the original Colab analysis (`DEGA.ipynb`). Each notebook corresponds to a top-level section (a heading) from the original notebook. 

> **Important:** Notebooks are exploratory and contain narrative explanations. When reproducing results, follow the recommended run order below. Run *each notebook from the first cell to the last cell* unless instructed otherwise inside the notebook.

---

## Recommended execution order
Run notebooks in the following order — this sequence preserves data flow (downloads → preprocessing → EDA → DE testing → visualizations → export):

 **01_INSTALL AND IMPORT LIBRARIES.ipynb**  
   - Run this first to install required packages in Colab and import modules. Contains `pip` install commands (if running in a fresh Colab runtime) and base imports.
 **02_Robust_Download_and_Parse_GSE200309.ipynb**  
   - Download the GEO Series (GSE200309) and parse count tables + metadata using `GEOparse`. If you already have counts on disk, skip the download and run the loading cells.
 **03_DOWNLOAD_AND_LOAD_DATA.ipynb**  
   - Load counts into memory, run basic checks (missing values), and persist filtered counts to  `output/` if needed.
 **04_Define_Sample_Groups.ipynb**  
   - Define sample groups (treated vs control) and create the `sample_metadata.csv`. This notebook produces the `treated_samples` and `control_samples` variables that later notebooks expect.
 **05_EXPLORATORY_DATA_ANALYSIS(EDA)_PREPARATION.ipynb**  
   - Filtering, gene-sum thresholds, sample QC, and preliminary plots (correlation heatmap, mean-variance). Use `src/qc.py` and `src/normalization.py` helpers where appropriate.
 **06_EXPLORATORY_DATA_ANALYSIS_(ROBUST_VERSION).ipynb**  
   - Robust EDA: PCA, PCA plots, additional normalization choices and exploratory plots. Save intermediate `expression_log2` and `expression_filtered` here.
 **07_ADVANCED_VISUALIZATIONS.ipynb**  
   - Volcano, MA, heatmaps, boxplots of top genes. Uses `visualization.py` functions.
 **08_CLUSTERING_AND_PATTERN_ANALYSIS.ipynb**  
   - K-means / hierarchical clustering of differentially expressed genes; cluster visualizations and summary plots.
 **09_COMPREHENSIVE_STATISTICAL_SUMMARY.ipynb**  
   - Generate `statistical_summary.txt`, lists of top up/down genes, and create `supplementary_all_genes_analysis.csv`.
 **10_Final_Validation_and_Quality_Assessment.ipynb**  
    - Quality assessment plots, Q-Q plots, residual checks, and final validation figures.
 **11_COMPREHENSIVE_TREATMENT_EFFECT_DEMONSTRATION.ipynb**  
    - Additional treatment-effect visualizations and more robust comparisons (sensitivity checks).
 **12_FINAL_SUMMARY_AND_EXPORT.ipynb**  
    - Write out `publication_ready_results.csv`, `comprehensive_deg_results.csv`, and create the `output/` files.
 **13_Download_Output.ipynb**  
    - Package `output/` into `.zip` and download from Colab (final step for sharing results).


