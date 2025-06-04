# Guanine Oxidation and Gene Expression in Stress and Circadian Disruption Analysis

A bioinformatics pipeline for analyzing guanine oxidation and gene expression data using statistical methods (ANOVA), correlation analysis, and machine learning approaches (PCA, t-SNE, UMAP, MOFA) efficiently. Results are organized by genomic regions (CpG islands, gene bodies, promoters) and analysis type.

This project aims to
1. Identify genome-wide changes in guanine oxidation associated with chronic stress level and time of day.
2. Determine whether guanine oxidation has an epigenetic role in gene expression regulation.

## Project Structure

```
.
├── .gitignore
├── README.md
├── bash/
│   ├── intersect_cpg.sh
│   ├── intersect_genebodies.sh
│   └── intersect_promoters.sh
├── images/
│   ├── anova_results/
│   │   ├── bin1000/
│   │   │   ├── combined_significance_heatmap.png
│   │   │   ├── combined_significance_heatmap_zscore.png
│   │   │   ├── combined_significance_heatmap_zscore_clustered.png
│   │   │   ├── heatmap_legend.png
│   │   │   ├── indiv/                                    # Individual plots for each significant bin
│   │   │   │   ├── interaction_top1-11_chr*.png        # (11 files)
│   │   │   │   ├── timepoint_top1-47_chr*.png          # (47 files)
│   │   │   │   └── treatment_top1-20_chr*.png          # (20 files)
│   │   │   ├── manhattan_plot_interaction.png
│   │   │   ├── manhattan_plot_timepoint.png
│   │   │   └── manhattan_plot_treatment.png
│   │   ├── bin10000/ 
│   │   │   ├── combined_significance_heatmap.png
│   │   │   ├── combined_significance_heatmap_zscore.png
│   │   │   ├── combined_significance_heatmap_zscore_clustered.png
│   │   │   ├── heatmap_legend.png
│   │   │   ├── indiv/                                    # Individual plots for each significant bin
│   │   │   │   └── timepoint_top1-5_chr*.png           # (5 files)
│   │   │   ├── manhattan_plot_interaction.png
│   │   │   ├── manhattan_plot_timepoint.png
│   │   │   └── manhattan_plot_treatment.png
│   │   ├── bin100000/
│   │   │   ├── combined_significance_heatmap.png
│   │   │   ├── combined_significance_heatmap_zscore.png
│   │   │   ├── combined_significance_heatmap_zscore_clustered.png
│   │   │   ├── heatmap_legend.png
│   │   │   ├── indiv/                                    # Individual plots for each significant bin
│   │   │   │   ├── timepoint_top1-10_chr*.png          # (11 files)
│   │   │   │   └── treatment_top1_chr*.png             # (1 file)
│   │   │   ├── manhattan_plot_interaction.png
│   │   │   ├── manhattan_plot_timepoint.png
│   │   │   └── manhattan_plot_treatment.png
│   │   ├── cpg/
│   │   │   ├── combined_significance_heatmap.png
│   │   │   ├── combined_significance_heatmap_zscore.png
│   │   │   ├── combined_significance_heatmap_zscore_clustered.png
│   │   │   ├── heatmap_legend.png
│   │   │   ├── indiv/                                    # Individual plots for each significant island
│   │   │   │   ├── interaction_top1_chr*.png           # (1 file)
│   │   │   │   ├── timepoint_top1-41_chr*.png          # (42 files)
│   │   │   │   └── treatment_top1-2_chr*.png           # (2 files)
│   │   │   ├── manhattan_plot_interaction.png
│   │   │   ├── manhattan_plot_timepoint.png
│   │   │   └── manhattan_plot_treatment.png
│   │   ├── gene_bodies/
│   │   │   ├── manhattan_plot_interaction.png
│   │   │   ├── manhattan_plot_timepoint.png
│   │   │   └── manhattan_plot_treatment.png
│   │   └── promoters/
│   │       ├── combined_significance_heatmap.png
│   │       ├── combined_significance_heatmap_zscore.png
│   │       ├── combined_significance_heatmap_zscore_clustered.png
│   │       ├── heatmap_legend.png
│   │       ├── indiv/                                    # Individual plots for each significant gene
│   │       │   └── timepoint_top1-4_ENSMUSG*.png       # (4 files)
│   │       ├── manhattan_plot_interaction.png
│   │       ├── manhattan_plot_timepoint.png
│   │       └── manhattan_plot_treatment.png
│   ├── expr_global/                                      # Global correlation analysis between oxidation and expression data
│   │   ├── cpg/
│   │   │   ├── raw_data_full.png
│   │   │   ├── raw_data_zoomed.png
│   │   │   └── spearman_correlation_by_sample.png
│   │   ├── gene_bodies/
│   │   │   ├── raw_data_full.png
│   │   │   ├── raw_data_zoomed.png
│   │   │   └── spearman_correlation_by_sample.png
│   │   └── promoters/
│   │       ├── raw_data_full.png
│   │       ├── raw_data_zoomed.png
│   │       └── spearman_correlation_by_sample.png
│   ├── gene_correlation_plots/                          # Local correlation analysis between oxidation and expression data for each significant gene
│   │   ├── ENSMUSG000000*.png                          # (16 files)
│   │   └── correlation_results.csv
│   ├── mofa_results/                                     # Global correlation analysis between oxidation and expression data using MOFA
│   │   ├── cpg/
│   │   │   ├── group_effects_analysis.png
│   │   │   └── mofa_analysis_overview.png
│   │   ├── gene_bodies/
│   │   │   ├── group_effects_analysis.png
│   │   │   └── mofa_analysis_overview.png
│   │   └── promoters/
│   │       ├── group_effects_analysis.png
│   │       └── mofa_analysis_overview.png
│   ├── pca_results/
│   │   ├── bin1000/
│   │   │   ├── pca_2d_plot_bin1000.png
│   │   │   ├── pca_3d_plot_bin1000.png
│   │   │   └── pca_scree_plot_bin1000.png
│   │   ├── bin10000/
│   │   │   ├── pca_2d_plot_bin10000.png
│   │   │   ├── pca_3d_plot_bin10000.png
│   │   │   └── pca_scree_plot_bin10000.png
│   │   ├── bin100000/
│   │   │   ├── pca_2d_plot_bin100000.png
│   │   │   ├── pca_3d_plot_bin100000.png
│   │   │   └── pca_scree_plot_bin100000.png
│   │   ├── cpg/
│   │   │   ├── pca_2d_plot.png
│   │   │   ├── pca_3d_plot.png
│   │   │   └── pca_scree_plot.png
│   │   ├── gene_bodies/
│   │   │   ├── pca_2d_plot.png
│   │   │   ├── pca_3d_plot.png
│   │   │   └── pca_scree_plot.png
│   │   └── promoters/
│   │       ├── pca_2d_plot.png
│   │       ├── pca_3d_plot.png
│   │       └── pca_scree_plot.png
│   ├── rf_results/                                       # Oxidation analysis using Random Forest
│   │   ├── bin1000/
│   │   │   └── model_performance_comparison_plot.png
│   │   ├── bin10000/
│   │   │   └── model_performance_comparison_plot.png
│   │   ├── bin100000/
│   │   │   └── model_performance_comparison_plot.png
│   │   ├── cpg/
│   │   │   └── model_performance_comparison_plot.png
│   │   ├── gene_bodies/
│   │   │   └── model_performance_comparison_plot.png
│   │   └── promoters/
│   │       └── model_performance_comparison_plot.png
│   ├── tsne_results/
│   │   ├── bin1000/
│   │   │   ├── tsne_2d_perplexity_comparison_bin1000.png
│   │   │   ├── tsne_2d_plot_bin1000.png
│   │   │   └── tsne_3d_perplexity_comparison_bin1000.png
│   │   ├── bin10000/
│   │   │   ├── tsne_2d_perplexity_comparison_bin10000.png
│   │   │   ├── tsne_2d_plot_bin10000.png
│   │   │   └── tsne_3d_perplexity_comparison_bin10000.png
│   │   ├── bin100000/
│   │   │   ├── tsne_2d_perplexity_comparison_bin100000.png
│   │   │   ├── tsne_2d_plot_bin100000.png
│   │   │   └── tsne_3d_perplexity_comparison_bin100000.png
│   │   ├── cpg/
│   │   │   ├── tsne_2d_perplexity_comparison.png
│   │   │   ├── tsne_2d_plot.png
│   │   │   └── tsne_3d_perplexity_comparison.png
│   │   ├── gene_bodies/
│   │   │   ├── tsne_2d_perplexity_comparison.png
│   │   │   ├── tsne_2d_plot.png
│   │   │   └── tsne_3d_perplexity_comparison.png
│   │   └── promoters/
│   │       ├── tsne_2d_perplexity_comparison.png
│   │       ├── tsne_2d_plot.png
│   │       └── tsne_3d_perplexity_comparison.png
│   ├── umap_results/ 
│   │   ├── bin1000/
│   │   │   ├── umap_2d_final_plot_bin1000.png
│   │   │   ├── umap_mindist_comparison_bin1000.png
│   │   │   └── umap_neighbors_comparison_bin1000.png
│   │   ├── bin10000/
│   │   │   ├── umap_2d_final_plot_bin10000.png
│   │   │   ├── umap_mindist_comparison_bin10000.png
│   │   │   └── umap_neighbors_comparison_bin10000.png
│   │   ├── bin100000/
│   │   │   ├── umap_2d_final_plot_bin100000.png
│   │   │   ├── umap_mindist_comparison_bin100000.png
│   │   │   └── umap_neighbors_comparison_bin100000.png
│   │   ├── cpg/
│   │   │   ├── umap_2d_final_plot.png
│   │   │   ├── umap_mindist_comparison.png
│   │   │   └── umap_neighbors_comparison.png
│   │   ├── gene_bodies/
│   │   │   ├── umap_2d_final_plot.png
│   │   │   ├── umap_mindist_comparison.png
│   │   │   └── umap_neighbors_comparison.png
│   │   └── promoters/
│   │       ├── umap_2d_final_plot.png
│   │       ├── umap_mindist_comparison.png
│   │       └── umap_neighbors_comparison.png
│   └── xg_results/                                       # Oxidation analysis using XGBoost
│       ├── bin1000/
│       │   └── xgboost_performance_comparison_plot.png
│       ├── bin10000/
│       │   └── xgboost_performance_comparison_plot.png
│       ├── bin100000/
│       │   └── xgboost_performance_comparison_plot.png
│       ├── cpg/
│       │   └── xgboost_performance_comparison_plot.png
│       ├── gene_bodies/
│       │   └── xgboost_performance_comparison_plot.png
│       └── promoters/
│           └── xgboost_performance_comparison_plot.png
└── jupyter_notebooks/
    ├── anova_bin.ipynb                                  # Two-way ANOVA for genomic bins
    ├── anova_rest.ipynb                                 # Two-way ANOVA for cpg islands, gene bodies and promoters
    ├── binning_genome_wide_data.ipynb                   # Getting oxidation data for genomic bins
    ├── cpg_expr.ipynb                                   # Getting expression data for cpg islands
    ├── cpg_oxid.ipynb                                   # Getting oxidation data for cpg islands
    ├── dim_reduction_bin.ipynb                          # Dimension reduction for genomic bins
    ├── dim_reduction_rest.ipynb                         # Dimension reduction for cpg islands, gene bodies and promoters
    ├── expr_global.ipynb                                # Global correlation analysis using Spearman Correlation
    ├── expr_local.ipynb                                 # Local correlation analysis for significant genes
    ├── fetch_data.ipynb                                 # Converting bed files to csv files for oxidation data per sample
    ├── get_oxid.ipynb                                   # Getting oxidation data for gene bodies and promoters
    ├── mofa.ipynb                                       # Global correlation analysis using MOFA
    ├── obtaining_gene_promoter_CpG_island_coordinates.ipynb
    ├── parsing_gene_bounds.ipynb
    ├── random_forest_bin.ipynb                          # Random Forest for genomic bins
    ├── random_forest_rest.ipynb                         # Random Forest for cpg islands, gene bodies and promoters
    ├── sig_intersection.ipynb                           # Intersecting significant features (genomic bins, cpg islands) with gene bodies and promoters
    ├── xgboost_bin.ipynb                                # XGBoost for genomic bins
    └── xgboost_rest.ipynb                               # XGBoost for cpg islands, gene bodies and promoters
```

## Summary 

