# AutoML4.0

AutoML4R packages a cleaned and testable `auto_ml_analysis()` workflow from `自动机器学习4.0.R` as an installable R package.

## Installation

```r
# install.packages("remotes")
remotes::install_github("Lxyhasadream/AutoML4.0")
```

Some modelling backends are optional runtime dependencies and may require special installation steps, especially `catboost` and `lightgbm`.

## Usage

```r
library(AutoML4R)

result <- auto_ml_analysis(
  hub_data = hub_data,
  group = group,
  positive_class = "Disease",
  output_dir = "ML_screening_results",
  methods = aml_methods(),
  top_n = 20
)
```

`hub_data` should be a data frame or matrix with samples in rows and features in columns. `group` should contain one binary class label per sample.

The main workflow writes:

- `feature_rankings_by_method.csv`
- `consensus_feature_votes.csv`
- one selected-gene TXT file per method
- basic consensus vote and method-support figures
- per-algorithm plots under `algorithm_plots/`, including ranking plots,
  glmnet CV/path plots, random-forest OOB/importance plots, Boruta history
  plots, boosting importance plots, and LightGBM SHAP summaries

## Publication-Style Plots

After feature selection, run `automl4r_beautiful_plots()` to generate the
full set of publication-oriented figures:

```r
beautiful <- automl4r_beautiful_plots(
  result,
  hub_data = hub_data,
  group = group,
  positive_class = "Disease"
)
```

This writes additional figures to `output_dir`, including:

- a Nature-style final consensus voting overview
- consensus vote lollipop plot
- method-by-gene consensus heatmap grouped by algorithm class
- algorithm selected-gene count plot
- one selected-gene score plot per algorithm
- one selected-gene expression violin/boxplot per algorithm
- algorithm-specific diagnostic plots where model objects support them,
  including LASSO/Ridge/ElasticNet CV and coefficient-path plots,
  random-forest OOB and importance plots, and model-specific ranking plots
  for tree, boosting, kernel, information, latent-variable, and heuristic
  search methods

Key output files include:

- `NatureStyle_Consensus_feature_selection_overview_FINAL_largefont.pdf`
- `NatureStyle_Consensus_gene_vote_lollipop_largefont.pdf/.png`
- `NatureStyle_Method_gene_consensus_heatmap_class_size_ordered_largefont.pdf/.png`
- `NatureStyle_Algorithm_selected_gene_count_class_size_ordered_largefont.pdf/.png`
- `NatureStyle_algorithm_selected_gene_plots/<method>/`

## PPI Module: Gene Set to STRING Network to Hub Genes

AutoML4R provides two PPI entry points.

The recommended entry point for most users is `auto_string_ppi_analysis()`.
The input is only a gene-symbol vector. If the local STRING database files are
not present, AutoML4R downloads the required STRING v12.0 human files into
`resource/`, then maps genes to STRING protein IDs, extracts the query-query
network, runs 11 cytoHubba-style topology algorithms, and returns the strict
intersection of the top genes across all 11 methods.

```r
library(AutoML4R)

genes <- c(
  "SERPINE1", "PF4", "HMOX1", "PPBP", "UCP2",
  "CXCR4", "HBB", "PTGS1", "GPX1", "HK2",
  "GZMB", "PRF1", "PRRT2"
)

ppi_result <- auto_string_ppi_analysis(
  genes = genes,
  resource_dir = "resource",
  output_dir = "STRING_PPI_analysis_results",
  species_id = "9606",
  score_threshold = 400,
  top_n = 10,
  download_resources = TRUE
)

# Gene-level STRING network used for topology analysis
ppi_result$edge_df

# Strict intersection of top-10 genes across all 11 algorithms
ppi_result$downstream_genes

# Broader consensus genes based on consensus_min_methods
ppi_result$consensus_hub_genes
```

After PPI hub screening, users can draw a publication-style UpSet plot for any
subset of algorithms and any top-n cutoff:

```r
upset <- auto_ppi_upset_plot(
  ppi_result,
  methods = c("Degree", "MCC", "MNC", "DMNC", "EPC"),
  top_n = 20,
  output_dir = "STRING_PPI_analysis_results"
)

upset$intersection_summary
```

The first run downloads three STRING resource files:

- `9606.protein.info.v12.0.txt.gz`
- `9606.protein.aliases.v12.0.txt.gz`
- `9606.protein.links.v12.0.txt.gz`

You can also download them explicitly before analysis:

```r
download_string_resources(
  resource_dir = "resource",
  species_id = "9606",
  string_version = "v12.0"
)
```

If downloading inside R is slow or blocked, download the three STRING files by
browser, server command line, or a shared local database, then pass their exact
paths directly. In this manual-file mode, all three paths must be supplied and
AutoML4R will not try to download anything:

```r
ppi_result <- auto_string_ppi_analysis(
  genes = genes,
  info_file = "/path/to/9606.protein.info.v12.0.txt.gz",
  aliases_file = "/path/to/9606.protein.aliases.v12.0.txt.gz",
  links_file = "/path/to/9606.protein.links.v12.0.txt.gz",
  download_resources = FALSE,
  output_dir = "STRING_PPI_analysis_results",
  score_threshold = 400,
  top_n = 10
)
```

For other species, replace `species_id` with the STRING taxonomic identifier
and use matching STRING resource files.

The one-click STRING PPI workflow writes:

- `STRING_offline_mapped_genes.csv`
- `STRING_offline_unmapped_genes.csv`
- `STRING_offline_raw_query_interactions.csv`
- `STRING_offline_edge_df_for_cytohubba.csv`
- `STRING_PPI_top<n>_intersection_all_11_methods.csv`
- `STRING_PPI_consensus_hub_genes.csv`
- all downstream 11-method PPI tables and figures from `auto_ppi_analysis()`,
  including `PPI_top<n>_intersection_all_11_methods.csv`

The lower-level entry point is `auto_ppi_analysis()`. Use it only when you
already have a PPI edge table:

```r
hub_result <- auto_ppi_analysis(
  edge_df = ppi_edges,
  from_col = "source",
  to_col = "target",
  output_dir = "PPI_screening_results",
  top_n = 10
)

# Strict intersection of top-10 nodes across all 11 PPI algorithms
hub_result$downstream_genes
```

`ppi_edges` should contain one PPI edge per row. The network is treated as an
undirected simple graph; duplicate edges and self-loops are removed.

The implemented PPI methods are:

```r
ppi_methods()
```

The PPI workflow writes:

- `PPI_cytohubba_11_scores.csv`
- per-method top-node tables
- integrated rank and consensus hub tables
- algorithm intersection tables from 2 to 11 methods
- PPI figures under `figures/`, including hub-frequency barplot,
  normalized 11-method score heatmap, and consensus-hub network plot

The 11 PPI topology algorithms are Degree, MNC, DMNC, MCC, BottleNeck,
EcCentricity, Closeness, Radiality, Betweenness, Stress, and EPC. The
`auto_string_ppi_analysis()` object stores the full `auto_ppi_analysis()` result
under `ppi_result$ppi`, so advanced users can inspect per-method rankings,
intersection tables, integrated rank, and figures.

Default methods currently include 25 feature-selection algorithms or algorithm variants:

```r
aml_methods()
```

For a runnable simulated-data test:

```sh
Rscript inst/examples/smoke_test.R
```

For a full all-method plot test:

```sh
Rscript inst/examples/full_methods_test.R
```

## GitHub Upload

After replacing the placeholder author email and GitHub username, initialize and push:

```sh
git init
git add .
git commit -m "Create AutoML4R package"
git branch -M main
git remote add origin https://github.com/YOUR_GITHUB_USERNAME/AutoML4R.git
git push -u origin main
```
