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

## PPI Hub Screening

AutoML4R also provides cytoHubba-style PPI hub screening with 11 topology
algorithms from the cytoHubba paper:

```r
ppi_result <- auto_ppi_analysis(
  edge_df = ppi_edges,
  from_col = "source",
  to_col = "target",
  output_dir = "PPI_screening_results",
  top_n = 10
)
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
