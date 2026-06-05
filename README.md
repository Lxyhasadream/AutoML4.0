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
  methods = c("univariate", "roc_auc", "lasso", "ridge", "elastic_net", "random_forest", "boruta", "xgboost"),
  top_n = 20
)
```

`hub_data` should be a data frame or matrix with samples in rows and features in columns. `group` should contain one binary class label per sample.

The workflow writes:

- `feature_rankings_by_method.csv`
- `consensus_feature_votes.csv`
- one selected-gene TXT file per method
- consensus vote and method-support figures

For a runnable simulated-data test:

```sh
Rscript inst/examples/smoke_test.R
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
