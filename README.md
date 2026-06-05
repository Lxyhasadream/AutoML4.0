# AutoML4R

AutoML4R packages the `auto_ml_analysis()` workflow from `自动机器学习4.0.R` as an installable R package.

## Installation

```r
# install.packages("remotes")
remotes::install_github("YOUR_GITHUB_USERNAME/AutoML4R")
```

Some modelling backends are optional runtime dependencies and may require special installation steps, especially `catboost` and `lightgbm`.

## Usage

```r
library(AutoML4R)

result <- auto_ml_analysis(
  hub_data = hub_data,
  group = group,
  output_dir = "ML_screening_results"
)
```

`hub_data` should be a data frame or matrix with samples in rows and features in columns. `group` should contain one class label per sample.

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
