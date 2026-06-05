if (requireNamespace("AutoML4R", quietly = TRUE)) {
  library(AutoML4R)
} else {
  source("R/auto_ml_analysis.R")
}

set.seed(2026)
n <- 60
p <- 15
group <- factor(rep(c("Control", "Disease"), each = n / 2))
hub_data <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(hub_data) <- paste0("Gene", seq_len(p))
hub_data[group == "Disease", 1:4] <- hub_data[group == "Disease", 1:4] + 1.5
hub_data <- as.data.frame(hub_data)

result <- auto_ml_analysis(
  hub_data = hub_data,
  group = group,
  positive_class = "Disease",
  output_dir = "ML_screening_results_smoke_test",
  methods = c("univariate", "roc_auc", "lasso", "ridge", "elastic_net", "random_forest", "boruta", "xgboost"),
  top_n = 8,
  seed = 2026
)

print(head(result$consensus, 10))
