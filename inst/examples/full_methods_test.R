if (file.exists("R/auto_ml_analysis.R") && file.exists("R/plots.R")) {
  source("R/auto_ml_analysis.R")
  source("R/plots.R")
} else {
  library(AutoML4R)
}

set.seed(2026)
n <- 40
p <- 12
group <- factor(rep(c("Control", "Disease"), each = n / 2))
hub_data <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(hub_data) <- paste0("Gene", seq_len(p))
hub_data[group == "Disease", 1:4] <- hub_data[group == "Disease", 1:4] + 1.3
hub_data <- as.data.frame(hub_data)

result <- auto_ml_analysis(
  hub_data = hub_data,
  group = group,
  positive_class = "Disease",
  output_dir = "ML_screening_results_full_methods_test",
  methods = aml_methods(),
  top_n = 5,
  seed = 2026,
  verbose = TRUE
)

print(names(result$selected_genes))
print(result$failures)
print(head(result$consensus, 10))

stopifnot(length(result$selected_genes) == length(aml_methods()))
stopifnot(length(result$failures) == 0)
