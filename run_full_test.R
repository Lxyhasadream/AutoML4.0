# AutoML4R full local test script
# Run from the project root:
#   source("run_full_test.R")
# or:
#   Rscript run_full_test.R

if (file.exists("R/auto_ml_analysis.R") && file.exists("R/plots.R")) {
  source("R/auto_ml_analysis.R")
  source("R/plots.R")
} else {
  library(AutoML4R)
}

set.seed(20260606)

n_samples <- 40
n_features <- 12

group <- factor(rep(c("Control", "Disease"), each = n_samples / 2))
hub_data <- matrix(rnorm(n_samples * n_features), nrow = n_samples, ncol = n_features)
colnames(hub_data) <- paste0("Gene", seq_len(n_features))

# Add real signal to the first four features so the methods have something to recover.
hub_data[group == "Disease", 1:4] <- hub_data[group == "Disease", 1:4] + 1.3
hub_data <- as.data.frame(hub_data)

output_dir <- "ML_screening_results_user_test"

cat("Available methods:\n")
print(aml_methods())
cat("\nRunning full AutoML test...\n")

result <- auto_ml_analysis(
  hub_data = hub_data,
  group = group,
  positive_class = "Disease",
  output_dir = output_dir,
  methods = aml_methods(),
  top_n = 5,
  seed = 20260606,
  write_outputs = TRUE,
  verbose = TRUE
)

cat("\nCompleted methods:\n")
print(names(result$selected_genes))

cat("\nFailed methods:\n")
print(result$failures)

cat("\nTop consensus features:\n")
print(head(result$consensus, 10))

stopifnot(length(result$selected_genes) == length(aml_methods()))
stopifnot(length(result$failures) == 0)
stopifnot(file.exists(file.path(output_dir, "consensus_feature_votes.csv")))
stopifnot(file.exists(file.path(output_dir, "feature_rankings_by_method.csv")))
stopifnot(file.exists(file.path(output_dir, "consensus_feature_votes.pdf")))
stopifnot(file.exists(file.path(output_dir, "method_support_heatmap.pdf")))

cat("\nTest passed. Output directory:\n")
cat(normalizePath(output_dir, mustWork = FALSE), "\n")
