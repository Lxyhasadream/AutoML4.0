test_that("auto_ml_analysis runs on simulated binary data", {
  set.seed(2026)
  n <- 48
  p <- 12
  group <- factor(rep(c("Control", "Disease"), each = n / 2))
  x <- matrix(stats::rnorm(n * p), nrow = n, ncol = p)
  colnames(x) <- paste0("Gene", seq_len(p))
  x[group == "Disease", 1:3] <- x[group == "Disease", 1:3] + 1.4
  x <- as.data.frame(x)

  out_dir <- file.path(tempdir(), "automl4r-smoke-test")
  result <- auto_ml_analysis(
    x,
    group,
    output_dir = out_dir,
    positive_class = "Disease",
    methods = c("univariate", "roc_auc", "lasso", "ridge", "random_forest"),
    top_n = 5,
    seed = 2026,
    write_outputs = TRUE,
    verbose = FALSE
  )

  expect_s3_class(result, "auto_ml_analysis")
  expect_true(nrow(result$consensus) > 0)
  expect_true(any(result$consensus$feature %in% c("Gene1", "Gene2", "Gene3")))
  expect_true(file.exists(file.path(out_dir, "consensus_feature_votes.csv")))
  expect_true(file.exists(file.path(out_dir, "feature_rankings_by_method.csv")))
})
