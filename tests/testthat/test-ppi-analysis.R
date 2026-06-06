test_that("auto_ppi_analysis ranks a small PPI network with 11 methods", {
  edge_df <- data.frame(
    source = c("A", "A", "A", "B", "B", "C", "D", "E"),
    target = c("B", "C", "D", "C", "E", "F", "E", "F"),
    stringsAsFactors = FALSE
  )

  result <- auto_ppi_analysis(
    edge_df,
    from_col = "source",
    to_col = "target",
    top_n = 3,
    epc_n_sim = 20,
    seed = 2026,
    write_outputs = FALSE,
    plot_outputs = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "auto_ppi_analysis")
  expect_length(ppi_methods(), 11)
  expect_true(all(ppi_methods() %in% colnames(result$scores)))
  expect_equal(names(result$top_tables), ppi_methods())
  expect_true(nrow(result$hub_frequency) > 0)
  expect_true(nrow(result$method_formula_check) == 11)
})
