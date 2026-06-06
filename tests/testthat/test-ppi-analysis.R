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
  expect_equal(result$intersection_genes_all_11_methods, sort(unique(Reduce(intersect, result$top_nodes))))
  expect_equal(result$downstream_genes, result$intersection_genes_all_11_methods)
  expect_true(nrow(result$hub_frequency) > 0)
  expect_true(nrow(result$method_formula_check) == 11)
})

test_that("auto_ppi_analysis writes the strict 11-method intersection file", {
  edge_df <- data.frame(
    source = c("A", "A", "A", "B", "B", "C", "D", "E"),
    target = c("B", "C", "D", "C", "E", "F", "E", "F"),
    stringsAsFactors = FALSE
  )
  out_dir <- file.path(tempdir(), "ppi_strict_intersection_output")

  result <- auto_ppi_analysis(
    edge_df,
    from_col = "source",
    to_col = "target",
    output_dir = out_dir,
    top_n = 3,
    epc_n_sim = 5,
    write_outputs = TRUE,
    plot_outputs = FALSE,
    verbose = FALSE
  )

  expect_equal(result$downstream_genes, result$intersection_genes_all_11_methods)
  expect_true(file.exists(file.path(out_dir, "PPI_top3_intersection_all_11_methods.csv")))
})

test_that("auto_string_ppi_analysis builds local STRING edges and ranks hubs", {
  resource_dir <- file.path(tempdir(), "string_resource")
  out_dir <- file.path(tempdir(), "string_ppi_output")
  dir.create(resource_dir, showWarnings = FALSE, recursive = TRUE)

  write_gz <- function(path, lines) {
    con <- gzfile(path, open = "wt")
    on.exit(close(con), add = TRUE)
    writeLines(lines, con)
  }

  write_gz(
    file.path(resource_dir, "9606.protein.info.v12.0.txt.gz"),
    c(
      "string_protein_id\tpreferred_name",
      "9606.P1\tTP53",
      "9606.P2\tEGFR",
      "9606.P3\tMYC"
    )
  )
  write_gz(
    file.path(resource_dir, "9606.protein.aliases.v12.0.txt.gz"),
    c(
      "string_protein_id\talias\tsource",
      "9606.P3\tALIAS_MYC\tEnsembl"
    )
  )
  write_gz(
    file.path(resource_dir, "9606.protein.links.v12.0.txt.gz"),
    c(
      "protein1 protein2 combined_score",
      "9606.P1 9606.P2 900",
      "9606.P2 9606.P1 900",
      "9606.P1 9606.P3 700",
      "9606.P3 9606.P1 700",
      "9606.P2 9606.P3 300"
    )
  )

  result <- auto_string_ppi_analysis(
    genes = c("TP53", "EGFR", "ALIAS_MYC"),
    resource_dir = file.path(tempdir(), "missing_standard_resource_dir"),
    info_file = file.path(resource_dir, "9606.protein.info.v12.0.txt.gz"),
    aliases_file = file.path(resource_dir, "9606.protein.aliases.v12.0.txt.gz"),
    links_file = file.path(resource_dir, "9606.protein.links.v12.0.txt.gz"),
    output_dir = out_dir,
    download_resources = FALSE,
    score_threshold = 400,
    top_n = 2,
    epc_n_sim = 5,
    write_outputs = TRUE,
    plot_outputs = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "auto_string_ppi_analysis")
  expect_s3_class(result$ppi, "auto_ppi_analysis")
  expect_equal(nrow(result$edge_df), 2)
  expect_true("ALIAS_MYC" %in% result$mapping$gene_input)
  expect_equal(result$downstream_genes, result$intersection_genes_all_11_methods)
  expect_equal(result$downstream_genes, result$ppi$intersection_genes_all_11_methods)
  expect_true(file.exists(file.path(out_dir, "STRING_PPI_top2_intersection_all_11_methods.csv")))
  expect_true(file.exists(file.path(out_dir, "PPI_top2_intersection_all_11_methods.csv")))
  expect_true(file.exists(file.path(out_dir, "STRING_offline_edge_df_for_cytohubba.csv")))
})

test_that("auto_ppi_upset_plot draws selected-method intersections", {
  edge_df <- data.frame(
    source = c("A", "A", "A", "B", "B", "C", "D", "E"),
    target = c("B", "C", "D", "C", "E", "F", "E", "F"),
    stringsAsFactors = FALSE
  )
  ppi <- auto_ppi_analysis(
    edge_df,
    from_col = "source",
    to_col = "target",
    top_n = 3,
    epc_n_sim = 5,
    write_outputs = FALSE,
    plot_outputs = FALSE,
    verbose = FALSE
  )
  out_dir <- file.path(tempdir(), "ppi_upset_output")
  upset <- auto_ppi_upset_plot(
    ppi,
    methods = c("Degree", "MCC", "MNC", "EPC"),
    top_n = 3,
    output_dir = out_dir,
    file_prefix = "test_upset",
    write_outputs = TRUE
  )

  expect_s3_class(upset, "auto_ppi_upset_plot")
  expect_true(nrow(upset$intersection_summary) > 0)
  expect_equal(names(upset$selected_sets), c("Degree", "MCC", "MNC", "EPC"))
  expect_true(file.exists(file.path(out_dir, "test_upset_intersection_summary.csv")))
  expect_true(file.exists(file.path(out_dir, "test_upset.png")))
})
