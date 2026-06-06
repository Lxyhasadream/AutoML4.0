library(AutoML4R)

# Replace these objects with your real data before running.
# hub_data <- read.table("hub_data.txt", header = TRUE, row.names = 1, sep = "\t")
# group <- read.table("group.txt", header = TRUE, sep = "\t")$group
#
# result <- auto_ml_analysis(
#   hub_data = hub_data,
#   group = group,
#   positive_class = "Disease",
#   output_dir = "ML_screening_results",
#   methods = aml_methods(),
#   top_n = 20
# )
#
# automl4r_beautiful_plots(
#   result,
#   hub_data = hub_data,
#   group = group,
#   positive_class = "Disease"
# )

# Optional PPI hub screening from an edge table.
# ppi_edges <- read.csv("ppi_edges.csv")
# ppi_result <- auto_ppi_analysis(
#   edge_df = ppi_edges,
#   from_col = "source",
#   to_col = "target",
#   output_dir = "PPI_screening_results",
#   top_n = 10
# )
