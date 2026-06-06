\name{auto_ppi_analysis}
\alias{auto_ppi_analysis}
\title{Run cytoHubba-style PPI hub screening}
\usage{
auto_ppi_analysis(
  edge_df,
  from_col = 1,
  to_col = 2,
  output_dir = "PPI_screening_results",
  top_n = 10,
  consensus_min_methods = NULL,
  epsilon = 1.7,
  epc_threshold = 0.1,
  epc_n_sim = 1000,
  seed = 123,
  include_self_radiality = TRUE,
  write_outputs = TRUE,
  plot_outputs = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{edge_df}{A data frame containing PPI network edges.}
\item{from_col}{Source-node column name or index.}
\item{to_col}{Target-node column name or index.}
\item{output_dir}{Directory where PPI tables and figures are written.}
\item{top_n}{Number of top nodes retained for each PPI algorithm.}
\item{consensus_min_methods}{Minimum number of top lists in which a node must appear to be called a consensus hub. If omitted, at least 30 percent of the 11 methods, with a minimum of 3, is used.}
\item{epsilon}{Exponent used by DMNC. The cytoHubba paper uses 1.7.}
\item{epc_threshold}{Edge percolation threshold. Edges with random values below this threshold are removed in EPC simulations.}
\item{epc_n_sim}{Number of EPC reduced-network simulations.}
\item{seed}{Random seed for EPC simulations and network layout.}
\item{include_self_radiality}{Whether to include the node itself in the radiality summation.}
\item{write_outputs}{Whether to write CSV result tables.}
\item{plot_outputs}{Whether to write PDF and PNG plots.}
\item{verbose}{Whether to print network summary messages.}
}
\description{
Ranks PPI network nodes using 11 cytoHubba-style topology algorithms and writes hub-screening tables, consensus hub summaries, intersection tables, and publication-oriented PPI figures.
}
\details{
The implemented algorithms are Degree, MNC, DMNC, MCC, BottleNeck, EcCentricity, Closeness, Radiality, Betweenness, Stress, and EPC, following the formulas described in the cytoHubba paper. The network is treated as an undirected simple graph; duplicate edges and self-loops are removed.

Main output files include \code{PPI_cytohubba_11_scores.csv}, per-method top-node tables, integrated rank tables, consensus hub tables, algorithm intersection tables, and figures under \code{figures/}.
}
\value{
An \code{auto_ppi_analysis} object containing the igraph graph, score table, top tables, top-node lists, integrated rank table, consensus hubs, intersection results, formula-check table, and output directory.
}
\examples{
\dontrun{
ppi_edges <- data.frame(
  source = c("A", "A", "B", "C"),
  target = c("B", "C", "C", "D")
)

ppi_result <- auto_ppi_analysis(
  ppi_edges,
  from_col = "source",
  to_col = "target",
  output_dir = "PPI_screening_results",
  top_n = 10
)
}
}
