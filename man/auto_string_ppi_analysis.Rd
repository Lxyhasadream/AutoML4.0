\name{auto_string_ppi_analysis}
\alias{auto_string_ppi_analysis}
\title{Offline STRING PPI Construction and Hub Screening}
\usage{
auto_string_ppi_analysis(
  genes,
  resource_dir = "resource",
  output_dir = "STRING_PPI_analysis_results",
  species_id = "9606",
  string_version = "v12.0",
  aliases_file = NULL,
  info_file = NULL,
  links_file = NULL,
  download_resources = TRUE,
  overwrite_resources = FALSE,
  score_threshold = 400,
  top_n = 10,
  consensus_min_methods = NULL,
  alias_chunk_size = 300000,
  links_chunk_size = 500000,
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
\item{genes}{Character vector of gene symbols or aliases to query.}
\item{resource_dir}{Directory containing local STRING resource files.}
\item{output_dir}{Directory where mapping, edge, score, and hub result files are written.}
\item{species_id}{STRING species identifier. The default is \code{"9606"} for human.}
\item{string_version}{STRING file version suffix. The default is \code{"v12.0"}.}
\item{aliases_file}{Optional full path to a local STRING protein aliases \code{.txt.gz} file.}
\item{info_file}{Optional full path to a local STRING protein info \code{.txt.gz} file.}
\item{links_file}{Optional full path to a local STRING protein links \code{.txt.gz} file.}
\item{download_resources}{Whether to download missing STRING resource files automatically before analysis.}
\item{overwrite_resources}{Whether to re-download and overwrite existing local STRING resource files.}
\item{score_threshold}{Minimum STRING \code{combined_score} to keep. Typical values are 150, 400, 700, and 900.}
\item{top_n}{Number of top genes retained for each topology method.}
\item{consensus_min_methods}{Minimum number of top lists a gene must appear in to be reported as a consensus hub.}
\item{alias_chunk_size}{Number of alias-file lines read per streaming chunk.}
\item{links_chunk_size}{Number of links-file lines read per streaming chunk.}
\item{epsilon}{Exponent used by the DMNC topology score.}
\item{epc_threshold}{Edge removal probability used by the EPC score.}
\item{epc_n_sim}{Number of EPC edge-percolation simulations.}
\item{seed}{Random seed used by EPC.}
\item{include_self_radiality}{Whether radiality includes the node itself.}
\item{write_outputs}{Whether CSV outputs are written.}
\item{plot_outputs}{Whether PPI figures are written.}
\item{verbose}{Whether progress messages are printed.}
}
\description{
Maps input genes to local STRING protein identifiers, extracts query-query
interactions from local STRING link files, collapses protein-level interactions
to gene-level edges, and runs \code{\link{auto_ppi_analysis}} on the resulting
network. If the required STRING files are absent and \code{download_resources}
is \code{TRUE}, the function downloads them into \code{resource_dir} first.
The final downstream gene set is the strict intersection of genes appearing in
the top-\code{n} lists from all 11 PPI topology algorithms.

If downloading through R is slow or blocked, users may manually download the
three STRING files and pass their exact locations through \code{aliases_file},
\code{info_file}, and \code{links_file}. When any manual file path is supplied,
all three file paths must be supplied, and \code{resource_dir} lookup/download is
skipped.
}
\details{
For a new user, the intended workflow is:
\enumerate{
\item provide a character vector of gene symbols;
\item let AutoML4R download the required STRING files into \code{resource_dir};
\item map gene symbols and aliases to STRING protein identifiers;
\item extract STRING query-query interactions above \code{score_threshold};
\item collapse protein-level links to gene-level edges;
\item run 11 cytoHubba-style algorithms;
\item return and write the strict 11-method intersection gene set.
}

The three required local STRING files are
\code{<species_id>.protein.info.<string_version>.txt.gz},
\code{<species_id>.protein.aliases.<string_version>.txt.gz}, and
\code{<species_id>.protein.links.<string_version>.txt.gz}.
}
\value{
An object of class \code{auto_string_ppi_analysis}, a list containing input
genes, mapping results, unmapped genes, raw STRING interactions, gene-level edge
table, the downstream \code{auto_ppi_analysis} result, consensus hub genes, and
\code{downstream_genes}, which is the strict top-\code{n} intersection across all
11 methods.
}
\examples{
\dontrun{
genes <- c("SERPINE1", "PF4", "HMOX1", "PPBP", "UCP2")
res <- auto_string_ppi_analysis(
  genes,
  resource_dir = "resource",
  output_dir = "STRING_PPI_analysis_results",
  score_threshold = 400,
  top_n = 10
)
res$edge_df
res$ppi$hub_frequency
res$downstream_genes

# Manual local-file mode for users who downloaded STRING files separately.
res2 <- auto_string_ppi_analysis(
  genes,
  info_file = "/path/to/9606.protein.info.v12.0.txt.gz",
  aliases_file = "/path/to/9606.protein.aliases.v12.0.txt.gz",
  links_file = "/path/to/9606.protein.links.v12.0.txt.gz",
  download_resources = FALSE
)
}
}
