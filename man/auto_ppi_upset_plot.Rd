\name{auto_ppi_upset_plot}
\alias{auto_ppi_upset_plot}
\title{Draw a Publication-Style UpSet Plot for PPI Hub Algorithms}
\usage{
auto_ppi_upset_plot(
  ppi_result,
  methods = ppi_methods(),
  top_n = 10,
  output_dir = NULL,
  file_prefix = NULL,
  max_intersections = 30,
  min_intersection_size = 1,
  width = 11,
  height = 8,
  dpi = 600,
  title = NULL,
  subtitle = NULL,
  write_outputs = TRUE
)
}
\arguments{
\item{ppi_result}{An \code{auto_ppi_analysis} or \code{auto_string_ppi_analysis} result.}
\item{methods}{Character vector of PPI algorithms to include. Values should come from \code{\link{ppi_methods}}.}
\item{top_n}{Number of top-ranked genes selected from each algorithm before intersection analysis.}
\item{output_dir}{Directory where the plot and intersection tables are written. If omitted, the PPI result output directory is used.}
\item{file_prefix}{Prefix for output PDF, PNG, and CSV files.}
\item{max_intersections}{Maximum number of intersections displayed, ordered by intersection size.}
\item{min_intersection_size}{Minimum intersection size retained in the plot.}
\item{width}{Plot width in inches.}
\item{height}{Plot height in inches.}
\item{dpi}{PNG resolution.}
\item{title}{Optional plot title.}
\item{subtitle}{Optional plot subtitle.}
\item{write_outputs}{Whether to write PDF, PNG, and CSV outputs.}
}
\description{
Builds a publication-style UpSet plot from selected PPI hub-screening
algorithms. The function lets users choose which algorithms to include and how
many top genes to take from each algorithm. It then summarizes exact
intersections among those top gene sets and draws a two-panel UpSet figure.
}
\details{
When the input result contains a full PPI score table, top genes are recalculated
from the requested \code{top_n}; therefore the plot is not limited by the
\code{top_n} used when \code{\link{auto_ppi_analysis}} was first run.

The function writes an intersection summary CSV and gene-membership CSV in
addition to PDF and PNG figures.
}
\value{
An object of class \code{auto_ppi_upset_plot} containing the patchwork plot,
selected top-gene sets, intersection summary, gene-membership table, matrix
plot data, and output directory.
}
\examples{
\dontrun{
ppi <- auto_ppi_analysis(ppi_edges, from_col = "source", to_col = "target")
up <- auto_ppi_upset_plot(
  ppi,
  methods = c("Degree", "MCC", "MNC", "EPC"),
  top_n = 20
)
up$intersection_summary
}
}
