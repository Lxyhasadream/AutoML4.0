\name{download_string_resources}
\alias{download_string_resources}
\title{Download Local STRING Resource Files}
\usage{
download_string_resources(
  resource_dir = "resource",
  species_id = "9606",
  string_version = "v12.0",
  overwrite = FALSE,
  timeout = 600,
  verbose = TRUE
)
}
\arguments{
\item{resource_dir}{Directory where STRING resource files are saved.}
\item{species_id}{STRING species identifier. The default is \code{"9606"} for human.}
\item{string_version}{STRING file version suffix. The default is \code{"v12.0"}.}
\item{overwrite}{Whether to re-download files that already exist locally.}
\item{timeout}{Download timeout in seconds.}
\item{verbose}{Whether to print download progress messages.}
}
\description{
Downloads the three local STRING files required by
\code{\link{auto_string_ppi_analysis}}: protein info, protein aliases, and
protein links. Users can run this once before PPI analysis, or rely on
\code{auto_string_ppi_analysis(download_resources = TRUE)} to download missing
files automatically.
}
\value{
Invisibly returns a named list of local file paths.
}
\examples{
\dontrun{
download_string_resources(resource_dir = "resource", species_id = "9606")
}
}
