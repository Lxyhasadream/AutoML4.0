\name{ppi_methods}
\alias{ppi_methods}
\title{List Implemented PPI Hub-Screening Methods}
\usage{
ppi_methods()
}
\description{
Returns the names of the 11 cytoHubba-style topology algorithms used by
\code{\link{auto_ppi_analysis}} and \code{\link{auto_string_ppi_analysis}}.
}
\value{
A character vector containing Degree, MNC, DMNC, MCC, BottleNeck,
EcCentricity, Closeness, Radiality, Betweenness, Stress, and EPC.
}
\examples{
ppi_methods()
}
