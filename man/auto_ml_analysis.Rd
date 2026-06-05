\name{auto_ml_analysis}
\alias{auto_ml_analysis}
\title{Run automated machine learning feature selection}
\usage{
auto_ml_analysis(hub_data, group, output_dir = "ML_screening_results")
}
\arguments{
\item{hub_data}{A data frame or matrix of candidate features, with samples in rows and features in columns.}
\item{group}{A vector or factor containing one class label per sample.}
\item{output_dir}{Directory where result tables and figures will be written.}
}
\description{
Runs a multi-algorithm automated machine learning feature-selection workflow and writes selected genes, model summaries, and consensus figures to disk.
}
\value{
The function is primarily called for its output files and console progress messages. It invisibly returns \code{NULL} unless the workflow code returns another value.
}
\examples{
\dontrun{
auto_ml_analysis(hub_data, group, output_dir = "ML_screening_results")
}
}
