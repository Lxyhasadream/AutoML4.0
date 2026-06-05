\name{auto_ml_analysis}
\alias{auto_ml_analysis}
\title{Run automated machine learning feature selection}
\usage{
auto_ml_analysis(
  hub_data,
  group,
  output_dir = "ML_screening_results",
  positive_class = NULL,
  methods = c("univariate", "roc_auc", "lasso", "ridge", "elastic_net",
    "random_forest", "boruta", "xgboost"),
  top_n = 20,
  seed = 123,
  write_outputs = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{hub_data}{A data frame or matrix of candidate features, with samples in rows and features in columns.}
\item{group}{A vector or factor containing one class label per sample.}
\item{output_dir}{Directory where result tables and figures will be written.}
\item{positive_class}{Optional positive class label. If omitted, the second factor level is used.}
\item{methods}{Feature-selection methods to run.}
\item{top_n}{Maximum number of top features selected per method.}
\item{seed}{Random seed for reproducible model fitting.}
\item{write_outputs}{Whether to write CSV, TXT, PDF, and PNG outputs.}
\item{verbose}{Whether to print progress messages.}
}
\description{
Runs a multi-algorithm binary-classification feature-selection workflow and writes selected genes, method rankings, consensus votes, and figures to disk.
}
\value{
An \code{auto_ml_analysis} object containing selected genes by method, per-method rankings, consensus votes, failed methods, the positive class, and the output directory.
}
\examples{
\dontrun{
auto_ml_analysis(hub_data, group, output_dir = "ML_screening_results", positive_class = "Disease")
}
}
