\name{automl4r_beautiful_plots}
\alias{automl4r_beautiful_plots}
\title{Generate publication-style plots for AutoML4R feature-selection results}
\usage{
automl4r_beautiful_plots(
  result,
  hub_data = NULL,
  group = NULL,
  output_dir = NULL,
  all_genes = NULL,
  top_n_vote = 30,
  heatmap_top_n = 40,
  top_n_per_method = 20,
  expression_top_n_per_method = 12,
  min_support_fraction = 0.3,
  positive_class = NULL
)
}
\arguments{
\item{result}{An \code{auto_ml_analysis} result object returned by \code{auto_ml_analysis()}.}
\item{hub_data}{Optional original feature matrix or data frame. Required for selected-gene expression violin/boxplots.}
\item{group}{Optional class labels matching rows of \code{hub_data}. Required for selected-gene expression violin/boxplots.}
\item{output_dir}{Output directory. If omitted, \code{result$output_dir} is used.}
\item{all_genes}{Optional named list of selected genes by algorithm. If omitted, \code{result$selected_genes} is used.}
\item{top_n_vote}{Number of top consensus genes shown in the vote lollipop plot.}
\item{heatmap_top_n}{Number of top consensus genes shown in the method-by-gene heatmap.}
\item{top_n_per_method}{Number of selected genes shown in per-algorithm score and ranking plots.}
\item{expression_top_n_per_method}{Number of selected genes shown in per-algorithm expression plots.}
\item{min_support_fraction}{Minimum fraction of algorithms used to define consensus genes. The threshold is at least 3 methods.}
\item{positive_class}{Optional positive class label, used to emphasize the disease/case group in expression plots.}
}
\description{
Generates publication-oriented visualizations from an \code{auto_ml_analysis()} result. The function is intended to be run after feature selection and writes Nature-style consensus figures, per-algorithm selected-gene plots, selected-gene expression plots, and model-specific diagnostic plots where model objects support them.
}
\details{
The function writes a final consensus overview figure, consensus vote lollipop plot, method-by-gene consensus heatmap grouped by algorithm class, algorithm output-size plot, and per-algorithm plots under \code{NatureStyle_algorithm_selected_gene_plots/}.

For LASSO, Ridge, and ElasticNet, coefficient-path and cross-validation plots are generated from the stored glmnet models. For random forest, feature-importance and OOB error plots are generated when the stored model is available. Other methods receive model-specific selected-gene ranking or importance plots based on their returned rankings.
}
\value{
Invisibly returns a list containing voting tables, method annotation tables, consensus genes, and ggplot objects for the main consensus figures.
}
\examples{
\dontrun{
result <- auto_ml_analysis(
  hub_data = hub_data,
  group = group,
  positive_class = "Disease",
  output_dir = "ML_screening_results",
  methods = aml_methods()
)

automl4r_beautiful_plots(
  result,
  hub_data = hub_data,
  group = group,
  positive_class = "Disease"
)
}
}
