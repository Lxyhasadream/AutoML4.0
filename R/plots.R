aml_plot_algorithm_panels <- function(results,
                                      output_dir,
                                      top_n = 20,
                                      x = NULL,
                                      y = NULL,
                                      positive_class = NULL) {
  plot_dir <- file.path(output_dir, "algorithm_plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  for (method in names(results)) {
    res <- results[[method]]
    method_dir <- file.path(plot_dir, method)
    dir.create(method_dir, showWarnings = FALSE, recursive = TRUE)

    try(aml_plot_ranking_bar(res, method_dir, top_n), silent = TRUE)
    try(aml_plot_selected_gene_scores(res, method_dir, top_n), silent = TRUE)
    try(aml_plot_selected_gene_expression(res, method_dir, x, y, top_n, positive_class), silent = TRUE)
    try(aml_plot_special_method(method, res, method_dir, top_n), silent = TRUE)
  }
}

aml_method_label <- function(method) {
  words <- strsplit(gsub("_", " ", method), " ", fixed = TRUE)[[1]]
  paste(toupper(substr(words, 1, 1)), substr(words, 2, nchar(words)), sep = "", collapse = " ")
}

aml_plot_palette <- function() {
  c("#2A6F97", "#F28E2B", "#4E79A7", "#59A14F", "#E15759", "#B07AA1", "#76B7B2")
}

aml_plot_typography <- function() {
  list(
    family = "Arial",
    axis_text = 18,
    axis_title = 20,
    title = 24
  )
}

aml_apply_plot_typography <- function(plot) {
  typo <- aml_plot_typography()
  plot +
    ggplot2::theme(
      text = ggplot2::element_text(family = typo$family),
      plot.title = ggplot2::element_text(family = typo$family, size = typo$title, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = ggplot2::element_text(family = typo$family, size = typo$title, color = "black"),
      axis.title = ggplot2::element_text(family = typo$family, size = typo$axis_title, face = "bold", color = "black"),
      axis.title.x = ggplot2::element_text(family = typo$family, size = typo$axis_title, face = "bold", color = "black"),
      axis.title.y = ggplot2::element_text(family = typo$family, size = typo$axis_title, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(family = typo$family, size = typo$axis_text, color = "black"),
      strip.text = ggplot2::element_text(family = typo$family, size = typo$title, face = "bold", color = "black"),
      legend.title = ggplot2::element_text(family = typo$family, size = typo$axis_text, face = "bold", color = "black"),
      legend.text = ggplot2::element_text(family = typo$family, size = typo$axis_text, color = "black")
    )
}

aml_base_plot_par <- function() {
  graphics::par(
    family = "Arial",
    cex.axis = 1.5,
    cex.lab = 1.67,
    cex.main = 2.0,
    font.lab = 2,
    font.main = 2
  )
}

aml_publication_theme <- function(base_size = 13) {
  aml_need_package("ggplot2")
  typo <- aml_plot_typography()
  ggplot2::theme_minimal(base_size = base_size, base_family = typo$family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = typo$family),
      plot.title = ggplot2::element_text(family = typo$family, size = typo$title, face = "bold", color = "#1F2933", margin = ggplot2::margin(b = 6)),
      plot.subtitle = ggplot2::element_text(family = typo$family, size = typo$title, color = "#52616B", margin = ggplot2::margin(b = 12)),
      axis.title = ggplot2::element_text(family = typo$family, size = typo$axis_title, color = "#1F2933"),
      axis.text = ggplot2::element_text(family = typo$family, size = typo$axis_text, color = "#334E68"),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "top",
      legend.title = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(family = typo$family, size = typo$title, face = "bold", color = "#1F2933"),
      strip.background = ggplot2::element_rect(fill = "#F4F7F9", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

aml_plot_ranking_bar <- function(res, method_dir, top_n) {
  aml_need_package("ggplot2")
  df <- res$rankings
  df <- df[order(df$rank), , drop = FALSE]
  df <- head(df, min(top_n, nrow(df)))
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }
  df$feature <- factor(df$feature, levels = rev(df$feature))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = feature)) +
    ggplot2::geom_col(fill = "#2A6F97", width = 0.68) +
    ggplot2::geom_text(ggplot2::aes(label = rank), x = 0, hjust = -0.35, color = "white", size = 3.2) +
    aml_publication_theme(base_size = 13) +
    ggplot2::labs(
      x = "Feature score",
      y = NULL,
      title = paste0(aml_method_label(unique(df$method)[1]), " feature ranking"),
      subtitle = paste0("Top ", nrow(df), " ranked features")
    ) +
    ggplot2::coord_cartesian(clip = "off")
  aml_save_plot(p, file.path(method_dir, "feature_ranking_barplot"), 8, max(5, 0.28 * nrow(df) + 2))
}

aml_plot_selected_gene_scores <- function(res, method_dir, top_n) {
  aml_need_package("ggplot2")
  selected <- head(unique(as.character(res$selected)), top_n)
  df <- res$rankings[res$rankings$feature %in% selected, , drop = FALSE]
  df <- df[order(df$rank), , drop = FALSE]
  if (nrow(df) == 0) {
    return(invisible(NULL))
  }
  df$feature <- factor(df$feature, levels = rev(df$feature))
  df$rank_label <- paste0("#", df$rank)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = feature)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = score, yend = feature),
      linewidth = 1.1,
      color = "#D5E4EE",
      lineend = "round"
    ) +
    ggplot2::geom_point(size = 4.2, color = "#F28E2B") +
    ggplot2::geom_text(ggplot2::aes(label = rank_label), hjust = -0.28, size = 3.3, color = "#1F2933") +
    aml_publication_theme(base_size = 13) +
    ggplot2::labs(
      x = "Selection score",
      y = NULL,
      title = paste0(aml_method_label(unique(df$method)[1]), " selected genes"),
      subtitle = "Selected genes ranked by algorithm-specific score"
    ) +
    ggplot2::coord_cartesian(clip = "off")
  aml_save_plot(p, file.path(method_dir, "selected_gene_scores"), 8.5, max(5, 0.34 * nrow(df) + 2.2))
}

aml_plot_selected_gene_expression <- function(res, method_dir, x, y, top_n, positive_class = NULL) {
  aml_need_package("ggplot2")
  if (is.null(x) || is.null(y)) {
    return(invisible(NULL))
  }
  x <- as.data.frame(x, check.names = FALSE)
  selected <- head(intersect(unique(as.character(res$selected)), colnames(x)), min(top_n, 12))
  if (length(selected) == 0) {
    return(invisible(NULL))
  }
  y <- droplevels(as.factor(y))
  feature_levels <- rev(selected)
  plot_df <- do.call(rbind, lapply(selected, function(feature) {
    data.frame(
      feature = feature,
      group = y,
      expression = x[[feature]],
      stringsAsFactors = FALSE
    )
  }))
  plot_df$feature <- factor(plot_df$feature, levels = feature_levels)
  plot_df$group <- factor(plot_df$group, levels = levels(y))

  fill_values <- setNames(rep(aml_plot_palette(), length.out = nlevels(y)), levels(y))
  if (!is.null(positive_class) && positive_class %in% names(fill_values)) {
    fill_values[positive_class] <- "#F28E2B"
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = expression, fill = group, color = group)) +
    ggplot2::geom_violin(width = 0.78, alpha = 0.18, linewidth = 0.25, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.22, alpha = 0.72, outlier.shape = NA, linewidth = 0.35) +
    ggplot2::geom_jitter(width = 0.12, alpha = 0.58, size = 1.35, stroke = 0) +
    ggplot2::facet_wrap(stats::as.formula("~ feature"), scales = "free_y", ncol = 4) +
    ggplot2::scale_fill_manual(values = fill_values) +
    ggplot2::scale_color_manual(values = fill_values) +
    aml_publication_theme(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(family = "Arial", size = 18, angle = 35, hjust = 1),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      x = NULL,
      y = "Expression",
      title = paste0(aml_method_label(unique(res$rankings$method)[1]), " selected-gene expression"),
      subtitle = paste0("Top ", length(selected), " selected genes across phenotype groups")
    )

  height <- max(5.5, 2.2 * ceiling(length(selected) / 4))
  aml_save_plot(p, file.path(method_dir, "selected_gene_expression"), 11, height)
}

aml_plot_special_method <- function(method, res, method_dir, top_n) {
  if (method %in% c("lasso", "ridge", "elastic_net")) {
    aml_plot_glmnet_process(res, method_dir)
  } else if (method == "random_forest") {
    aml_plot_random_forest_process(res, method_dir)
  } else if (method == "ranger") {
    aml_plot_ranger_process(res, method_dir)
  } else if (method == "boruta") {
    aml_plot_boruta_process(res, method_dir)
  } else if (method %in% c("xgboost", "lightgbm", "catboost")) {
    aml_plot_boosting_process(method, res, method_dir, top_n)
  } else if (method == "svm") {
    aml_plot_caret_train_process(res, method_dir, "SVM tuning performance")
  } else if (method == "svm_rfe") {
    aml_plot_svm_rfe_process(res, method_dir)
  } else if (method == "splsda") {
    aml_plot_splsda_process(res, method_dir, top_n)
  } else if (method == "pam") {
    aml_plot_pam_process(res, method_dir)
  } else if (method == "stability_selection") {
    aml_plot_stability_process(res, method_dir, top_n)
  } else if (method %in% c("genetic_algorithm", "simulated_annealing")) {
    aml_plot_search_process(method, res, method_dir, top_n)
  }
}

aml_plot_glmnet_process <- function(res, method_dir) {
  cv <- res$details$cv
  if (is.null(cv)) {
    return(invisible(NULL))
  }
  grDevices::cairo_pdf(file.path(method_dir, "cv_curve.pdf"), width = 7, height = 6, family = "Arial")
  aml_base_plot_par()
  plot(cv)
  graphics::title("Cross-validation curve")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "cv_curve.png"), width = 1800, height = 1500, res = 240)
  aml_base_plot_par()
  plot(cv)
  graphics::title("Cross-validation curve")
  grDevices::dev.off()

  grDevices::cairo_pdf(file.path(method_dir, "coefficient_path.pdf"), width = 8, height = 6, family = "Arial")
  aml_base_plot_par()
  plot(cv$glmnet.fit, xvar = "lambda", label = FALSE)
  graphics::abline(v = log(cv$lambda.min), lty = 2, col = "#D95F02", lwd = 2)
  graphics::abline(v = log(cv$lambda.1se), lty = 3, col = "#1B9E77", lwd = 2)
  graphics::legend("topright", c("lambda.min", "lambda.1se"), lty = c(2, 3), col = c("#D95F02", "#1B9E77"), bty = "n")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "coefficient_path.png"), width = 2000, height = 1500, res = 240)
  aml_base_plot_par()
  plot(cv$glmnet.fit, xvar = "lambda", label = FALSE)
  graphics::abline(v = log(cv$lambda.min), lty = 2, col = "#D95F02", lwd = 2)
  graphics::abline(v = log(cv$lambda.1se), lty = 3, col = "#1B9E77", lwd = 2)
  graphics::legend("topright", c("lambda.min", "lambda.1se"), lty = c(2, 3), col = c("#D95F02", "#1B9E77"), bty = "n")
  grDevices::dev.off()
}

aml_plot_random_forest_process <- function(res, method_dir) {
  fit <- res$details$model
  if (is.null(fit)) {
    return(invisible(NULL))
  }
  grDevices::cairo_pdf(file.path(method_dir, "oob_error_curve.pdf"), width = 7, height = 6, family = "Arial")
  aml_base_plot_par()
  plot(fit, main = "Random forest OOB error")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "oob_error_curve.png"), width = 1800, height = 1500, res = 240)
  aml_base_plot_par()
  plot(fit, main = "Random forest OOB error")
  grDevices::dev.off()

  grDevices::cairo_pdf(file.path(method_dir, "variable_importance_dotchart.pdf"), width = 8, height = 7, family = "Arial")
  aml_base_plot_par()
  randomForest::varImpPlot(fit, main = "Random forest variable importance")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "variable_importance_dotchart.png"), width = 1800, height = 1600, res = 240)
  aml_base_plot_par()
  randomForest::varImpPlot(fit, main = "Random forest variable importance")
  grDevices::dev.off()
}

aml_plot_ranger_process <- function(res, method_dir) {
  imp <- res$details$model$variable.importance
  if (is.null(imp)) {
    return(invisible(NULL))
  }
  aml_plot_named_score(imp, "Ranger permutation importance", method_dir, "permutation_importance")
}

aml_plot_boruta_process <- function(res, method_dir) {
  fixed <- res$details$fixed_model
  if (is.null(fixed)) {
    return(invisible(NULL))
  }
  grDevices::cairo_pdf(file.path(method_dir, "boruta_importance_boxplot.pdf"), width = 10, height = 7, family = "Arial")
  aml_base_plot_par()
  plot(fixed, las = 2, cex.axis = 0.8, main = "Boruta feature importance")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "boruta_importance_boxplot.png"), width = 2400, height = 1800, res = 240)
  aml_base_plot_par()
  plot(fixed, las = 2, cex.axis = 0.8, main = "Boruta feature importance")
  grDevices::dev.off()

  grDevices::cairo_pdf(file.path(method_dir, "boruta_importance_history.pdf"), width = 10, height = 7, family = "Arial")
  aml_base_plot_par()
  Boruta::plotImpHistory(fixed, main = "Boruta importance history")
  grDevices::dev.off()
}

aml_plot_boosting_process <- function(method, res, method_dir, top_n) {
  if (method == "xgboost" && !is.null(res$details$importance)) {
    imp <- res$details$importance
    score <- imp$Gain
    names(score) <- imp$Feature
    aml_plot_named_score(score, "XGBoost gain importance", method_dir, "gain_importance")
  }

  if (method == "lightgbm") {
    if (!is.null(res$details$importance)) {
      imp <- res$details$importance
      score <- imp$Gain
      names(score) <- imp$Feature
      aml_plot_named_score(score, "LightGBM gain importance", method_dir, "gain_importance")
    }
    aml_plot_lightgbm_shap(res, method_dir, top_n)
  }

  if (method == "catboost") {
    aml_plot_ranking_bar(res, method_dir, top_n)
  }
}

aml_plot_lightgbm_shap <- function(res, method_dir, top_n) {
  fit <- res$details$model
  x <- res$details$x
  if (is.null(fit) || is.null(x)) {
    return(invisible(NULL))
  }
  contrib <- tryCatch(
    stats::predict(fit, as.matrix(x), type = "contrib"),
    error = function(e) {
      tryCatch(
        stats::predict(fit, as.matrix(x), predcontrib = TRUE),
        error = function(e2) NULL
      )
    }
  )
  if (is.null(contrib)) {
    return(invisible(NULL))
  }
  if (is.null(colnames(contrib))) {
    if (ncol(contrib) == ncol(x) + 1) {
      colnames(contrib) <- c(colnames(x), "BIAS")
    } else if (ncol(contrib) == ncol(x)) {
      colnames(contrib) <- colnames(x)
    }
  }
  contrib <- as.data.frame(contrib, check.names = FALSE)
  contrib <- contrib[, colnames(contrib) != "BIAS", drop = FALSE]
  shap_score <- vapply(contrib, function(z) mean(abs(z), na.rm = TRUE), numeric(1))
  aml_plot_named_score(shap_score, "LightGBM mean absolute SHAP contribution", method_dir, "shap_mean_abs")

  top_features <- names(sort(shap_score, decreasing = TRUE))[seq_len(min(top_n, length(shap_score)))]
  long <- do.call(rbind, lapply(top_features, function(feature) {
    data.frame(feature = feature, shap_value = contrib[[feature]], expression = x[[feature]], stringsAsFactors = FALSE)
  }))
  long$feature <- factor(long$feature, levels = rev(top_features))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = shap_value, y = feature, color = expression)) +
    ggplot2::geom_jitter(height = 0.18, width = 0, alpha = 0.75, size = 1.8) +
    ggplot2::scale_color_gradient(low = "#2C7FB8", high = "#D95F02") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = "SHAP contribution", y = NULL, color = "Expression", title = "LightGBM SHAP summary")
  aml_save_plot(p, file.path(method_dir, "shap_summary"), 8, max(5, 0.28 * length(top_features) + 2))
}

aml_plot_caret_train_process <- function(res, method_dir, title) {
  fit <- res$details$model
  if (is.null(fit) || is.null(fit$results)) {
    return(invisible(NULL))
  }
  df <- fit$results
  metric <- if ("Accuracy" %in% names(df)) "Accuracy" else names(df)[vapply(df, is.numeric, logical(1))][1]
  xcol <- setdiff(names(df)[vapply(df, is.numeric, logical(1))], c(metric, paste0(metric, "SD")))[1]
  if (is.na(xcol)) {
    df$.index <- seq_len(nrow(df))
    xcol <- ".index"
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[xcol]], y = .data[[metric]])) +
    ggplot2::geom_line(color = "#2C7FB8") +
    ggplot2::geom_point(size = 2.2, color = "#D95F02") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = xcol, y = metric, title = title)
  aml_save_plot(p, file.path(method_dir, "tuning_performance"), 7, 5)
}

aml_plot_svm_rfe_process <- function(res, method_dir) {
  fit <- res$details$model
  if (is.null(fit) || is.null(fit$results)) {
    return(invisible(NULL))
  }
  df <- fit$results
  metric <- if ("Accuracy" %in% names(df)) "Accuracy" else names(df)[vapply(df, is.numeric, logical(1))][1]
  p <- ggplot2::ggplot(df, ggplot2::aes(x = Variables, y = .data[[metric]])) +
    ggplot2::geom_line(color = "#2C7FB8") +
    ggplot2::geom_point(size = 2.2, color = "#D95F02") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = "Number of variables", y = metric, title = "SVM-RFE performance by subset size")
  aml_save_plot(p, file.path(method_dir, "rfe_subset_performance"), 7, 5)
}

aml_plot_splsda_process <- function(res, method_dir, top_n) {
  loadings <- res$details$loadings
  if (is.null(loadings)) {
    return(invisible(NULL))
  }
  score <- abs(loadings[, 1])
  names(score) <- rownames(loadings)
  aml_plot_named_score(score, "sPLS-DA component 1 loadings", method_dir, "component1_loadings", top_n)
}

aml_plot_pam_process <- function(res, method_dir) {
  fit <- res$details$model
  if (is.null(fit) || is.null(fit$threshold)) {
    return(invisible(NULL))
  }
  df <- data.frame(index = seq_along(fit$threshold), threshold = fit$threshold)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = index, y = threshold)) +
    ggplot2::geom_line(color = "#2C7FB8") +
    ggplot2::geom_point(size = 1.8, color = "#D95F02") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = "Threshold index", y = "Threshold", title = "PAM threshold path")
  aml_save_plot(p, file.path(method_dir, "threshold_path"), 7, 5)
}

aml_plot_stability_process <- function(res, method_dir, top_n) {
  freq <- res$details$selection_frequency
  if (is.null(freq)) {
    return(invisible(NULL))
  }
  aml_plot_named_score(freq, "Stability selection frequency", method_dir, "selection_frequency", top_n)
}

aml_plot_search_process <- function(method, res, method_dir, top_n) {
  score <- res$details$base_score
  if (is.null(score)) {
    return(invisible(NULL))
  }
  title <- if (method == "genetic_algorithm") "Genetic-style search base scores" else "Simulated annealing base scores"
  aml_plot_named_score(score, title, method_dir, "search_base_scores", top_n)
}

aml_plot_named_score <- function(score, title, method_dir, file_stem, top_n = 20) {
  aml_need_package("ggplot2")
  score <- score[is.finite(score)]
  score <- sort(score, decreasing = TRUE)
  score <- head(score, min(top_n, length(score)))
  if (length(score) == 0) {
    return(invisible(NULL))
  }
  df <- data.frame(feature = names(score), score = as.numeric(score), stringsAsFactors = FALSE)
  df$feature <- factor(df$feature, levels = rev(df$feature))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = score, y = feature)) +
    ggplot2::geom_col(fill = "#2C7FB8", width = 0.72) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = "Score", y = NULL, title = title)
  aml_save_plot(p, file.path(method_dir, file_stem), 8, max(5, 0.28 * nrow(df) + 2))
}

aml_save_plot <- function(plot, file_base, width, height) {
  plot <- aml_apply_plot_typography(plot)
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
}
