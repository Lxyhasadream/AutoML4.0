aml_plot_algorithm_panels <- function(results, output_dir, top_n = 20) {
  plot_dir <- file.path(output_dir, "algorithm_plots")
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

  for (method in names(results)) {
    res <- results[[method]]
    method_dir <- file.path(plot_dir, method)
    dir.create(method_dir, showWarnings = FALSE, recursive = TRUE)

    try(aml_plot_ranking_bar(res, method_dir, top_n), silent = TRUE)
    try(aml_plot_special_method(method, res, method_dir, top_n), silent = TRUE)
  }
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
    ggplot2::geom_col(fill = "#2C7FB8", width = 0.72) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(
      x = "Feature score",
      y = NULL,
      title = paste0(unique(df$method)[1], " feature ranking")
    )
  aml_save_plot(p, file.path(method_dir, "feature_ranking_barplot"), 8, max(5, 0.28 * nrow(df) + 2))
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
  grDevices::pdf(file.path(method_dir, "cv_curve.pdf"), width = 7, height = 6)
  plot(cv)
  graphics::title("Cross-validation curve")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "cv_curve.png"), width = 1800, height = 1500, res = 240)
  plot(cv)
  graphics::title("Cross-validation curve")
  grDevices::dev.off()

  grDevices::pdf(file.path(method_dir, "coefficient_path.pdf"), width = 8, height = 6)
  plot(cv$glmnet.fit, xvar = "lambda", label = FALSE)
  graphics::abline(v = log(cv$lambda.min), lty = 2, col = "#D95F02", lwd = 2)
  graphics::abline(v = log(cv$lambda.1se), lty = 3, col = "#1B9E77", lwd = 2)
  graphics::legend("topright", c("lambda.min", "lambda.1se"), lty = c(2, 3), col = c("#D95F02", "#1B9E77"), bty = "n")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "coefficient_path.png"), width = 2000, height = 1500, res = 240)
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
  grDevices::pdf(file.path(method_dir, "oob_error_curve.pdf"), width = 7, height = 6)
  plot(fit, main = "Random forest OOB error")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "oob_error_curve.png"), width = 1800, height = 1500, res = 240)
  plot(fit, main = "Random forest OOB error")
  grDevices::dev.off()

  grDevices::pdf(file.path(method_dir, "variable_importance_dotchart.pdf"), width = 8, height = 7)
  randomForest::varImpPlot(fit, main = "Random forest variable importance")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "variable_importance_dotchart.png"), width = 1800, height = 1600, res = 240)
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
  grDevices::pdf(file.path(method_dir, "boruta_importance_boxplot.pdf"), width = 10, height = 7)
  plot(fixed, las = 2, cex.axis = 0.8, main = "Boruta feature importance")
  grDevices::dev.off()

  grDevices::png(file.path(method_dir, "boruta_importance_boxplot.png"), width = 2400, height = 1800, res = 240)
  plot(fixed, las = 2, cex.axis = 0.8, main = "Boruta feature importance")
  grDevices::dev.off()

  grDevices::pdf(file.path(method_dir, "boruta_importance_history.pdf"), width = 10, height = 7)
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
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height, bg = "white")
  ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height, dpi = 300, bg = "white")
}
