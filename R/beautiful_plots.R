# AutoML4R beautiful plotting module
# Run after auto_ml_analysis().
#
# Example:
#   library(AutoML4R)
#   result <- auto_ml_analysis(hub_data, group, methods = aml_methods())
#   automl4r_beautiful_plots(result, hub_data = hub_data, group = group)

automl4r_require_plot_pkgs <- function() {
  pkgs <- c("dplyr", "ggplot2", "stringr", "scales", "patchwork", "grid")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Please install required plotting packages first: ",
      paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

automl4r_method_label <- function(method) {
  labels <- c(
    univariate = "Wilcoxon",
    roc_auc = "ROC AUC",
    lasso = "LASSO",
    ridge = "Ridge",
    elastic_net = "ElasticNet",
    stability_selection = "StabilitySelection glmnet",
    random_forest = "RandomForest",
    ranger = "Ranger RF Permutation",
    boruta = "Boruta",
    vsurf_interpretation = "VSURF interpretation",
    vsurf_prediction = "VSURF prediction",
    xgboost = "XGBoost",
    lightgbm = "LightGBM",
    catboost = "CatBoost",
    mrmr = "mRMR",
    splsda = "sPLSDA",
    pam = "PAM",
    svm = "SVM",
    svm_rfe = "SVM RFE",
    relief = "ReliefF",
    information_gain = "InformationGain",
    gain_ratio = "GainRatio",
    symmetrical_uncertainty = "SymmetricalUncertainty",
    genetic_algorithm = "GeneticAlgorithm",
    simulated_annealing = "SimulatedAnnealing"
  )
  out <- labels[method]
  out[is.na(out)] <- method[is.na(out)]
  unname(out)
}

automl4r_method_class <- function(method) {
  dplyr::case_when(
    stringr::str_detect(method, stringr::regex("SHAP", ignore_case = TRUE)) ~ "Model interpretation",
    stringr::str_detect(method, stringr::regex("Wilcoxon|ROC|AUC|Ttest|DEG", ignore_case = TRUE)) ~ "Statistics",
    stringr::str_detect(method, stringr::regex("LASSO|Ridge|Enet|Elastic|Stability|PAM|glmnet", ignore_case = TRUE)) ~ "Regression/Stability",
    stringr::str_detect(method, stringr::regex("RF|Random|Boruta|Ranger|VSURF|XGBoost|LightGBM|CatBoost|GBM", ignore_case = TRUE)) ~ "Tree/Boosting",
    stringr::str_detect(method, stringr::regex("SVM|Svm|RFE", ignore_case = TRUE)) ~ "Kernel/RFE",
    stringr::str_detect(method, stringr::regex("mRMR|Information|Gain|Symmetrical|ReliefF|Relief", ignore_case = TRUE)) ~ "Information/Neighbor",
    stringr::str_detect(method, stringr::regex("sPLS|PLS|DIABLO", ignore_case = TRUE)) ~ "Latent-variable",
    stringr::str_detect(method, stringr::regex("Genetic|GA|Simulated|Annealing|SA", ignore_case = TRUE)) ~ "Heuristic search",
    TRUE ~ "Other"
  )
}

automl4r_theme_nature <- function(base_size = 15,
                                  font_panel_title = 24,
                                  font_axis_title = 20,
                                  font_axis_text = 18,
                                  font_legend_title = 18,
                                  font_legend_text = 18,
                                  font_facet = 24,
                                  font_family = "Arial") {
  ggplot2::theme_classic(base_size = base_size, base_family = font_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family, color = "black"),
      plot.title = ggplot2::element_text(family = font_family, size = font_panel_title, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = ggplot2::element_text(family = font_family, size = font_panel_title, color = "black"),
      axis.title = ggplot2::element_text(family = font_family, size = font_axis_title, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(family = font_family, size = font_axis_text, color = "black"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.55),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.45),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.65),
      legend.title = ggplot2::element_text(family = font_family, size = font_legend_title, face = "bold", color = "black"),
      legend.text = ggplot2::element_text(family = font_family, size = font_legend_text, color = "black"),
      strip.background = ggplot2::element_rect(fill = "grey97", color = "grey40", linewidth = 0.45),
      strip.text = ggplot2::element_text(family = font_family, size = font_facet, face = "bold", color = "black"),
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

automl4r_apply_plot_typography <- function(plot,
                                           font_family = "Arial",
                                           axis_text_size = 18,
                                           axis_title_size = 20,
                                           title_size = 24) {
  plot +
    ggplot2::theme(
      text = ggplot2::element_text(family = font_family),
      plot.title = ggplot2::element_text(family = font_family, size = title_size, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = ggplot2::element_text(family = font_family, size = title_size, color = "black"),
      axis.title = ggplot2::element_text(family = font_family, size = axis_title_size, face = "bold", color = "black"),
      axis.title.x = ggplot2::element_text(family = font_family, size = axis_title_size, face = "bold", color = "black"),
      axis.title.y = ggplot2::element_text(family = font_family, size = axis_title_size, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(family = font_family, size = axis_text_size, color = "black"),
      strip.text = ggplot2::element_text(family = font_family, size = title_size, face = "bold", color = "black"),
      legend.title = ggplot2::element_text(family = font_family, size = axis_text_size, face = "bold", color = "black"),
      legend.text = ggplot2::element_text(family = font_family, size = axis_text_size, color = "black"),
      plot.tag = ggplot2::element_text(family = font_family, size = title_size, face = "bold", color = "black")
    )
}

automl4r_extract_genes_from_object <- function(x) {
  if (is.null(x)) {
    return(character(0))
  }
  if (is.data.frame(x)) {
    possible_cols <- c("Genes", "genes", "Gene", "gene", "Feature", "feature", "features",
                       "selected.fea", "selected_features", "ID", "id", "all", "All")
    gene_col <- intersect(possible_cols, colnames(x))[1]
    genes <- if (!is.na(gene_col)) x[[gene_col]] else x[[1]]
  } else if (is.list(x)) {
    if ("all" %in% names(x)) {
      genes <- x$all
    } else if ("Genes" %in% names(x)) {
      genes <- x$Genes
    } else if ("genes" %in% names(x)) {
      genes <- x$genes
    } else if ("gene" %in% names(x)) {
      genes <- x$gene
    } else {
      genes <- unlist(x)
    }
  } else {
    genes <- unlist(x)
  }
  genes <- unique(na.omit(as.character(genes)))
  genes <- genes[genes != ""]
  genes <- genes[genes != "NA"]
  genes <- genes[genes != "NULL"]
  genes <- genes[!grepl("^\\s*$", genes)]
  genes
}

automl4r_result_to_all_genes <- function(result) {
  if (is.null(result$selected_genes) || !is.list(result$selected_genes)) {
    stop("result must be an auto_ml_analysis result containing selected_genes.", call. = FALSE)
  }
  all_genes <- result$selected_genes
  names(all_genes) <- automl4r_method_label(names(all_genes))
  all_genes
}

automl4r_save_plot <- function(plot, file_base, width, height, dpi = 600) {
  plot <- automl4r_apply_plot_typography(plot)
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
}

automl4r_plot_method_selected_scores <- function(result, output_dir, top_n_per_method = 20) {
  if (is.null(result$method_results)) {
    return(invisible(NULL))
  }
  method_dir_root <- file.path(output_dir, "NatureStyle_algorithm_selected_gene_plots")
  dir.create(method_dir_root, recursive = TRUE, showWarnings = FALSE)

  for (method in names(result$method_results)) {
    res <- result$method_results[[method]]
    method_label <- automl4r_method_label(method)
    method_dir <- file.path(method_dir_root, method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

    selected <- head(automl4r_extract_genes_from_object(res$selected), top_n_per_method)
    if (length(selected) == 0) {
      next
    }

    if (!is.null(res$rankings) && all(c("feature", "score") %in% colnames(res$rankings))) {
      plot_df <- res$rankings[res$rankings$feature %in% selected, , drop = FALSE]
      plot_df <- plot_df[order(plot_df$rank), , drop = FALSE]
    } else {
      plot_df <- data.frame(feature = selected, score = rev(seq_along(selected)), rank = seq_along(selected))
    }
    if (nrow(plot_df) == 0) {
      next
    }
    plot_df$feature <- factor(plot_df$feature, levels = rev(plot_df$feature))
    plot_df$score <- as.numeric(plot_df$score)

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = score, y = feature)) +
      ggplot2::geom_segment(
        ggplot2::aes(x = 0, xend = score, yend = feature),
        linewidth = 1.25,
        color = "#A6CEE3",
        lineend = "round"
      ) +
      ggplot2::geom_point(shape = 21, fill = "#1F78B4", color = "black", size = 4.2, stroke = 0.45) +
      ggplot2::labs(
        x = "Algorithm-specific selection score",
        y = NULL,
        title = paste0(method_label, " selected genes")
      ) +
      automl4r_theme_nature(base_size = 13, font_panel_title = 24, font_axis_title = 20, font_axis_text = 18) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.30),
        panel.grid.major.y = ggplot2::element_blank(),
        legend.position = "none"
      )

    height <- max(5.5, 0.34 * nrow(plot_df) + 2.2)
    automl4r_save_plot(p, file.path(method_dir, "NatureStyle_Selected_gene_score_lollipop"), 9.5, height)
  }
  invisible(NULL)
}

automl4r_plot_method_selected_expression <- function(result,
                                                     hub_data,
                                                     group,
                                                     output_dir,
                                                     top_n_per_method = 12,
                                                     positive_class = NULL) {
  if (is.null(hub_data) || is.null(group) || is.null(result$selected_genes)) {
    return(invisible(NULL))
  }

  x <- as.data.frame(hub_data, check.names = FALSE)
  y <- droplevels(as.factor(group))
  if (length(y) != nrow(x)) {
    stop("length(group) must equal nrow(hub_data).", call. = FALSE)
  }

  method_dir_root <- file.path(output_dir, "NatureStyle_algorithm_selected_gene_plots")
  dir.create(method_dir_root, recursive = TRUE, showWarnings = FALSE)
  fill_values <- setNames(c("#4C78A8", "#E15759", "#59A14F", "#F28E2B")[seq_len(nlevels(y))], levels(y))
  if (!is.null(positive_class) && positive_class %in% names(fill_values)) {
    fill_values[positive_class] <- "#E15759"
  }

  for (method in names(result$selected_genes)) {
    selected <- head(intersect(automl4r_extract_genes_from_object(result$selected_genes[[method]]), colnames(x)), top_n_per_method)
    if (length(selected) == 0) {
      next
    }
    method_label <- automl4r_method_label(method)
    method_dir <- file.path(method_dir_root, method)
    dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

    plot_df <- do.call(rbind, lapply(selected, function(gene) {
      data.frame(gene = gene, group = y, expression = x[[gene]], stringsAsFactors = FALSE)
    }))
    plot_df$gene <- factor(plot_df$gene, levels = selected)
    plot_df$group <- factor(plot_df$group, levels = levels(y))

    p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = group, y = expression, fill = group, color = group)) +
      ggplot2::geom_violin(width = 0.82, alpha = 0.18, linewidth = 0.35, trim = FALSE) +
      ggplot2::geom_boxplot(width = 0.22, alpha = 0.78, outlier.shape = NA, linewidth = 0.35) +
      ggplot2::geom_jitter(width = 0.12, alpha = 0.62, size = 1.25, stroke = 0) +
      ggplot2::facet_wrap(stats::as.formula("~ gene"), scales = "free_y", ncol = 4) +
      ggplot2::scale_fill_manual(values = fill_values) +
      ggplot2::scale_color_manual(values = fill_values) +
      ggplot2::labs(
        x = NULL,
        y = "Expression",
        title = paste0(method_label, " selected-gene expression")
      ) +
      automl4r_theme_nature(base_size = 12, font_panel_title = 24, font_axis_title = 20, font_axis_text = 18, font_facet = 24) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(family = "Arial", size = 18, angle = 35, hjust = 1),
        panel.grid.major.x = ggplot2::element_blank(),
        legend.position = "top"
      )

    height <- max(6, 2.35 * ceiling(length(selected) / 4))
    automl4r_save_plot(p, file.path(method_dir, "NatureStyle_Selected_gene_expression_boxviolin"), 12, height)
  }
  invisible(NULL)
}

automl4r_plot_lasso_optimized <- function(result, output_dir) {
  lasso_res <- result$method_results$lasso
  if (is.null(lasso_res) || is.null(lasso_res$details$cv)) {
    return(invisible(NULL))
  }

  cvfit_lasso <- lasso_res$details$cv
  fit_lasso <- cvfit_lasso$glmnet.fit
  method_dir <- file.path(output_dir, "NatureStyle_algorithm_selected_gene_plots", "lasso")
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  coef_min <- stats::coef(cvfit_lasso, s = "lambda.min")
  coef_mat_min <- as.matrix(coef_min)
  lasso_coef_df <- data.frame(
    Gene = rownames(coef_mat_min),
    Coefficient = as.numeric(coef_mat_min[, 1]),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(Gene != "(Intercept)") |>
    dplyr::mutate(AbsCoefficient = abs(Coefficient)) |>
    dplyr::arrange(dplyr::desc(AbsCoefficient))

  lasso_selected_df <- lasso_coef_df |>
    dplyr::filter(Coefficient != 0) |>
    dplyr::arrange(dplyr::desc(AbsCoefficient))

  lasso_genes <- lasso_selected_df$Gene
  utils::write.table(
    lasso_coef_df,
    file.path(method_dir, "LASSO_All_Coefficients.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  utils::write.table(
    lasso_selected_df,
    file.path(method_dir, "LASSO_Selected_Coefficients.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )

  coef_path_mat <- as.matrix(fit_lasso$beta)
  coef_path_df <- data.frame(
    log_lambda = rep(log(fit_lasso$lambda), times = nrow(coef_path_mat)),
    Gene = rep(rownames(coef_path_mat), each = length(fit_lasso$lambda)),
    Coefficient = as.vector(t(coef_path_mat)),
    stringsAsFactors = FALSE
  )
  coef_path_df$Selected <- ifelse(coef_path_df$Gene %in% lasso_genes, "Selected", "Other")
  coef_path_selected <- coef_path_df[coef_path_df$Gene %in% lasso_genes, , drop = FALSE]

  p_lasso_path <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = coef_path_df,
      ggplot2::aes(x = log_lambda, y = Coefficient, group = Gene),
      color = "grey75",
      linewidth = 0.4,
      alpha = 0.8
    ) +
    ggplot2::geom_line(
      data = coef_path_selected,
      ggplot2::aes(x = log_lambda, y = Coefficient, color = Gene, group = Gene),
      linewidth = 1.0
    ) +
    ggplot2::geom_vline(xintercept = log(cvfit_lasso$lambda.min), linetype = "dashed", linewidth = 0.8, color = "#E41A1C") +
    ggplot2::geom_vline(xintercept = log(cvfit_lasso$lambda.1se), linetype = "dashed", linewidth = 0.8, color = "#377EB8") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(x = "Log(lambda)", y = "Coefficients", title = "LASSO coefficient path") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = 13, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black"),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11, color = "black"),
      legend.position = ifelse(length(lasso_genes) <= 15, "right", "none"),
      plot.margin = ggplot2::margin(10, 20, 10, 10, unit = "pt")
    )

  automl4r_save_plot(p_lasso_path, file.path(method_dir, "Model_LASSO_Coefficient_Path_optimized"), 10, 8)
  automl4r_save_plot(p_lasso_path + ggplot2::theme(legend.position = "none"), file.path(method_dir, "Model_LASSO_Coefficient_Path_NoLegend"), 9, 7)

  cv_df <- data.frame(
    lambda = cvfit_lasso$lambda,
    log_lambda = log(cvfit_lasso$lambda),
    cvm = cvfit_lasso$cvm,
    cvsd = cvfit_lasso$cvsd,
    nzero = cvfit_lasso$nzero
  ) |>
    dplyr::mutate(ymin = cvm - cvsd, ymax = cvm + cvsd)

  cv_y_max <- max(cv_df$ymax, na.rm = TRUE)
  cv_y_min <- min(cv_df$ymin, na.rm = TRUE)
  cv_y_range <- cv_y_max - cv_y_min
  cv_label_y1 <- cv_y_max - 0.08 * cv_y_range
  cv_label_y2 <- cv_y_max - 0.18 * cv_y_range
  y_axis_label <- ifelse(!is.null(cvfit_lasso$name), cvfit_lasso$name, "Cross-validation error")

  p_lasso_cv <- ggplot2::ggplot(cv_df, ggplot2::aes(x = log_lambda, y = cvm)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.05, linewidth = 0.5, color = "grey50") +
    ggplot2::geom_point(size = 2.2, color = "black") +
    ggplot2::geom_vline(xintercept = log(cvfit_lasso$lambda.min), linetype = "dashed", linewidth = 0.9, color = "#E41A1C") +
    ggplot2::geom_vline(xintercept = log(cvfit_lasso$lambda.1se), linetype = "dashed", linewidth = 0.9, color = "#377EB8") +
    ggplot2::annotate("text", x = log(cvfit_lasso$lambda.min), y = cv_label_y1, label = paste0("lambda.min\n", signif(cvfit_lasso$lambda.min, 4)), color = "#E41A1C", size = 4, hjust = -0.05) +
    ggplot2::annotate("text", x = log(cvfit_lasso$lambda.1se), y = cv_label_y2, label = paste0("lambda.1se\n", signif(cvfit_lasso$lambda.1se, 4)), color = "#377EB8", size = 4, hjust = -0.05) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::labs(x = "Log(lambda)", y = y_axis_label, title = "Cross-validation for LASSO") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = 13, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black"),
      plot.margin = ggplot2::margin(10, 20, 10, 10, unit = "pt")
    )

  automl4r_save_plot(p_lasso_cv, file.path(method_dir, "Model_LASSO_CV_Curve_optimized"), 8, 7)

  if (nrow(lasso_selected_df) > 0) {
    p_lasso_bar <- ggplot2::ggplot(
      lasso_selected_df,
      ggplot2::aes(x = reorder(Gene, AbsCoefficient), y = AbsCoefficient)
    ) +
      ggplot2::geom_col(width = 0.75, fill = "#D55E00") +
      ggplot2::coord_flip() +
      ggplot2::theme_classic(base_size = 16) +
      ggplot2::labs(x = NULL, y = "Absolute LASSO coefficient", title = "LASSO-selected genes") +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
        axis.text.y = ggplot2::element_text(family = "Arial", size = 18, color = "black"),
        axis.text.x = ggplot2::element_text(family = "Arial", size = 18, color = "black"),
        axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
        axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black"),
        plot.margin = ggplot2::margin(10, 20, 10, 10, unit = "pt")
      )

    automl4r_save_plot(
      p_lasso_bar,
      file.path(method_dir, "Model_LASSO_Selected_Genes_Barplot"),
      10,
      max(6, nrow(lasso_selected_df) * 0.35)
    )
  }

  invisible(NULL)
}

automl4r_plot_ranked_algorithm_bar <- function(rankings,
                                               selected = NULL,
                                               output_dir,
                                               method,
                                               title,
                                               xlab = "Feature score",
                                               file_stem = "Model_Selected_Features",
                                               fill = "#1F78B4",
                                               top_n = 30) {
  if (is.null(rankings) || !all(c("feature", "score") %in% colnames(rankings))) {
    return(invisible(NULL))
  }
  plot_df <- rankings
  if (!is.null(selected)) {
    selected <- automl4r_extract_genes_from_object(selected)
    keep <- plot_df$feature %in% selected
    if (any(keep)) {
      plot_df <- plot_df[keep, , drop = FALSE]
    }
  }
  plot_df <- plot_df[is.finite(plot_df$score), , drop = FALSE]
  plot_df <- plot_df[order(-abs(plot_df$score), plot_df$feature), , drop = FALSE]
  plot_df <- head(plot_df, min(top_n, nrow(plot_df)))
  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }
  plot_df$feature <- factor(plot_df$feature, levels = rev(plot_df$feature))
  plot_df$score_abs <- abs(plot_df$score)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = feature, y = score_abs)) +
    ggplot2::geom_col(width = 0.72, fill = fill, color = "black", linewidth = 0.25) +
    ggplot2::geom_point(color = "#B2182B", size = 2.6) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = xlab, title = title) +
    automl4r_theme_nature(base_size = 13, font_panel_title = 24, font_axis_title = 20, font_axis_text = 18) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.30),
      panel.grid.major.y = ggplot2::element_blank(),
      legend.position = "none"
    )

  method_dir <- file.path(output_dir, "NatureStyle_algorithm_selected_gene_plots", method)
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)
  automl4r_save_plot(p, file.path(method_dir, file_stem), 10, max(6, 0.35 * nrow(plot_df)))
  invisible(p)
}

automl4r_plot_glmnet_optimized <- function(result, output_dir, method) {
  res <- result$method_results[[method]]
  if (is.null(res) || is.null(res$details$cv)) {
    return(invisible(NULL))
  }
  cvfit <- res$details$cv
  fit <- cvfit$glmnet.fit
  method_label <- automl4r_method_label(method)
  method_dir <- file.path(output_dir, "NatureStyle_algorithm_selected_gene_plots", method)
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  coef_min <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
  coef_df <- data.frame(
    Gene = rownames(coef_min),
    Coefficient = as.numeric(coef_min[, 1]),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(Gene != "(Intercept)") |>
    dplyr::mutate(AbsCoefficient = abs(Coefficient)) |>
    dplyr::arrange(dplyr::desc(AbsCoefficient))
  selected_df <- coef_df |>
    dplyr::filter(Coefficient != 0) |>
    dplyr::arrange(dplyr::desc(AbsCoefficient))
  if (nrow(selected_df) == 0) {
    selected_df <- head(coef_df, min(30, nrow(coef_df)))
  }
  selected_genes <- selected_df$Gene

  utils::write.table(coef_df, file.path(method_dir, paste0(method_label, "_All_Coefficients.txt")), row.names = FALSE, quote = FALSE, sep = "\t")
  utils::write.table(selected_df, file.path(method_dir, paste0(method_label, "_Selected_Coefficients.txt")), row.names = FALSE, quote = FALSE, sep = "\t")

  coef_path_mat <- as.matrix(fit$beta)
  coef_path_df <- data.frame(
    log_lambda = rep(log(fit$lambda), times = nrow(coef_path_mat)),
    Gene = rep(rownames(coef_path_mat), each = length(fit$lambda)),
    Coefficient = as.vector(t(coef_path_mat)),
    stringsAsFactors = FALSE
  )
  coef_path_selected <- coef_path_df[coef_path_df$Gene %in% selected_genes, , drop = FALSE]

  p_path <- ggplot2::ggplot() +
    ggplot2::geom_line(
      data = coef_path_df,
      ggplot2::aes(x = log_lambda, y = Coefficient, group = Gene),
      color = "grey76",
      linewidth = 0.35,
      alpha = 0.80
    ) +
    ggplot2::geom_line(
      data = coef_path_selected,
      ggplot2::aes(x = log_lambda, y = Coefficient, color = Gene, group = Gene),
      linewidth = 0.95
    ) +
    ggplot2::geom_vline(xintercept = log(cvfit$lambda.min), linetype = "dashed", linewidth = 0.85, color = "#E41A1C") +
    ggplot2::geom_vline(xintercept = log(cvfit$lambda.1se), linetype = "dashed", linewidth = 0.85, color = "#377EB8") +
    ggplot2::labs(x = "Log(lambda)", y = "Coefficients", title = paste0(method_label, " coefficient path")) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = 13, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black"),
      legend.title = ggplot2::element_blank(),
      legend.position = ifelse(length(selected_genes) <= 15, "right", "none")
    )
  automl4r_save_plot(p_path, file.path(method_dir, paste0("Model_", method_label, "_Coefficient_Path_optimized")), 10, 8)

  cv_df <- data.frame(
    log_lambda = log(cvfit$lambda),
    cvm = cvfit$cvm,
    cvsd = cvfit$cvsd,
    stringsAsFactors = FALSE
  ) |> dplyr::mutate(ymin = cvm - cvsd, ymax = cvm + cvsd)
  cv_y_max <- max(cv_df$ymax, na.rm = TRUE)
  cv_y_min <- min(cv_df$ymin, na.rm = TRUE)
  cv_y_range <- cv_y_max - cv_y_min
  y_axis_label <- ifelse(!is.null(cvfit$name), cvfit$name, "Cross-validation error")
  p_cv <- ggplot2::ggplot(cv_df, ggplot2::aes(x = log_lambda, y = cvm)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = ymin, ymax = ymax), width = 0.05, linewidth = 0.5, color = "grey50") +
    ggplot2::geom_point(size = 2.2, color = "black") +
    ggplot2::geom_vline(xintercept = log(cvfit$lambda.min), linetype = "dashed", linewidth = 0.9, color = "#E41A1C") +
    ggplot2::geom_vline(xintercept = log(cvfit$lambda.1se), linetype = "dashed", linewidth = 0.9, color = "#377EB8") +
    ggplot2::annotate("text", x = log(cvfit$lambda.min), y = cv_y_max - 0.08 * cv_y_range, label = paste0("lambda.min\n", signif(cvfit$lambda.min, 4)), color = "#E41A1C", size = 4, hjust = -0.05) +
    ggplot2::annotate("text", x = log(cvfit$lambda.1se), y = cv_y_max - 0.18 * cv_y_range, label = paste0("lambda.1se\n", signif(cvfit$lambda.1se, 4)), color = "#377EB8", size = 4, hjust = -0.05) +
    ggplot2::labs(x = "Log(lambda)", y = y_axis_label, title = paste0("Cross-validation for ", method_label)) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = 13, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black")
    )
  automl4r_save_plot(p_cv, file.path(method_dir, paste0("Model_", method_label, "_CV_Curve_optimized")), 8, 7)

  p_abs <- ggplot2::ggplot(selected_df, ggplot2::aes(x = reorder(Gene, AbsCoefficient), y = AbsCoefficient)) +
    ggplot2::geom_col(width = 0.75, fill = ifelse(method == "ridge", "#0072B2", "#009E73")) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = paste0("Absolute ", method_label, " coefficient"), title = paste0(method_label, "-selected genes")) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = 13, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black")
    )
  automl4r_save_plot(p_abs, file.path(method_dir, paste0("Model_", method_label, "_Selected_Genes_AbsCoefficient")), 10, max(6, nrow(selected_df) * 0.35))

  p_signed <- ggplot2::ggplot(selected_df, ggplot2::aes(x = reorder(Gene, Coefficient), y = Coefficient, fill = Coefficient > 0)) +
    ggplot2::geom_col(width = 0.75, color = "black", linewidth = 0.2) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "#377EB8", `TRUE` = "#D55E00")) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = paste0(method_label, " coefficient"), title = paste0(method_label, " signed coefficients")) +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = 13, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black"),
      legend.position = "none"
    )
  automl4r_save_plot(p_signed, file.path(method_dir, paste0("Model_", method_label, "_Selected_Genes_SignedCoefficient")), 10, max(6, nrow(selected_df) * 0.35))
  invisible(NULL)
}

automl4r_plot_random_forest_optimized <- function(result, output_dir) {
  res <- result$method_results$random_forest
  fit <- res$details$model
  if (is.null(res) || is.null(fit) || is.null(res$rankings)) {
    return(invisible(NULL))
  }
  method_dir <- file.path(output_dir, "NatureStyle_algorithm_selected_gene_plots", "random_forest")
  dir.create(method_dir, recursive = TRUE, showWarnings = FALSE)

  rf_data_plot <- res$rankings |>
    dplyr::mutate(relative_imp = score / max(score, na.rm = TRUE), ID = feature) |>
    dplyr::arrange(relative_imp)
  rf_data_plot$ID <- factor(rf_data_plot$ID, levels = rf_data_plot$ID)
  importance_cutoff <- 0.5
  p_rf_imp <- ggplot2::ggplot(rf_data_plot, ggplot2::aes(x = ID, y = relative_imp)) +
    ggplot2::geom_col(width = 0.65, fill = "#1E90FF", color = "black", linewidth = 0.2) +
    ggplot2::geom_point(color = "#B22222", size = 2.8) +
    ggplot2::geom_hline(yintercept = importance_cutoff, linetype = "dashed", color = "red", linewidth = 0.7) +
    ggplot2::coord_flip() +
    ggplot2::labs(x = NULL, y = "Relative importance", title = "Random forest feature importance") +
    ggplot2::theme_classic(base_size = 16) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
      axis.title.x = ggplot2::element_text(size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(size = 13, color = "black"),
      axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
      axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black")
    )
  automl4r_save_plot(p_rf_imp, file.path(method_dir, "Model_RF_Feature_Importance"), 10, max(6, nrow(rf_data_plot) * 0.32))

  if (!is.null(fit$err.rate)) {
    oob_col <- if ("OOB" %in% colnames(fit$err.rate)) "OOB" else colnames(fit$err.rate)[1]
    oob_df <- data.frame(trees = seq_len(nrow(fit$err.rate)), error = fit$err.rate[, oob_col])
    best_tree <- oob_df$trees[which.min(oob_df$error)]
    best_error <- min(oob_df$error, na.rm = TRUE)
    p_oob <- ggplot2::ggplot(oob_df, ggplot2::aes(x = trees, y = error)) +
      ggplot2::geom_line(linewidth = 0.7, color = "gray35") +
      ggplot2::geom_vline(xintercept = best_tree, linetype = "dashed", color = "red", linewidth = 0.8) +
      ggplot2::geom_point(data = oob_df[oob_df$trees == best_tree, , drop = FALSE], color = "red", size = 3) +
      ggplot2::annotate("text", x = best_tree, y = best_error, label = paste0("Best tree = ", best_tree), color = "red", size = 4, hjust = -0.05, vjust = -0.8) +
      ggplot2::labs(x = "Number of trees", y = "OOB error", title = "Random forest OOB error curve") +
      ggplot2::theme_classic(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, color = "black"),
        axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
        axis.text = ggplot2::element_text(size = 13, color = "black"),
        axis.line = ggplot2::element_line(linewidth = 0.7, color = "black"),
        axis.ticks = ggplot2::element_line(linewidth = 0.7, color = "black")
      )
    automl4r_save_plot(p_oob, file.path(method_dir, "Model_RF_OOB_Error_Curve"), 9, 7)
  }
  invisible(NULL)
}

automl4r_plot_model_specific_figures <- function(result, output_dir, top_n = 30) {
  if (is.null(result$method_results)) {
    return(invisible(NULL))
  }
  try(automl4r_plot_lasso_optimized(result, output_dir), silent = TRUE)
  for (method in intersect(c("ridge", "elastic_net"), names(result$method_results))) {
    try(automl4r_plot_glmnet_optimized(result, output_dir, method), silent = TRUE)
  }
  try(automl4r_plot_random_forest_optimized(result, output_dir), silent = TRUE)

  plot_specs <- list(
    ranger = list(title = "Ranger random forest permutation importance", xlab = "Permutation importance", fill = "#0072B2", stem = "Model_Ranger_RF_Permutation_Importance"),
    boruta = list(title = "Boruta confirmed feature importance", xlab = "Median importance", fill = "#D55E00", stem = "Model_Boruta_Selected_Features_Barplot"),
    xgboost = list(title = "XGBoost gain importance", xlab = "Gain importance", fill = "#1F78B4", stem = "Model_XGBoost_Gain_Importance_optimized"),
    lightgbm = list(title = "LightGBM gain importance", xlab = "Gain importance", fill = "#59A14F", stem = "Model_LightGBM_Gain"),
    catboost = list(title = "CatBoost feature importance", xlab = "Feature importance", fill = "#F28E2B", stem = "Model_CatBoost_Feature_Importance_optimized"),
    mrmr = list(title = "mRMR selection rank", xlab = "Selection score", fill = "#76B7B2", stem = "Model_mRMR_Selection_Rank_optimized"),
    splsda = list(title = "sPLS-DA selected gene loadings", xlab = "Loading score", fill = "#9C755F", stem = "Model_sPLSDA_Selected_Gene_Loadings"),
    pam = list(title = "PAM selected gene score", xlab = "PAM score", fill = "#E15759", stem = "Model_PAM_Selected_Genes_Score"),
    svm = list(title = "SVM feature importance", xlab = "Variable importance", fill = "#B07AA1", stem = "Model_SVM_Selected_Features"),
    svm_rfe = list(title = "SVM-RFE selected genes", xlab = "Selection score", fill = "#B07AA1", stem = "Model_SVM_RFE_Selected_Genes_optimized"),
    relief = list(title = "ReliefF feature importance", xlab = "ReliefF score", fill = "#76B7B2", stem = "Model_ReliefF_Feature_Importance"),
    information_gain = list(title = "Information gain feature ranking", xlab = "Information gain", fill = "#4C78A8", stem = "InformationGain_ranking_plot"),
    gain_ratio = list(title = "Gain ratio feature ranking", xlab = "Gain ratio", fill = "#4C78A8", stem = "GainRatio_ranking_plot"),
    symmetrical_uncertainty = list(title = "Symmetrical uncertainty feature ranking", xlab = "Symmetrical uncertainty", fill = "#4C78A8", stem = "SymmetricalUncertainty_ranking_plot"),
    stability_selection = list(title = "Stability selection probability ranking", xlab = "Selection probability", fill = "#59A14F", stem = "Model_StabilitySelection_Probability_Ranking"),
    genetic_algorithm = list(title = "Genetic algorithm-selected feature importance", xlab = "Selection score", fill = "#EDC948", stem = "Model_GA_Selected_Features"),
    simulated_annealing = list(title = "Simulated annealing selected features", xlab = "Selection score", fill = "#EDC948", stem = "Model_SA_Selected_Features_Barplot_Normal"),
    vsurf_interpretation = list(title = "VSURF interpretation feature importance", xlab = "Selection score", fill = "#F28E2B", stem = "Model_VSURF_interpretation_Feature_Importance"),
    vsurf_prediction = list(title = "VSURF prediction selected features", xlab = "Selection score", fill = "#F28E2B", stem = "Model_VSURF_prediction_Feature_Importance"),
    roc_auc = list(title = "ROC AUC top-gene ranking", xlab = "Absolute AUC score", fill = "#4C78A8", stem = "Model_ROC_AUC_Selected_Genes"),
    univariate = list(title = "Wilcoxon differential feature ranking", xlab = "-log10 P-value", fill = "#4C78A8", stem = "Model_Wilcoxon_Pvalue_Ranking")
  )

  for (method in intersect(names(plot_specs), names(result$method_results))) {
    spec <- plot_specs[[method]]
    res <- result$method_results[[method]]
    try(
      automl4r_plot_ranked_algorithm_bar(
        rankings = res$rankings,
        selected = res$selected,
        output_dir = output_dir,
        method = method,
        title = spec$title,
        xlab = spec$xlab,
        file_stem = spec$stem,
        fill = spec$fill,
        top_n = top_n
      ),
      silent = TRUE
    )
  }

  invisible(NULL)
}

automl4r_plot_consensus_nature <- function(all_genes,
                                           output_dir,
                                           top_n_vote = 30,
                                           heatmap_top_n = 40,
                                           min_support_fraction = 0.30) {
  if (is.data.frame(all_genes)) {
    all_genes <- as.list(all_genes)
  }
  if (!is.list(all_genes)) {
    stop("all_genes must be a named list or data.frame.", call. = FALSE)
  }

  method_names <- names(all_genes)
  if (is.null(method_names)) {
    method_names <- paste0("Method_", seq_along(all_genes))
  }
  method_names <- make.unique(method_names)

  vote_long <- lapply(seq_along(all_genes), function(i) {
    genes <- automl4r_extract_genes_from_object(all_genes[[i]])
    if (length(genes) == 0) {
      return(NULL)
    }
    data.frame(method = method_names[i], gene = genes, stringsAsFactors = FALSE)
  })
  vote_long <- dplyr::bind_rows(vote_long)
  vote_long <- as.data.frame(vote_long, stringsAsFactors = FALSE)
  if (nrow(vote_long) == 0) {
    stop("No valid gene records were extracted from all_genes.", call. = FALSE)
  }

  vote_long <- vote_long |>
    dplyr::filter(!is.na(method), method != "", !is.na(gene), gene != "", gene != "NA", gene != "NULL") |>
    dplyr::distinct(method, gene)

  n_methods_total <- dplyr::n_distinct(vote_long$method)
  utils::write.csv(vote_long, file.path(output_dir, "All_methods_gene_long_table.csv"), row.names = FALSE)

  vote_table <- vote_long |>
    dplyr::group_by(gene) |>
    dplyr::summarise(
      n_methods = dplyr::n_distinct(method),
      support_rate = n_methods / n_methods_total,
      methods = paste(sort(unique(method)), collapse = "; "),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_methods), gene)

  utils::write.csv(vote_table, file.path(output_dir, "Consensus_gene_vote_table.csv"), row.names = FALSE)

  min_support_methods <- min(n_methods_total, max(3, ceiling(min_support_fraction * n_methods_total)))
  consensus_genes <- vote_table |>
    dplyr::filter(n_methods >= min_support_methods) |>
    dplyr::pull(gene)
  utils::write.table(
    data.frame(Genes = consensus_genes),
    file.path(output_dir, "Consensus_Genes_selected_by_voting.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )

  method_info <- vote_long |>
    dplyr::distinct(method) |>
    dplyr::mutate(method_class = automl4r_method_class(method))
  utils::write.csv(method_info, file.path(output_dir, "Feature_selection_method_class_annotation.csv"), row.names = FALSE)

  font_main_title <- 24
  font_panel_tag <- 24
  font_panel_title <- 24
  font_axis_title <- 20
  font_axis_text <- 18
  font_gene_A <- 18
  font_gene_B <- 18
  font_gene_C <- 18
  font_heatmap_x <- 18
  font_legend_title <- 18
  font_legend_text <- 18

  pal_low <- "#A6CEE3"
  pal_high <- "#1F78B4"
  pal_cutoff <- "#B2182B"

  method_class_palette <- c(
    "Statistics" = "#4C78A8",
    "Regression/Stability" = "#59A14F",
    "Tree/Boosting" = "#F28E2B",
    "Kernel/RFE" = "#B07AA1",
    "Model interpretation" = "#E15759",
    "Information/Neighbor" = "#76B7B2",
    "Latent-variable" = "#9C755F",
    "Heuristic search" = "#EDC948",
    "Other" = "#8C8C8C",
    "Not selected" = "grey94"
  )

  vote_plot_df <- vote_table |>
    dplyr::arrange(dplyr::desc(n_methods), gene) |>
    dplyr::slice_head(n = min(top_n_vote, nrow(vote_table))) |>
    dplyr::distinct(gene, .keep_all = TRUE) |>
    dplyr::arrange(n_methods, gene) |>
    dplyr::mutate(gene = factor(gene, levels = gene))

  x_break_by <- max(1, ceiling(n_methods_total / 8))
  p_vote <- ggplot2::ggplot(vote_plot_df, ggplot2::aes(x = n_methods, y = gene)) +
    ggplot2::geom_segment(
      ggplot2::aes(x = 0, xend = n_methods, y = gene, yend = gene, color = support_rate),
      linewidth = 1.35,
      lineend = "round"
    ) +
    ggplot2::geom_point(ggplot2::aes(fill = support_rate), shape = 21, color = "black", size = 4.0, stroke = 0.45) +
    ggplot2::geom_vline(xintercept = min_support_methods, linetype = "dashed", color = pal_cutoff, linewidth = 0.75) +
    ggplot2::scale_color_gradient(low = pal_low, high = pal_high, labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_gradient(low = pal_low, high = pal_high, labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_x_continuous(
      breaks = seq(0, n_methods_total, by = x_break_by),
      limits = c(0, n_methods_total + 0.6),
      expand = ggplot2::expansion(mult = c(0, 0.02))
    ) +
    ggplot2::labs(
      x = "Number of supporting algorithms",
      y = NULL,
      color = "Support rate",
      fill = "Support rate",
      title = "Consensus ranking of candidate genes"
    ) +
    automl4r_theme_nature(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(family = "Arial", size = font_axis_text, color = "black"),
      axis.text.y = ggplot2::element_text(family = "Arial", size = font_gene_A, color = "black"),
      legend.position = "right",
      legend.key.height = grid::unit(0.70, "cm"),
      legend.key.width = grid::unit(0.45, "cm"),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.30),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::guides(
      color = "none",
      fill = ggplot2::guide_colorbar(title.position = "top", barwidth = grid::unit(0.45, "cm"), barheight = grid::unit(3.2, "cm"))
    )

  method_count <- vote_long |>
    dplyr::count(method, name = "n_selected_genes") |>
    dplyr::left_join(method_info, by = "method")
  utils::write.csv(method_count, file.path(output_dir, "Algorithm_selected_gene_count_table.csv"), row.names = FALSE)

  class_order <- method_count |>
    dplyr::group_by(method_class) |>
    dplyr::summarise(
      n_algorithms = dplyr::n_distinct(method),
      total_selected_genes = sum(n_selected_genes, na.rm = TRUE),
      mean_selected_genes = mean(n_selected_genes, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_algorithms), dplyr::desc(total_selected_genes), dplyr::desc(mean_selected_genes), method_class) |>
    dplyr::pull(method_class)

  method_count <- method_count |>
    dplyr::mutate(method_class = factor(method_class, levels = class_order)) |>
    dplyr::arrange(method_class, dplyr::desc(n_selected_genes), method)
  method_order <- method_count$method

  heatmap_genes <- vote_table |>
    dplyr::arrange(dplyr::desc(n_methods), gene) |>
    dplyr::pull(gene)
  heatmap_genes <- heatmap_genes[seq_len(min(heatmap_top_n, length(heatmap_genes)))]

  gene_order <- vote_table |>
    dplyr::filter(gene %in% heatmap_genes) |>
    dplyr::arrange(dplyr::desc(n_methods), gene) |>
    dplyr::pull(gene)

  heatmap_long <- expand.grid(gene = heatmap_genes, method = method_order, stringsAsFactors = FALSE) |>
    dplyr::left_join(vote_long |> dplyr::mutate(selected = 1), by = c("gene", "method")) |>
    dplyr::mutate(selected = ifelse(is.na(selected), 0, selected)) |>
    dplyr::left_join(method_info, by = "method") |>
    dplyr::mutate(
      method_class = factor(method_class, levels = class_order),
      tile_fill = ifelse(selected == 1, as.character(method_class), "Not selected"),
      gene = factor(gene, levels = rev(gene_order)),
      method = factor(method, levels = method_order)
    )

  p_heatmap <- ggplot2::ggplot(heatmap_long, ggplot2::aes(x = method, y = gene, fill = tile_fill)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::facet_grid(. ~ method_class, scales = "free_x", space = "free_x") +
    ggplot2::scale_fill_manual(values = method_class_palette, breaks = class_order, name = "Algorithm class") +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
    ggplot2::labs(x = "Feature-selection algorithms grouped by algorithm class", y = NULL, title = "Algorithm-level support for candidate genes") +
    automl4r_theme_nature(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(family = "Arial", size = font_panel_title, face = "bold", hjust = 0.5, color = "black"),
      axis.title.x = ggplot2::element_text(family = "Arial", size = font_axis_title, face = "bold", color = "black"),
      axis.text.x = ggplot2::element_text(family = "Arial", size = font_heatmap_x, angle = 50, hjust = 1, vjust = 1, color = "black"),
      axis.text.y = ggplot2::element_text(family = "Arial", size = font_gene_B, color = "black"),
      strip.text.x = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.border = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(family = "Arial", size = font_legend_title, face = "bold", color = "black"),
      legend.text = ggplot2::element_text(family = "Arial", size = font_legend_text, color = "black"),
      legend.key.height = grid::unit(0.60, "cm"),
      legend.key.width = grid::unit(0.60, "cm"),
      panel.spacing.x = grid::unit(0.08, "lines"),
      plot.margin = ggplot2::margin(8, 8, 8, 8)
    )

  method_count_plot <- method_count |>
    dplyr::mutate(method = factor(method, levels = rev(method_order)))

  p_method_count <- ggplot2::ggplot(method_count_plot, ggplot2::aes(x = n_selected_genes, y = method, fill = method_class)) +
    ggplot2::geom_col(width = 0.70, color = "black", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = method_class_palette, drop = FALSE, name = "Algorithm class") +
    ggplot2::scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 26)) +
    ggplot2::labs(x = "Number of selected genes", y = NULL, title = "Feature output size of individual algorithms") +
    automl4r_theme_nature(base_size = 13) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(family = "Arial", size = font_axis_text, color = "black"),
      axis.text.y = ggplot2::element_text(family = "Arial", size = font_gene_C, color = "black"),
      legend.position = "right",
      legend.key.height = grid::unit(0.60, "cm"),
      legend.key.width = grid::unit(0.60, "cm"),
      panel.grid.major.x = ggplot2::element_line(color = "grey90", linewidth = 0.30),
      panel.grid.major.y = ggplot2::element_blank()
    )

  n_heatmap_genes <- length(heatmap_genes)
  vote_height <- max(5.5, 0.30 * nrow(vote_plot_df) + 1.8)
  heatmap_width <- max(24, min(46, 0.90 * n_methods_total + 12))
  heatmap_height <- max(7.5, 0.30 * n_heatmap_genes + 3.0)
  method_count_width <- 11
  method_count_height <- max(6.5, 0.30 * nrow(method_count_plot) + 2.0)
  combined_width <- max(24, min(48, 0.92 * n_methods_total + 13))
  combined_height <- max(22, 0.38 * n_heatmap_genes + 12)

  automl4r_save_plot(p_vote, file.path(output_dir, "NatureStyle_Consensus_gene_vote_lollipop_largefont"), 10, vote_height)
  automl4r_save_plot(p_heatmap, file.path(output_dir, "NatureStyle_Method_gene_consensus_heatmap_class_size_ordered_largefont"), heatmap_width, heatmap_height)
  automl4r_save_plot(p_method_count, file.path(output_dir, "NatureStyle_Algorithm_selected_gene_count_class_size_ordered_largefont"), method_count_width, method_count_height)

  p_combined_final <- p_vote / p_heatmap / p_method_count +
    patchwork::plot_layout(heights = c(1.00, 1.90, 1.5)) +
    patchwork::plot_annotation(
      tag_levels = "A",
      title = "Multi-algorithm consensus feature-selection overview",
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(family = "Arial", color = "black", size = font_main_title, face = "bold", hjust = 0.5),
        plot.tag = ggplot2::element_text(family = "Arial", color = "black", size = font_panel_tag, face = "bold")
      )
    )

  ggplot2::ggsave(
    filename = file.path(output_dir, "NatureStyle_Consensus_feature_selection_overview_FINAL_largefont.pdf"),
    plot = automl4r_apply_plot_typography(p_combined_final),
    width = combined_width + 12,
    height = combined_height + 20,
    bg = "white",
    device = grDevices::cairo_pdf
  )

  list(
    vote_long = vote_long,
    vote_table = vote_table,
    method_info = method_info,
    method_count = method_count,
    consensus_genes = consensus_genes,
    plots = list(
      vote = p_vote,
      heatmap = p_heatmap,
      method_count = p_method_count,
      combined = p_combined_final
    )
  )
}

automl4r_beautiful_plots <- function(result,
                                     hub_data = NULL,
                                     group = NULL,
                                     output_dir = NULL,
                                     all_genes = NULL,
                                     top_n_vote = 30,
                                     heatmap_top_n = 40,
                                     top_n_per_method = 20,
                                     expression_top_n_per_method = 12,
                                     min_support_fraction = 0.30,
                                     positive_class = NULL) {
  automl4r_require_plot_pkgs()

  if (is.null(output_dir)) {
    output_dir <- result$output_dir
  }
  if (is.null(output_dir) || length(output_dir) == 0 || is.na(output_dir)) {
    output_dir <- "ML_screening_results"
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  if (is.null(all_genes)) {
    all_genes <- automl4r_result_to_all_genes(result)
  }

  consensus_out <- automl4r_plot_consensus_nature(
    all_genes = all_genes,
    output_dir = output_dir,
    top_n_vote = top_n_vote,
    heatmap_top_n = heatmap_top_n,
    min_support_fraction = min_support_fraction
  )

  automl4r_plot_method_selected_scores(result, output_dir, top_n_per_method = top_n_per_method)
  automl4r_plot_model_specific_figures(result, output_dir, top_n = top_n_per_method)
  automl4r_plot_method_selected_expression(
    result,
    hub_data = hub_data,
    group = group,
    output_dir = output_dir,
    top_n_per_method = expression_top_n_per_method,
    positive_class = positive_class
  )

  cat("Beautiful AutoML4R plots saved in:\n")
  cat(normalizePath(output_dir, mustWork = FALSE), "\n")
  cat("Main output file:\n")
  cat(file.path(output_dir, "NatureStyle_Consensus_feature_selection_overview_FINAL_largefont.pdf"), "\n")

  invisible(consensus_out)
}
