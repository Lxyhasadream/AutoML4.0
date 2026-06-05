# AutoML4R: robust automated feature-selection workflow.

auto_ml_analysis <- function(hub_data,
                             group,
                             output_dir = "ML_screening_results",
                             positive_class = NULL,
                             methods = c(
                               "univariate",
                               "roc_auc",
                               "lasso",
                               "ridge",
                               "elastic_net",
                               "random_forest",
                               "boruta",
                               "xgboost"
                             ),
                             top_n = 20,
                             seed = 123,
                             write_outputs = TRUE,
                             verbose = TRUE) {
  set.seed(seed)
  methods <- match.arg(
    methods,
    choices = c(
      "univariate",
      "roc_auc",
      "lasso",
      "ridge",
      "elastic_net",
      "random_forest",
      "boruta",
      "xgboost"
    ),
    several.ok = TRUE
  )
  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }

  prepared <- aml_prepare_input(hub_data, group, positive_class)
  x <- prepared$x
  y <- prepared$y
  y01 <- prepared$y01
  positive_class <- prepared$positive_class

  if (write_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  method_jobs <- list(
    univariate = function() aml_run_univariate(x, y, top_n),
    roc_auc = function() aml_run_roc_auc(x, y01, top_n),
    lasso = function() aml_run_glmnet(x, y, alpha = 1, top_n = top_n, seed = seed),
    ridge = function() aml_run_glmnet(x, y, alpha = 0, top_n = top_n, seed = seed),
    elastic_net = function() aml_run_glmnet(x, y, alpha = 0.5, top_n = top_n, seed = seed),
    random_forest = function() aml_run_random_forest(x, y, top_n, seed),
    boruta = function() aml_run_boruta(x, y, top_n, seed),
    xgboost = function() aml_run_xgboost(x, y01, top_n, seed)
  )

  results <- list()
  failures <- list()

  for (method in methods) {
    if (verbose) {
      message("Running ", method, "...")
    }
    out <- tryCatch(
      method_jobs[[method]](),
      error = function(e) e
    )
    if (inherits(out, "error")) {
      failures[[method]] <- conditionMessage(out)
      if (verbose) {
        message("Skipping ", method, ": ", conditionMessage(out))
      }
    } else {
      results[[method]] <- out
    }
  }

  if (length(results) == 0) {
    stop("No feature-selection method completed successfully.", call. = FALSE)
  }

  rankings <- aml_bind_rankings(results)
  consensus <- aml_make_consensus(results, n_methods = length(results))

  if (write_outputs) {
    utils::write.csv(rankings, file.path(output_dir, "feature_rankings_by_method.csv"), row.names = FALSE)
    utils::write.csv(consensus, file.path(output_dir, "consensus_feature_votes.csv"), row.names = FALSE)
    aml_write_selected_genes(results, output_dir)
    aml_plot_consensus(consensus, output_dir, top_n)
    aml_plot_method_heatmap(results, output_dir, top_n)
  }

  structure(
    list(
      selected_genes = lapply(results, function(z) z$selected),
      rankings = rankings,
      consensus = consensus,
      failures = failures,
      positive_class = positive_class,
      output_dir = if (write_outputs) normalizePath(output_dir, mustWork = FALSE) else NA_character_
    ),
    class = "auto_ml_analysis"
  )
}

aml_prepare_input <- function(hub_data, group, positive_class = NULL) {
  if (missing(hub_data) || missing(group)) {
    stop("hub_data and group are required.", call. = FALSE)
  }

  x <- as.data.frame(hub_data, check.names = FALSE)
  if (nrow(x) < 6 || ncol(x) < 2) {
    stop("hub_data must contain at least 6 samples and 2 features.", call. = FALSE)
  }
  if (length(group) != nrow(x)) {
    stop("length(group) must equal nrow(hub_data).", call. = FALSE)
  }
  if (anyNA(group)) {
    stop("group cannot contain NA values.", call. = FALSE)
  }

  y <- droplevels(as.factor(group))
  if (nlevels(y) != 2) {
    stop("This workflow currently supports binary classification only.", call. = FALSE)
  }
  if (is.null(positive_class)) {
    positive_class <- levels(y)[2]
  }
  if (!positive_class %in% levels(y)) {
    stop("positive_class must be one of: ", paste(levels(y), collapse = ", "), call. = FALSE)
  }

  x[] <- lapply(x, function(col) {
    if (is.factor(col) || is.character(col) || is.logical(col)) {
      col <- suppressWarnings(as.numeric(as.character(col)))
    }
    suppressWarnings(as.numeric(col))
  })
  bad_cols <- names(x)[!vapply(x, function(col) all(is.finite(col) | is.na(col)), logical(1))]
  if (length(bad_cols) > 0) {
    stop("Non-numeric or infinite feature columns detected: ", paste(bad_cols, collapse = ", "), call. = FALSE)
  }
  x <- aml_impute_numeric(x)
  zero_var <- vapply(x, function(col) stats::sd(col) == 0, logical(1))
  if (any(zero_var)) {
    x <- x[, !zero_var, drop = FALSE]
  }
  if (ncol(x) < 2) {
    stop("At least 2 non-constant numeric features are required.", call. = FALSE)
  }

  list(
    x = x,
    y = y,
    y01 = as.integer(y == positive_class),
    positive_class = positive_class
  )
}

aml_impute_numeric <- function(x) {
  x[] <- lapply(x, function(col) {
    if (anyNA(col)) {
      med <- stats::median(col, na.rm = TRUE)
      if (!is.finite(med)) {
        med <- 0
      }
      col[is.na(col)] <- med
    }
    col
  })
  x
}

aml_rank_df <- function(method, feature, score, selected = NULL, extra = list()) {
  out <- data.frame(
    method = method,
    feature = as.character(feature),
    score = as.numeric(score),
    stringsAsFactors = FALSE
  )
  out <- out[is.finite(out$score), , drop = FALSE]
  out <- out[order(-out$score, out$feature), , drop = FALSE]
  out$rank <- seq_len(nrow(out))
  if (is.null(selected)) {
    selected <- out$feature
  }
  for (nm in names(extra)) {
    out[[nm]] <- extra[[nm]]
  }
  list(rankings = out, selected = unique(as.character(selected)))
}

aml_run_univariate <- function(x, y, top_n) {
  lev <- levels(y)
  scores <- vapply(x, function(col) {
    p <- tryCatch(stats::wilcox.test(col[y == lev[1]], col[y == lev[2]])$p.value, error = function(e) NA_real_)
    if (!is.finite(p) || p <= 0) {
      return(0)
    }
    -log10(p)
  }, numeric(1))
  ranked <- aml_rank_df("univariate", names(scores), scores)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked
}

aml_run_roc_auc <- function(x, y01, top_n) {
  aml_need_package("pROC")
  scores <- vapply(x, function(col) {
    auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(y01, col, quiet = TRUE, direction = "auto"))),
      error = function(e) NA_real_
    )
    if (!is.finite(auc)) {
      return(0)
    }
    abs(auc - 0.5) + 0.5
  }, numeric(1))
  ranked <- aml_rank_df("roc_auc", names(scores), scores)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked
}

aml_run_glmnet <- function(x, y, alpha, top_n, seed) {
  aml_need_package("glmnet")
  aml_need_package("Matrix")
  set.seed(seed)
  xx <- as.matrix(x)
  cv <- glmnet::cv.glmnet(
    x = xx,
    y = y,
    family = "binomial",
    alpha = alpha,
    nfolds = min(5, min(table(y))),
    standardize = TRUE
  )
  coefs <- as.matrix(stats::coef(cv, s = "lambda.min"))
  coefs <- coefs[rownames(coefs) != "(Intercept)", , drop = FALSE]
  score <- abs(as.numeric(coefs[, 1]))
  names(score) <- rownames(coefs)
  method <- if (alpha == 1) "lasso" else if (alpha == 0) "ridge" else "elastic_net"
  ranked <- aml_rank_df(method, names(score), score)
  selected <- ranked$rankings$feature[ranked$rankings$score > 0]
  if (length(selected) == 0) {
    selected <- head(ranked$rankings$feature, top_n)
  }
  ranked$selected <- head(selected, top_n)
  ranked
}

aml_run_random_forest <- function(x, y, top_n, seed) {
  aml_need_package("randomForest")
  set.seed(seed)
  fit <- randomForest::randomForest(x = x, y = y, ntree = 300, importance = TRUE)
  imp <- randomForest::importance(fit, type = 2)
  score <- as.numeric(imp[, 1])
  names(score) <- rownames(imp)
  ranked <- aml_rank_df("random_forest", names(score), score)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked
}

aml_run_boruta <- function(x, y, top_n, seed) {
  aml_need_package("Boruta")
  set.seed(seed)
  df <- data.frame(group = y, x, check.names = FALSE)
  fit <- Boruta::Boruta(group ~ ., data = df, maxRuns = 50, doTrace = 0)
  fixed <- Boruta::TentativeRoughFix(fit)
  stats <- Boruta::attStats(fixed)
  score <- stats$medianImp
  names(score) <- rownames(stats)
  selected <- rownames(stats)[stats$decision == "Confirmed"]
  if (length(selected) == 0) {
    selected <- head(names(sort(score, decreasing = TRUE)), top_n)
  }
  ranked <- aml_rank_df("boruta", names(score), score)
  ranked$rankings$decision <- stats[ranked$rankings$feature, "decision"]
  ranked$selected <- head(selected, top_n)
  ranked
}

aml_run_xgboost <- function(x, y01, top_n, seed) {
  aml_need_package("xgboost")
  set.seed(seed)
  dtrain <- xgboost::xgb.DMatrix(data = as.matrix(x), label = y01)
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 2,
    eta = 0.1,
    subsample = 0.9,
    colsample_bytree = 0.9,
    nthread = 1
  )
  fit <- xgboost::xgb.train(params = params, data = dtrain, nrounds = 40, verbose = 0)
  imp <- xgboost::xgb.importance(feature_names = colnames(x), model = fit)
  if (nrow(imp) == 0) {
    stop("xgboost returned no feature importance.", call. = FALSE)
  }
  ranked <- aml_rank_df("xgboost", imp$Feature, imp$Gain)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked
}

aml_need_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for this method.", call. = FALSE)
  }
}

aml_bind_rankings <- function(results) {
  frames <- lapply(results, function(z) z$rankings)
  all_cols <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(df) {
    missing_cols <- setdiff(all_cols, names(df))
    for (col in missing_cols) {
      df[[col]] <- NA
    }
    df[, all_cols, drop = FALSE]
  })
  do.call(rbind, frames)
}

aml_make_consensus <- function(results, n_methods) {
  selected <- lapply(names(results), function(method) {
    data.frame(method = method, feature = results[[method]]$selected, stringsAsFactors = FALSE)
  })
  selected <- do.call(rbind, selected)
  vote <- stats::aggregate(method ~ feature, selected, function(z) paste(sort(unique(z)), collapse = ";"))
  count <- stats::aggregate(method ~ feature, selected, function(z) length(unique(z)))
  names(vote)[2] <- "methods"
  names(count)[2] <- "n_methods"
  out <- merge(count, vote, by = "feature")
  out$support_rate <- out$n_methods / n_methods
  out <- out[order(-out$n_methods, out$feature), , drop = FALSE]
  rownames(out) <- NULL
  out
}

aml_write_selected_genes <- function(results, output_dir) {
  for (method in names(results)) {
    utils::write.table(
      data.frame(Genes = results[[method]]$selected),
      file = file.path(output_dir, paste0(method, "_selected_genes.txt")),
      row.names = FALSE,
      col.names = TRUE,
      quote = FALSE,
      sep = "\t"
    )
  }
}

aml_plot_consensus <- function(consensus, output_dir, top_n) {
  aml_need_package("ggplot2")
  plot_df <- head(consensus, top_n)
  plot_df$feature <- factor(plot_df$feature, levels = rev(plot_df$feature))
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = n_methods, y = feature)) +
    ggplot2::geom_col(fill = "#2C7FB8", width = 0.72) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::labs(x = "Number of supporting methods", y = NULL, title = "Consensus selected features")
  ggplot2::ggsave(file.path(output_dir, "consensus_feature_votes.pdf"), p, width = 8, height = max(5, 0.28 * nrow(plot_df) + 2))
  ggplot2::ggsave(file.path(output_dir, "consensus_feature_votes.png"), p, width = 8, height = max(5, 0.28 * nrow(plot_df) + 2), dpi = 300)
}

aml_plot_method_heatmap <- function(results, output_dir, top_n) {
  aml_need_package("ggplot2")
  consensus <- aml_make_consensus(results, length(results))
  features <- head(consensus$feature, top_n)
  heat <- do.call(rbind, lapply(names(results), function(method) {
    data.frame(
      method = method,
      feature = features,
      selected = features %in% results[[method]]$selected,
      stringsAsFactors = FALSE
    )
  }))
  heat$feature <- factor(heat$feature, levels = rev(features))
  heat$method <- factor(heat$method, levels = names(results))
  p <- ggplot2::ggplot(heat, ggplot2::aes(x = method, y = feature, fill = selected)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.35) +
    ggplot2::scale_fill_manual(values = c(`FALSE` = "grey92", `TRUE` = "#D95F02")) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "none") +
    ggplot2::labs(x = NULL, y = NULL, title = "Method support heatmap")
  ggplot2::ggsave(file.path(output_dir, "method_support_heatmap.pdf"), p, width = max(7, 1.1 * length(results)), height = max(5, 0.28 * length(features) + 2))
  ggplot2::ggsave(file.path(output_dir, "method_support_heatmap.png"), p, width = max(7, 1.1 * length(results)), height = max(5, 0.28 * length(features) + 2), dpi = 300)
}
