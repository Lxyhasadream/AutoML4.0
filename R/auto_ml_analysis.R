# AutoML4R: robust automated feature-selection workflow.

auto_ml_analysis <- function(hub_data,
                             group,
                             output_dir = "ML_screening_results",
                             positive_class = NULL,
                             methods = aml_methods(),
                             top_n = 20,
                             seed = 123,
                             write_outputs = TRUE,
                             verbose = TRUE) {
  set.seed(seed)
  methods <- match.arg(
    methods,
    choices = aml_methods(),
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
    stability_selection = function() aml_run_stability_selection(x, y01, top_n, seed),
    random_forest = function() aml_run_random_forest(x, y, top_n, seed),
    ranger = function() aml_run_ranger(x, y, top_n, seed),
    boruta = function() aml_run_boruta(x, y, top_n, seed),
    vsurf_interpretation = function() aml_run_vsurf(x, y, top_n, seed, stage = "interpretation"),
    vsurf_prediction = function() aml_run_vsurf(x, y, top_n, seed, stage = "prediction"),
    xgboost = function() aml_run_xgboost(x, y01, top_n, seed),
    lightgbm = function() aml_run_lightgbm(x, y01, top_n, seed),
    catboost = function() aml_run_catboost(x, y01, top_n, seed),
    mrmr = function() aml_run_mrmr(x, y01, top_n),
    splsda = function() aml_run_splsda(x, y, top_n, seed),
    pam = function() aml_run_pam(x, y, top_n, seed),
    svm = function() aml_run_svm(x, y, top_n, seed),
    svm_rfe = function() aml_run_svm_rfe(x, y, top_n, seed),
    relief = function() aml_run_relief_family(x, y, top_n, method = "ReliefFequalK"),
    information_gain = function() aml_run_fselector(x, y, top_n, method = "information_gain"),
    gain_ratio = function() aml_run_fselector(x, y, top_n, method = "gain_ratio"),
    symmetrical_uncertainty = function() aml_run_fselector(x, y, top_n, method = "symmetrical_uncertainty"),
    genetic_algorithm = function() aml_run_ga(x, y, top_n, seed),
    simulated_annealing = function() aml_run_sa(x, y, top_n, seed)
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
    aml_plot_algorithm_panels(results, output_dir, top_n, x = x, y = y, positive_class = positive_class)
    aml_plot_consensus(consensus, output_dir, top_n)
    aml_plot_method_heatmap(results, output_dir, top_n)
  }

  structure(
    list(
      selected_genes = lapply(results, function(z) z$selected),
      rankings = rankings,
      consensus = consensus,
      method_results = results,
      failures = failures,
      positive_class = positive_class,
      output_dir = if (write_outputs) normalizePath(output_dir, mustWork = FALSE) else NA_character_
    ),
    class = "auto_ml_analysis"
  )
}

aml_methods <- function() {
  c(
    "univariate",
    "roc_auc",
    "lasso",
    "ridge",
    "elastic_net",
    "stability_selection",
    "random_forest",
    "ranger",
    "boruta",
    "vsurf_interpretation",
    "vsurf_prediction",
    "xgboost",
    "lightgbm",
    "catboost",
    "mrmr",
    "splsda",
    "pam",
    "svm",
    "svm_rfe",
    "relief",
    "information_gain",
    "gain_ratio",
    "symmetrical_uncertainty",
    "genetic_algorithm",
    "simulated_annealing"
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
  ranked$details <- list(test = "wilcox")
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
  ranked$details <- list(metric = "absolute_auc")
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
  ranked$details <- list(cv = cv, alpha = alpha)
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
  ranked$details <- list(model = fit)
  ranked
}

aml_run_ranger <- function(x, y, top_n, seed) {
  aml_need_package("ranger")
  set.seed(seed)
  df <- data.frame(group = y, x, check.names = FALSE)
  fit <- ranger::ranger(
    group ~ .,
    data = df,
    num.trees = 300,
    importance = "permutation",
    probability = FALSE,
    seed = seed
  )
  score <- fit$variable.importance
  ranked <- aml_rank_df("ranger", names(score), score)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked$details <- list(model = fit)
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
  ranked$details <- list(model = fit, fixed_model = fixed, stats = stats)
  ranked
}

aml_run_vsurf <- function(x, y, top_n, seed, stage = c("interpretation", "prediction")) {
  aml_need_package("VSURF")
  stage <- match.arg(stage)
  set.seed(seed)
  fit <- VSURF::VSURF(
    x = x,
    y = y,
    ntree = 500,
    nfor.thres = 50,
    nfor.interp = 50,
    nfor.pred = 50,
    parallel = FALSE
  )
  idx <- if (stage == "interpretation") fit$varselect.interp else fit$varselect.pred
  selected <- colnames(x)[idx]
  if (length(selected) == 0) {
    selected <- head(colnames(x), top_n)
  }
  score <- rep(0, ncol(x))
  names(score) <- colnames(x)
  score[selected] <- rev(seq_along(selected))
  ranked <- aml_rank_df(paste0("vsurf_", stage), names(score), score)
  ranked$selected <- head(selected, top_n)
  ranked$details <- list(model = fit, stage = stage)
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
  ranked$details <- list(model = fit, importance = imp)
  ranked
}

aml_run_lightgbm <- function(x, y01, top_n, seed) {
  aml_need_package("lightgbm")
  set.seed(seed)
  dtrain <- lightgbm::lgb.Dataset(data = as.matrix(x), label = y01)
  params <- list(
    objective = "binary",
    metric = "binary_logloss",
    learning_rate = 0.05,
    num_leaves = 7,
    feature_fraction = 0.9,
    bagging_fraction = 0.9,
    bagging_freq = 1,
    num_threads = 1,
    verbosity = -1,
    seed = seed
  )
  fit <- lightgbm::lgb.train(params = params, data = dtrain, nrounds = 40, verbose = -1)
  imp <- lightgbm::lgb.importance(fit, percentage = TRUE)
  if (nrow(imp) == 0) {
    score <- aml_univariate_auc_score(x, factor(y01, levels = c(0, 1)))
    ranked <- aml_rank_df("lightgbm", names(score), score)
    ranked$selected <- head(ranked$rankings$feature, top_n)
    ranked$details <- list(model = fit, x = x, y01 = y01, fallback = TRUE)
    return(ranked)
  }
  ranked <- aml_rank_df("lightgbm", imp$Feature, imp$Gain)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked$details <- list(model = fit, importance = imp, x = x, y01 = y01, fallback = FALSE)
  ranked
}

aml_run_catboost <- function(x, y01, top_n, seed) {
  aml_need_package("catboost")
  set.seed(seed)
  pool <- catboost::catboost.load_pool(data = x, label = y01)
  fit <- catboost::catboost.train(
    learn_pool = pool,
    params = list(
      loss_function = "Logloss",
      iterations = 60,
      depth = 3,
      learning_rate = 0.05,
      random_seed = seed,
      logging_level = "Silent"
    )
  )
  score <- catboost::catboost.get_feature_importance(fit, pool = pool, type = "FeatureImportance")
  names(score) <- colnames(x)
  ranked <- aml_rank_df("catboost", names(score), score)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked$details <- list(model = fit, pool = pool)
  ranked
}

aml_run_mrmr <- function(x, y01, top_n) {
  aml_need_package("mRMRe")
  df <- data.frame(target = y01, x, check.names = FALSE)
  df[] <- lapply(df, as.numeric)
  mrmr_data <- mRMRe::mRMR.data(data = df)
  n_select <- min(top_n, ncol(x))
  fit <- mRMRe::mRMR.classic(data = mrmr_data, target_indices = 1, feature_count = n_select)
  idx <- mRMRe::solutions(fit)[[1]]
  names_all <- mRMRe::featureNames(mrmr_data)
  selected <- names_all[idx]
  selected <- selected[selected != "target"]
  score <- rep(0, ncol(x))
  names(score) <- colnames(x)
  score[selected] <- rev(seq_along(selected))
  ranked <- aml_rank_df("mrmr", names(score), score)
  ranked$selected <- head(selected, top_n)
  ranked$details <- list(model = fit)
  ranked
}

aml_run_splsda <- function(x, y, top_n, seed) {
  aml_need_package("mixOmics")
  set.seed(seed)
  keepx <- min(top_n, ncol(x))
  fit <- mixOmics::splsda(X = as.matrix(x), Y = y, ncomp = 1, keepX = keepx)
  loadings <- mixOmics::selectVar(fit, comp = 1)$value
  if (is.null(loadings) || nrow(loadings) == 0) {
    stop("sPLS-DA returned no selected variables.", call. = FALSE)
  }
  feature <- rownames(loadings)
  score <- abs(loadings[, 1])
  ranked <- aml_rank_df("splsda", feature, score)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked$details <- list(model = fit, loadings = loadings)
  ranked
}

aml_run_pam <- function(x, y, top_n, seed) {
  aml_need_package("pamr")
  set.seed(seed)
  dat <- list(
    x = t(as.matrix(x)),
    y = y,
    geneid = colnames(x),
    genenames = colnames(x)
  )
  fit <- pamr::pamr.train(dat)
  thresh <- fit$threshold[which.min(abs(fit$threshold - stats::median(fit$threshold)))]
  genes <- pamr::pamr.listgenes(fit, dat, threshold = thresh)
  if (is.null(genes) || nrow(genes) == 0) {
    selected <- head(colnames(x), top_n)
    score <- rep(0, ncol(x))
    names(score) <- colnames(x)
  } else {
    gene_col <- intersect(c("id", "geneid", "genenames", "name"), colnames(genes))[1]
    if (is.na(gene_col)) {
      gene_col <- 1
    }
    selected <- as.character(genes[, gene_col])
    score <- rep(0, ncol(x))
    names(score) <- colnames(x)
    selected <- intersect(selected, names(score))
    score[selected] <- rev(seq_along(selected))
  }
  ranked <- aml_rank_df("pam", names(score), score)
  ranked$selected <- head(selected, top_n)
  ranked$details <- list(model = fit, threshold = thresh, genes = genes)
  ranked
}

aml_run_svm <- function(x, y, top_n, seed) {
  aml_need_package("caret")
  aml_need_package("e1071")
  set.seed(seed)
  ctrl <- caret::trainControl(method = "cv", number = min(3, min(table(y))))
  fit <- caret::train(
    x = x,
    y = y,
    method = "svmLinear",
    trControl = ctrl,
    preProcess = c("center", "scale"),
    tuneLength = 3
  )
  imp <- caret::varImp(fit, scale = FALSE)$importance
  score <- rowMeans(imp, na.rm = TRUE)
  ranked <- aml_rank_df("svm", names(score), score)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked$details <- list(model = fit)
  ranked
}

aml_run_svm_rfe <- function(x, y, top_n, seed) {
  aml_need_package("caret")
  aml_need_package("e1071")
  set.seed(seed)
  sizes <- unique(pmax(1, pmin(ncol(x), c(1, 2, 5, 10, top_n))))
  ctrl <- caret::rfeControl(functions = caret::caretFuncs, method = "cv", number = min(3, min(table(y))))
  fit <- caret::rfe(
    x = x,
    y = y,
    sizes = sizes,
    rfeControl = ctrl,
    method = "svmLinear",
    preProcess = c("center", "scale"),
    tuneLength = 3
  )
  selected <- caret::predictors(fit)
  score <- rep(0, ncol(x))
  names(score) <- colnames(x)
  score[selected] <- rev(seq_along(selected))
  ranked <- aml_rank_df("svm_rfe", names(score), score)
  ranked$selected <- head(selected, top_n)
  ranked$details <- list(model = fit)
  ranked
}

aml_run_stability_selection <- function(x, y01, top_n, seed) {
  aml_need_package("glmnet")
  set.seed(seed)
  xx <- as.matrix(x)
  B <- 40
  counts <- numeric(ncol(x))
  names(counts) <- colnames(x)
  for (i in seq_len(B)) {
    idx <- unlist(lapply(split(seq_along(y01), y01), function(z) {
      sample(z, max(2, floor(length(z) / 2)))
    }))
    cv <- glmnet::cv.glmnet(
      x = xx[idx, , drop = FALSE],
      y = y01[idx],
      family = "binomial",
      alpha = 1,
      nfolds = min(3, min(table(y01[idx]))),
      standardize = TRUE
    )
    beta <- as.vector(stats::coef(cv, s = "lambda.min"))[-1]
    selected_idx <- which(abs(beta) > 0)
    if (length(selected_idx) > 0) {
      counts[selected_idx] <- counts[selected_idx] + 1
    }
  }
  score <- counts / B
  names(score) <- colnames(x)
  selected <- names(score)[score >= 0.6]
  if (length(selected) == 0) {
    selected <- names(sort(score, decreasing = TRUE))[seq_len(min(top_n, length(score)))]
  }
  ranked <- aml_rank_df("stability_selection", names(score), score)
  ranked$selected <- head(selected, top_n)
  ranked$details <- list(selection_frequency = score, n_resamples = B)
  ranked
}

aml_run_relief_family <- function(x, y, top_n, method = "ReliefFequalK") {
  aml_need_package("CORElearn")
  df <- data.frame(group = y, x, check.names = FALSE)
  form <- stats::as.formula("group ~ .")
  score <- CORElearn::attrEval(form, data = df, estimator = method)
  ranked <- aml_rank_df("relief", names(score), score)
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked$details <- list(estimator = method)
  ranked
}

aml_run_fselector <- function(x, y, top_n, method) {
  if (method == "information_gain") {
    aml_need_package("FSelectorRcpp")
    df <- data.frame(group = y, x, check.names = FALSE)
    form <- stats::as.formula("group ~ .")
    scores <- FSelectorRcpp::information_gain(form, df)
    score_col <- setdiff(colnames(scores), "attributes")[1]
    ranked <- aml_rank_df(method, scores$attributes, scores[[score_col]])
  } else {
    score <- vapply(x, function(col) aml_information_score(col, y, method), numeric(1))
    ranked <- aml_rank_df(method, names(score), score)
  }
  ranked$selected <- head(ranked$rankings$feature, top_n)
  ranked$details <- list(metric = method)
  ranked
}

aml_run_ga <- function(x, y, top_n, seed) {
  set.seed(seed)
  base_score <- aml_univariate_auc_score(x, y)
  selected <- aml_stochastic_subset_search(x, y, top_n, seed, iterations = 40, temperature = 0, base_score)
  score <- rep(0, ncol(x))
  names(score) <- colnames(x)
  score[selected] <- rev(seq_along(selected))
  ranked <- aml_rank_df("genetic_algorithm", names(score), score)
  ranked$selected <- head(selected, top_n)
  ranked$details <- list(search = "stochastic_genetic_style", base_score = base_score)
  ranked
}

aml_run_sa <- function(x, y, top_n, seed) {
  set.seed(seed)
  base_score <- aml_univariate_auc_score(x, y)
  selected <- aml_stochastic_subset_search(x, y, top_n, seed, iterations = 60, temperature = 0.08, base_score)
  score <- rep(0, ncol(x))
  names(score) <- colnames(x)
  score[selected] <- rev(seq_along(selected))
  ranked <- aml_rank_df("simulated_annealing", names(score), score)
  ranked$selected <- head(selected, top_n)
  ranked$details <- list(search = "simulated_annealing_style", base_score = base_score)
  ranked
}

aml_discretize <- function(col, bins = 5) {
  probs <- seq(0, 1, length.out = bins + 1)
  cuts <- unique(stats::quantile(col, probs = probs, na.rm = TRUE, type = 7))
  if (length(cuts) < 3) {
    return(factor(rep("single", length(col))))
  }
  cut(col, breaks = cuts, include.lowest = TRUE, ordered_result = TRUE)
}

aml_entropy <- function(z) {
  tab <- table(z)
  p <- as.numeric(tab) / sum(tab)
  p <- p[p > 0]
  -sum(p * log2(p))
}

aml_information_score <- function(col, y, method) {
  xdisc <- aml_discretize(col)
  hy <- aml_entropy(y)
  hx <- aml_entropy(xdisc)
  joint <- interaction(xdisc, y, drop = TRUE)
  hxy <- aml_entropy(joint)
  ig <- hx + hy - hxy
  if (!is.finite(ig)) {
    return(0)
  }
  if (method == "gain_ratio") {
    if (!is.finite(hx) || hx <= 0) {
      return(0)
    }
    return(ig / hx)
  }
  if (method == "symmetrical_uncertainty") {
    denom <- hx + hy
    if (!is.finite(denom) || denom <= 0) {
      return(0)
    }
    return(2 * ig / denom)
  }
  ig
}

aml_univariate_auc_score <- function(x, y) {
  y01 <- as.integer(y == levels(y)[2])
  vapply(x, function(col) {
    auc <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(y01, col, quiet = TRUE, direction = "auto"))),
      error = function(e) NA_real_
    )
    if (!is.finite(auc)) {
      return(0)
    }
    abs(auc - 0.5) + 0.5
  }, numeric(1))
}

aml_subset_score <- function(features, base_score) {
  if (length(features) == 0) {
    return(0)
  }
  mean(base_score[features], na.rm = TRUE) + 0.01 * sqrt(length(features))
}

aml_stochastic_subset_search <- function(x, y, top_n, seed, iterations, temperature, base_score) {
  set.seed(seed)
  features <- names(sort(base_score, decreasing = TRUE))
  current <- head(features, min(top_n, length(features)))
  current_score <- aml_subset_score(current, base_score)
  best <- current
  best_score <- current_score

  for (i in seq_len(iterations)) {
    candidate <- current
    if (runif(1) < 0.5 && length(candidate) > 1) {
      candidate <- setdiff(candidate, sample(candidate, 1))
    }
    pool <- setdiff(features, candidate)
    if (length(pool) > 0 && length(candidate) < top_n) {
      weights <- base_score[pool] + 1e-6
      candidate <- unique(c(candidate, sample(pool, 1, prob = weights)))
    }
    if (length(candidate) == 0) {
      candidate <- sample(features, 1)
    }

    candidate_score <- aml_subset_score(candidate, base_score)
    accept <- candidate_score >= current_score
    if (!accept && temperature > 0) {
      accept <- runif(1) < exp((candidate_score - current_score) / temperature)
    }
    if (accept) {
      current <- candidate
      current_score <- candidate_score
    }
    if (candidate_score > best_score) {
      best <- candidate
      best_score <- candidate_score
    }
  }

  remaining <- setdiff(features, best)
  head(c(best, remaining), top_n)
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
