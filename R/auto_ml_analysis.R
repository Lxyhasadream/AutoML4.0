
# Auto Machine Learning Script - Optimized by Xinyu Li on 2025-01-20
# Auto Machine Learning Script - Updated by Xinyu Li on March 4, 2025
# YUsheng Li  
# Ensure you have hub_data and group files prepared in advance.
# hub_data: a matrix of preliminary selected features.
# group: a file containing grouping information.
#again update at 2025-05-26

auto_ml_analysis <- function(hub_data, group, output_dir = "ML_screening_results") {
  runtime_pkgs <- c(
    "Boruta", "ggplot2", "tibble", "ggsci", "glmnet", "randomForest",
    "mlbench", "caret", "xgboost", "Matrix", "PRROC", "shapviz", "dplyr",
    "tidyverse", "lightgbm", "tidyr", "ggrepel", "catboost", "mRMRe",
    "mixOmics", "VSURF", "CORElearn", "ranger", "pROC", "stabs", "pamr",
    "FSelectorRcpp", "stringr", "scales", "patchwork"
  )
  missing_pkgs <- runtime_pkgs[!vapply(runtime_pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required runtime packages: ",
      paste(missing_pkgs, collapse = ", "),
      ". Install them before running auto_ml_analysis().",
      call. = FALSE
    )
  }
 # Load required libraries
  library(Boruta)
  library(ggplot2)
  library(tibble)
  library(ggsci)
  library(glmnet)
  library(randomForest)
  library(mlbench)
  library(caret)
  library(xgboost)
  library(Matrix)
  library(PRROC)
  library(shapviz)
  library(mlbench)
  library(caret)
  library('xgboost')  
  library("Matrix")
  library('PRROC')
  library('ggplot2')
  library('dplyr')
  library("tidyverse")
  
  hub_data <- as.data.frame(hub_data, check.names = FALSE)
  group <- group
# hub_data <- as.data.frame(t(hub_data))
  mycolors <- c('#E64A35','#4DBBD4' ,'#01A187'  ,'#6BD66B','#3C5588'  ,'#F29F80'  ,
                  '#8491B6','#91D0C1','#7F5F48','#AF9E85','#4F4FFF',
                  '#739B57','#EFE685','#446983')
  group <- as.factor(group)
  # Check output directory
  if (missing(output_dir) || output_dir == "") {
    stop("Output directory must be specified.")
  }
  set.seed(123)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  # Validate inputs
  if (missing(hub_data) || missing(group)) {
    stop("hub_data and group must be provided as data frames.")
  }
  # Initialize a list to store selected genes
  all_genes <- list()
  iter.times=1000
  all_genes <- list()
  

#Boruta Model------
  cat("Running Boruta model...\n")
  
  Var.Selec<-Boruta(
    group~.,
    hub_data,
    pValue = 0.01, #confidence level. 
    mcAdj = TRUE, #Bonferroni
    #if set to TRUE, a multiple comparisons adjustment using the Bonferroni method will be applied. 
    maxRuns = 100, #
    doTrace = 0,#0-30
    holdHistory = TRUE, #TRUEImpHistory
    getImp = getImpRfZ # getImpRfZ ranger  Z 
  )

  library(Boruta)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  
  # Output directory
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Fix tentative features
  
  Var.Selec.final <- TentativeRoughFix(Var.Selec)
  
  # Extract Boruta results
  
  Boruta_gene <- getSelectedAttributes(
    Var.Selec.final,
    withTentative = FALSE
  )
  
  Boruta_stats <- attStats(Var.Selec.final) %>%
    rownames_to_column("Gene") %>%
    arrange(desc(medianImp))
  
  Boruta_selected_stats <- Boruta_stats %>%
    filter(decision == "Confirmed") %>%
    arrange(desc(medianImp))
  
  # Export selected genes
  
  write.table(
    Boruta_gene,
    file = file.path(output_dir, "Boruta_Genes.txt"),
    row.names = FALSE,
    col.names = "Genes",
    quote = FALSE,
    sep = "\t"
  )
  
  cat("Boruta Selected Genes:\n", paste(Boruta_gene, collapse = ", "), "\n\n")
  
  # Export full Boruta statistics
  
  write.table(
    Boruta_stats,
    file = file.path(output_dir, "Boruta_All_Feature_Statistics.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  # Boruta importance history plot
  
  pdf(
    file = file.path(output_dir, "Model_Boruta_ImpHistory_optimized.pdf"),
    width = 12,
    height = 8
  )
  
  par(
    mar = c(7, 6, 4, 2),
    cex.lab = 1.5,
    cex.axis = 1.2,
    cex.main = 1.6,
    font.lab = 2
  )
  
  plotImpHistory(
    Var.Selec.final,
    whichShadow = c(TRUE, TRUE, TRUE),
    ylab = "Z-scores",
    xlab = "Boruta iterations",
    las = 1,
    main = "Boruta Importance History"
  )
  
  dev.off()
  
  # Boruta all-feature boxplot
  
  n_features <- length(Var.Selec.final$finalDecision)
  full_width <- max(14, min(80, n_features * 0.22))
  
  pdf(
    file = file.path(output_dir, "Model_Boruta_All_Features_optimized.pdf"),
    width = full_width,
    height = 9
  )
  
  par(
    mar = c(12, 6, 4, 2),
    cex.lab = 1.6,
    cex.axis = 1.0,
    cex.main = 1.6,
    font.lab = 2
  )
  
  plot(
    Var.Selec.final,
    whichShadow = c(FALSE, FALSE, FALSE),
    xlab = "",
    ylab = "Z-scores",
    las = 2,
    cex.axis = 0.9,
    main = "Boruta Feature Importance"
  )
  
  dev.off()
  
  # Boruta selected-feature boxplot
  
  if (length(Boruta_gene) > 0) {
    
    selected_width <- max(10, min(30, length(Boruta_gene) * 0.45))
    
    pdf(
      file = file.path(output_dir, "Model_Boruta_Selected_Features_optimized.pdf"),
      width = selected_width,
      height = 8
    )
    
    par(
      mar = c(12, 6, 4, 2),
      cex.lab = 1.6,
      cex.axis = 1.1,
      cex.main = 1.6,
      font.lab = 2
    )
    
    plot(
      Var.Selec.final,
      whichShadow = c(FALSE, FALSE, FALSE),
      xlab = "",
      ylab = "Z-scores",
      las = 2,
      cex.axis = 1.0,
      main = "Boruta-Confirmed Features"
    )
    
    dev.off()
  }
  
  # Boruta selected-feature barplot
  
  if (nrow(Boruta_selected_stats) > 0) {
    
    p_boruta <- ggplot(
      Boruta_selected_stats,
      aes(x = reorder(Gene, medianImp), y = medianImp)
    ) +
      geom_col(
        width = 0.75,
        fill = "#D55E00"
      ) +
      coord_flip() +
      theme_classic(base_size = 16) +
      labs(
        x = NULL,
        y = "Median Boruta importance",
        title = "Boruta-confirmed feature importance"
      ) +
      theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          color = "black"
        ),
        axis.text.y = element_text(
          size = 13,
          color = "black"
        ),
        axis.text.x = element_text(
          size = 13,
          color = "black"
        ),
        axis.title.x = element_text(
          size = 16,
          face = "bold",
          color = "black"
        ),
        axis.line = element_line(
          linewidth = 0.7,
          color = "black"
        ),
        axis.ticks = element_line(
          linewidth = 0.7,
          color = "black"
        ),
        plot.margin = ggplot2::margin(
          10, 20, 10, 10,
          unit = "pt"
        )
      )
    
    ggsave(
      filename = file.path(output_dir, "Model_Boruta_Selected_Features_Barplot.pdf"),
      plot = p_boruta,
      width = 10,
      height = max(6, nrow(Boruta_selected_stats) * 0.35),
      limitsize = FALSE
    )
    
    ggsave(
      filename = file.path(output_dir, "Model_Boruta_Selected_Features_Barplot.png"),
      plot = p_boruta,
      width = 10,
      height = max(6, nrow(Boruta_selected_stats) * 0.35),
      dpi = 600,
      limitsize = FALSE
    )
  }
  
  # Compatibility version for old ggplot2
  
  if (nrow(Boruta_selected_stats) > 0) {
    
    p_boruta_old <- ggplot(
      Boruta_selected_stats,
      aes(x = reorder(Gene, medianImp), y = medianImp)
    ) +
      geom_col(
        width = 0.75,
        fill = "#D55E00"
      ) +
      coord_flip() +
      theme_classic(base_size = 16) +
      labs(
        x = NULL,
        y = "Median Boruta importance",
        title = "Boruta-confirmed feature importance"
      ) +
      theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          color = "black"
        ),
        axis.text.y = element_text(
          size = 13,
          color = "black"
        ),
        axis.text.x = element_text(
          size = 13,
          color = "black"
        ),
        axis.title.x = element_text(
          size = 16,
          face = "bold",
          color = "black"
        ),
        axis.line = element_line(
          size = 0.7,
          color = "black"
        ),
        axis.ticks = element_line(
          size = 0.7,
          color = "black"
        ),
        plot.margin = ggplot2::margin(
          10, 20, 10, 10,
          unit = "pt"
        )
      )
  }
  
  
  
  
#Ridge-----
  cat("Running Ridge model...\n")
  x <- as.matrix(hub_data)
  y <- group
  ridge <- glmnet(x, y, family = 'binomial', nlambda = 1000, alpha = 0)
  cvfit <- cv.glmnet(x, y, nfolds = 10, family = "binomial", type.measure = "deviance")
  
  get
  cat(paste0("#lambda","\n\n",cvfit$lambda.min))
  cat(paste0("lambda.1se","\n\n",cvfit$lambda.1se))
  coefficients<-coef(cvfit,s=cvfit$lambda.min)#
  Active.Index<-which(coefficients!=0)
  Active.coefficients<- coefficients[Active.Index]
  coefficients
  Ridge_genes<- colnames(x)[(Active.Index-1)[-1]]
  Ridge_genes 
  all_genes$Ridge_genes <- Ridge_genes
  
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  
  # Output directory
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Prepare data
  
  cat("Running Ridge model...\n")
  
  x <- as.matrix(hub_data)
  y <- group
  
  if (is.character(y)) {
    y <- as.factor(y)
  }
  
  if (is.factor(y)) {
    y <- droplevels(y)
  }
  
  if (length(unique(y)) != 2) {
    stop("Ridge logistic regression requires exactly two groups.")
  }
  
  if (!exists("all_genes") || !is.list(all_genes)) {
    all_genes <- list()
  }
  
  set.seed(123)
  
  # Fit Ridge model with cross-validation
  
  cvfit <- cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 0,
    nfolds = 10,
    type.measure = "deviance",
    nlambda = 1000,
    standardize = TRUE
  )
  
  ridge_fit <- cvfit$glmnet.fit
  
  cat("lambda.min:\n", cvfit$lambda.min, "\n\n")
  cat("lambda.1se:\n", cvfit$lambda.1se, "\n\n")
  
  # Extract Ridge coefficients
  
  coef_min <- coef(cvfit, s = "lambda.min")
  coef_mat <- as.matrix(coef_min)
  
  Ridge_coef_df <- data.frame(
    Gene = rownames(coef_mat),
    Coefficient = as.numeric(coef_mat[, 1]),
    stringsAsFactors = FALSE
  ) %>%
    filter(Gene != "(Intercept)") %>%
    mutate(
      AbsCoefficient = abs(Coefficient)
    ) %>%
    arrange(desc(AbsCoefficient))
  
  # Select top Ridge-ranked genes
  
  top_n <- min(20, nrow(Ridge_coef_df))
  
  Ridge_selected_df <- Ridge_coef_df %>%
    slice_head(n = top_n)
  
  Ridge_genes <- Ridge_selected_df$Gene
  all_genes$Ridge_genes <- Ridge_genes
  
  cat("Top Ridge-ranked genes:\n", paste(Ridge_genes, collapse = ", "), "\n\n")
  
  # Export Ridge results
  
  write.table(
    Ridge_genes,
    file = file.path(output_dir, "Ridge_Top_Ranked_Genes.txt"),
    row.names = FALSE,
    col.names = "Genes",
    quote = FALSE,
    sep = "\t"
  )
  
  write.table(
    Ridge_coef_df,
    file = file.path(output_dir, "Ridge_All_Coefficients.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  write.table(
    Ridge_selected_df,
    file = file.path(output_dir, "Ridge_Top_Ranked_Coefficients.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  # Ridge cross-validation curve
  
  pdf(
    file = file.path(output_dir, "Model_Ridge_CV_Curve.pdf"),
    width = 8,
    height = 7
  )
  
  par(
    mar = c(6, 6, 4, 2),
    cex.lab = 1.5,
    cex.axis = 1.2,
    cex.main = 1.5,
    font.lab = 2
  )
  
  plot(cvfit)
  
  title(
    main = "Ridge Cross-Validation Curve",
    line = 2.5,
    cex.main = 1.5,
    font.main = 2
  )
  
  dev.off()
  
  # Ridge coefficient path
  
  pdf(
    file = file.path(output_dir, "Model_Ridge_Coefficient_Path.pdf"),
    width = 9,
    height = 7
  )
  
  par(
    mar = c(6, 6, 4, 2),
    cex.lab = 1.5,
    cex.axis = 1.2,
    cex.main = 1.5,
    font.lab = 2
  )
  
  plot(
    ridge_fit,
    xvar = "lambda",
    label = FALSE,
    xlab = "Log(lambda)",
    ylab = "Coefficients"
  )
  
  abline(
    v = log(cvfit$lambda.min),
    lty = 2,
    lwd = 2
  )
  
  abline(
    v = log(cvfit$lambda.1se),
    lty = 3,
    lwd = 2
  )
  
  title(
    main = "Ridge Coefficient Path",
    line = 2.5,
    cex.main = 1.5,
    font.main = 2
  )
  
  legend(
    "topright",
    legend = c("lambda.min", "lambda.1se"),
    lty = c(2, 3),
    lwd = 2,
    bty = "n",
    cex = 1.1
  )
  
  dev.off()
  
  # Ridge absolute coefficient ranking plot
  
  p_ridge_abs <- ggplot(
    Ridge_selected_df,
    aes(x = reorder(Gene, AbsCoefficient), y = AbsCoefficient)
  ) +
    geom_col(
      width = 0.75,
      fill = "#0072B2"
    ) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "Absolute Ridge coefficient",
      title = paste0("Top ", top_n, " Ridge-ranked genes")
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_Ridge_Top_Ranked_Genes_AbsCoefficient.pdf"),
    plot = p_ridge_abs,
    width = 10,
    height = max(6, top_n * 0.35),
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_Ridge_Top_Ranked_Genes_AbsCoefficient.png"),
    plot = p_ridge_abs,
    width = 10,
    height = max(6, top_n * 0.35),
    dpi = 600,
    limitsize = FALSE
  )
  
  # Ridge signed coefficient ranking plot
  
  p_ridge_signed <- ggplot(
    Ridge_selected_df,
    aes(x = reorder(Gene, Coefficient), y = Coefficient)
  ) +
    geom_col(
      width = 0.75,
      fill = "#0072B2"
    ) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "Ridge coefficient",
      title = paste0("Top ", top_n, " Ridge-ranked genes")
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_Ridge_Top_Ranked_Genes_SignedCoefficient.pdf"),
    plot = p_ridge_signed,
    width = 10,
    height = max(6, top_n * 0.35),
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_Ridge_Top_Ranked_Genes_SignedCoefficient.png"),
    plot = p_ridge_signed,
    width = 10,
    height = max(6, top_n * 0.35),
    dpi = 600,
    limitsize = FALSE
  )
  
  
#LASSO Model,1000------------
  cat("Running LASSO model...\n")
  x <- as.matrix(hub_data)
  y <- group
  lasso_fea_list <- list()
  list.of.seed <- 1:1000
  y <- ifelse(group=='Ctrl',0,1)#Ctrl0
  ###
  for (i in seq_along(list.of.seed)) {
    seed <- list.of.seed[i]  # 
    set.seed(seed)  # 
    
    #  cv.glmnet 
    cvfit <- cv.glmnet(
      x = as.matrix(x),  #  x 
      y = y,             # 
      nfolds = 10,       # 10-fold 
      alpha = 1,         # Lasso 
      family = "binomial", # 
      maxit = 1000       # 
    )
    
    # 
    coefs <- coef(cvfit, s = "lambda.min")  #  lambda 
    fea <- rownames(coefs)[coefs[, 1] != 0]  # 
    
    # 
    if ("(Intercept)" %in% fea) {
      lasso_fea <- sort(fea[fea != "(Intercept)"])  # 
    } else {
      lasso_fea <- sort(fea)
    }
    
    # 
    lasso_fea_list[[i]] <- lasso_fea
  }
  
  # 
  lasso_res <- NULL
  for (i in 1:1000) {
    lasso_res <- rbind.data.frame(lasso_res,
                                  data.frame(
                                    iteration = i,
                                    n.gene = length(lasso_fea_list[[i]]),
                                    genelist = paste0(lasso_fea_list[[i]], collapse = " | "),
                                    stringsAsFactors = F
                                  ),
                                  stringsAsFactors = F
    )
  }
 ##
  genes <- sort(table(unlist(lasso_fea_list)), decreasing = T) # 
  freq.cutoff <- 1000 * 0.05
  genes <- names(genes[genes > freq.cutoff]) # 50lasso. 95%
  
  
  result <- data.frame(
    method = c(rep("Lasso", length(genes))),
    selected.fea = genes
  )

  

  

##Lasoo------
  lasso <- glmnet(x, y, family = 'binomial', nlambda = 1000, alpha = 1)
  cvfit <- cv.glmnet(x, y, nfolds = 10, family = "binomial", type.measure = "deviance")
  
  #get
  cat(paste0("#lambda","\n\n",cvfit$lambda.min))
  cat(paste0("lambda.1se","\n\n",cvfit$lambda.1se))
  coefficients<-coef(cvfit,s=cvfit$lambda.min)#
  Active.Index<-which(coefficients!=0)
  Active.coefficients<- coefficients[Active.Index]
  coefficients
  lasso_genes<- colnames(x)[(Active.Index-1)[-1]]
  lasso_genes 
  # Save LASSO genes and print
  write.table(lasso_genes, file = file.path(output_dir, "LASSO_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$LASSO <- lasso_genes
  cat("LASSO Selected Genes:\n", paste(lasso_genes, collapse = ", "), "\n\n")
  
  #
  
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Output directory
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (!exists("all_genes") || !is.list(all_genes)) {
    all_genes <- list()
  }
  
  # Prepare LASSO input data
  
  cat("Running LASSO model...\n")
  
  x <- as.matrix(hub_data)
  y <- group
  
  if (is.character(y)) {
    y <- as.factor(y)
  }
  
  if (is.factor(y)) {
    y <- droplevels(y)
  }
  
  if (length(unique(y)) != 2) {
    stop("LASSO logistic regression requires exactly two groups.")
  }
  
  set.seed(123)
  
  # Fit LASSO model
  
  cvfit_lasso <- cv.glmnet(
    x = x,
    y = y,
    family = "binomial",
    alpha = 1,
    nfolds = 10,
    type.measure = "deviance",
    nlambda = 1000,
    standardize = TRUE
  )
  
  fit_lasso <- cvfit_lasso$glmnet.fit
  
  cat("lambda.min:\n", cvfit_lasso$lambda.min, "\n\n")
  cat("lambda.1se:\n", cvfit_lasso$lambda.1se, "\n\n")
  
  # Extract LASSO coefficients at lambda.min
  
  coef_min <- coef(cvfit_lasso, s = "lambda.min")
  coef_mat_min <- as.matrix(coef_min)
  
  lasso_coef_df <- data.frame(
    Gene = rownames(coef_mat_min),
    Coefficient = as.numeric(coef_mat_min[, 1]),
    stringsAsFactors = FALSE
  ) %>%
    filter(Gene != "(Intercept)") %>%
    mutate(
      AbsCoefficient = abs(Coefficient)
    ) %>%
    arrange(desc(AbsCoefficient))
  
  lasso_selected_df <- lasso_coef_df %>%
    filter(Coefficient != 0) %>%
    arrange(desc(AbsCoefficient))
  
  lasso_genes <- lasso_selected_df$Gene
  
  all_genes$LASSO <- lasso_genes
  
  cat("LASSO Selected Genes:\n", paste(lasso_genes, collapse = ", "), "\n\n")
  
  # Export LASSO results
  
  write.table(
    lasso_genes,
    file = file.path(output_dir, "LASSO_Genes.txt"),
    row.names = FALSE,
    col.names = "Genes",
    quote = FALSE,
    sep = "\t"
  )
  
  write.table(
    lasso_coef_df,
    file = file.path(output_dir, "LASSO_All_Coefficients.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  write.table(
    lasso_selected_df,
    file = file.path(output_dir, "LASSO_Selected_Coefficients.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  # Prepare coefficient path data from the same LASSO model
  
  coef_path_mat <- as.matrix(fit_lasso$beta)
  
  coef_path_df <- as.data.frame(
    t(coef_path_mat),
    check.names = FALSE
  )
  
  coef_path_df$lambda <- fit_lasso$lambda
  coef_path_df$log_lambda <- log(fit_lasso$lambda)
  
  coef_path_long <- coef_path_df %>%
    pivot_longer(
      cols = -c(lambda, log_lambda),
      names_to = "Gene",
      values_to = "Coefficient"
    ) %>%
    mutate(
      Selected = ifelse(Gene %in% lasso_genes, "Selected", "Other")
    )
  
  coef_path_selected <- coef_path_long %>%
    filter(Gene %in% lasso_genes)
  
  # Check lambda consistency
  
  cat("Coefficient-path log(lambda) range:\n")
  print(range(coef_path_long$log_lambda, na.rm = TRUE))
  
  cat("log(lambda.min):\n")
  print(log(cvfit_lasso$lambda.min))
  
  cat("log(lambda.1se):\n")
  print(log(cvfit_lasso$lambda.1se))
  
  # LASSO coefficient path
  
  p_lasso_path <- ggplot() +
    geom_line(
      data = coef_path_long,
      aes(x = log_lambda, y = Coefficient, group = Gene),
      color = "grey75",
      linewidth = 0.4,
      alpha = 0.8
    ) +
    geom_line(
      data = coef_path_selected,
      aes(x = log_lambda, y = Coefficient, color = Gene, group = Gene),
      linewidth = 1.0
    ) +
    geom_vline(
      xintercept = log(cvfit_lasso$lambda.min),
      linetype = "dashed",
      linewidth = 0.8,
      color = "#E41A1C"
    ) +
    geom_vline(
      xintercept = log(cvfit_lasso$lambda.1se),
      linetype = "dashed",
      linewidth = 0.8,
      color = "#377EB8"
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = "Log(lambda)",
      y = "Coefficients",
      title = "LASSO coefficient path"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        size = 13,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 11,
        color = "black"
      ),
      legend.position = ifelse(length(lasso_genes) <= 15, "right", "none"),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_LASSO_Coefficient_Path_optimized.pdf"),
    plot = p_lasso_path,
    width = 10,
    height = 8,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_LASSO_Coefficient_Path_optimized.png"),
    plot = p_lasso_path,
    width = 10,
    height = 8,
    dpi = 600,
    limitsize = FALSE
  )
  
  # LASSO coefficient path without legend
  
  p_lasso_path_nolegend <- p_lasso_path +
    theme(
      legend.position = "none"
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_LASSO_Coefficient_Path_NoLegend.pdf"),
    plot = p_lasso_path_nolegend,
    width = 9,
    height = 7,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_LASSO_Coefficient_Path_NoLegend.png"),
    plot = p_lasso_path_nolegend,
    width = 9,
    height = 7,
    dpi = 600,
    limitsize = FALSE
  )
  
  # Prepare cross-validation data
  
  cv_df <- data.frame(
    lambda = cvfit_lasso$lambda,
    log_lambda = log(cvfit_lasso$lambda),
    cvm = cvfit_lasso$cvm,
    cvsd = cvfit_lasso$cvsd,
    nzero = cvfit_lasso$nzero
  ) %>%
    mutate(
      ymin = cvm - cvsd,
      ymax = cvm + cvsd
    )
  
  cv_y_max <- max(cv_df$ymax, na.rm = TRUE)
  cv_y_min <- min(cv_df$ymin, na.rm = TRUE)
  cv_y_range <- cv_y_max - cv_y_min
  
  cv_label_y1 <- cv_y_max - 0.08 * cv_y_range
  cv_label_y2 <- cv_y_max - 0.18 * cv_y_range
  
  y_axis_label <- ifelse(
    !is.null(cvfit_lasso$name),
    cvfit_lasso$name,
    "Cross-validation error"
  )
  
  # LASSO cross-validation curve
  
  p_lasso_cv <- ggplot(
    cv_df,
    aes(x = log_lambda, y = cvm)
  ) +
    geom_errorbar(
      aes(ymin = ymin, ymax = ymax),
      width = 0.05,
      linewidth = 0.5,
      color = "grey50"
    ) +
    geom_point(
      size = 2.2,
      color = "black"
    ) +
    geom_vline(
      xintercept = log(cvfit_lasso$lambda.min),
      linetype = "dashed",
      linewidth = 0.9,
      color = "#E41A1C"
    ) +
    geom_vline(
      xintercept = log(cvfit_lasso$lambda.1se),
      linetype = "dashed",
      linewidth = 0.9,
      color = "#377EB8"
    ) +
    annotate(
      geom = "text",
      x = log(cvfit_lasso$lambda.min),
      y = cv_label_y1,
      label = paste0("lambda.min\n", signif(cvfit_lasso$lambda.min, 4)),
      color = "#E41A1C",
      size = 4,
      hjust = -0.05
    ) +
    annotate(
      geom = "text",
      x = log(cvfit_lasso$lambda.1se),
      y = cv_label_y2,
      label = paste0("lambda.1se\n", signif(cvfit_lasso$lambda.1se, 4)),
      color = "#377EB8",
      size = 4,
      hjust = -0.05
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = "Log(lambda)",
      y = y_axis_label,
      title = "10-fold cross-validation for LASSO"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        size = 13,
        color = "black"
      ),
      axis.line = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        linewidth = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_LASSO_CV_Curve_optimized.pdf"),
    plot = p_lasso_cv,
    width = 8,
    height = 7,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_LASSO_CV_Curve_optimized.png"),
    plot = p_lasso_cv,
    width = 8,
    height = 7,
    dpi = 600,
    limitsize = FALSE
  )
  
  # LASSO selected-gene coefficient barplot
  
  if (nrow(lasso_selected_df) > 0) {
    
    p_lasso_bar <- ggplot(
      lasso_selected_df,
      aes(x = reorder(Gene, AbsCoefficient), y = AbsCoefficient)
    ) +
      geom_col(
        width = 0.75,
        fill = "#D55E00"
      ) +
      coord_flip() +
      theme_classic(base_size = 16) +
      labs(
        x = NULL,
        y = "Absolute LASSO coefficient",
        title = "LASSO-selected genes"
      ) +
      theme(
        plot.title = element_text(
          size = 18,
          face = "bold",
          hjust = 0.5,
          color = "black"
        ),
        axis.text.y = element_text(
          size = 13,
          color = "black"
        ),
        axis.text.x = element_text(
          size = 13,
          color = "black"
        ),
        axis.title.x = element_text(
          size = 16,
          face = "bold",
          color = "black"
        ),
        axis.line = element_line(
          linewidth = 0.7,
          color = "black"
        ),
        axis.ticks = element_line(
          linewidth = 0.7,
          color = "black"
        ),
        plot.margin = ggplot2::margin(
          10, 20, 10, 10,
          unit = "pt"
        )
      )
    
    ggsave(
      filename = file.path(output_dir, "Model_LASSO_Selected_Genes_Barplot.pdf"),
      plot = p_lasso_bar,
      width = 10,
      height = max(6, nrow(lasso_selected_df) * 0.35),
      limitsize = FALSE
    )
    
    ggsave(
      filename = file.path(output_dir, "Model_LASSO_Selected_Genes_Barplot.png"),
      plot = p_lasso_bar,
      width = 10,
      height = max(6, nrow(lasso_selected_df) * 0.35),
      dpi = 600,
      limitsize = FALSE
    )
  }
  
  
#Enet------
  x <- as.matrix(hub_data)
  y <- group
  alpha_fea_list <- list()
  alpha_fea_list_all <-list() 
  for (alpha in seq(0.1, 0.9, 0.1)) {
    set.seed(seed)  # 
    #  alpha 
    #  list.of.seed
    for (i in 1:1000) {
      seed <- list.of.seed[i]  # 
      set.seed(seed)  # 
      
      #  cv.glmnet 
      cvfit <- cv.glmnet(
        x = as.matrix(x),  #  x1 
        y = y,             # 
        nfolds = 10,        # 10-fold 
        alpha = alpha,      #  alpha 
        family = "binomial", # 
        maxit = 1000        # 
      )
      
      # 
      coefs <- coef(cvfit, s = "lambda.min")  #  lambda 
      fea <- rownames(coefs)[coefs[, 1] != 0]  # 
      
      # 
      if ("(Intercept)" %in% fea) {
        lasso_fea <- sort(fea[fea != "(Intercept)"])  # 
      } else {
        lasso_fea <- sort(fea)
      }
      
      #  alpha 
      alpha_fea_list[[i]] <- lasso_fea
    }
    
    #  alpha  fea_list 
    alpha_fea_list_all[[as.character(alpha)]] <- alpha_fea_list
    
  }
  

alpha_fea_list_all[["0.2"]]
alpha_fea_list_all[["0.3"]]
alpha_fea_list_all[["0.4"]]  
alpha_fea_list_all[["0.5"]]  
tmp <- c()
all <- c()
final <- list()  # 
for (i in c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) {
  tmp <- c()
  all <- c()
  for (j in 1:1000) {
    tmp <- alpha_fea_list_all[[i]][[j]]
    all <- c(all, tmp)
  }
  
  freq_table <- as.data.frame(table(all))
  freq_table <- as.data.frame(table(all)) %>%
    arrange(desc(Freq))  #  Freq 
  final[[i]] <- freq_table  #  list 
}

for (i in c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) {
  all_genes[[paste0("Enet_", i)]] <- final[[i]]  #  all_genes
  cat("all_genes$", paste0("Enet_", i), ":\n", sep = "")
  print(all_genes[[paste0("Enet_", i)]]$all)  # 
  cat("\n")  # 
}
#
all_genes$Enet_0.1 <- final$"0.1"#0.1
all_genes$Enet_0.2 <- final$"0.2"#0.2
all_genes$Enet_0.3 <- final$"0.3"#0.3
all_genes$Enet_0.4 <- final$"0.4"#0.4
all_genes$Enet_0.5 <- final$"0.5"#0.5
all_genes$Enet_0.6 <- final$"0.6"#0.6
all_genes$Enet_0.7 <- final$"0.7"#0.7
all_genes$Enet_0.8 <- final$"0.8"#0.8
all_genes$Enet_0.9 <- final$"0.9"#0.9


library(ggplot2)
#  0.1  0.9  9 
library(ggplot2)

for (i in c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9")) {
  # 
  df <- all_genes[[paste0("Enet_", i)]]
  # 
  df$Freq <- df$Freq / 1000
  # 
  p <- ggplot(df, aes(x = reorder(all, Freq), y = Freq, fill = Freq)) +
    geom_bar(stat = "identity", width = 0.7, color = "black") +  # 
    scale_fill_gradient(low = "lightblue", high = "#F4A7B9") +  # 
    labs(
      title = paste("Top Features (Enet_", i, ")", sep = ""),  # 
      x = "Feature",  # x
      y = "Relative Frequency"  # y
    ) +
    theme_minimal() +  # 
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  # 
      axis.title.x = element_text(size = 14, face = "bold"),  # x
      axis.title.y = element_text(size = 14, face = "bold"),  # y
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),  # x
      axis.text.y = element_text(size = 12)  # y
    ) +
    coord_flip()  # 
  
  # 
  ggsave(filename = paste0("Enet_", i, "_barplot.pdf"), plot = p, width = 8, height = 6, dpi = 300)
  
  # 
  print(p)
}



#Random Forest Model------
  cat("Running Random Forest model...\n")
  RF <- randomForest(group ~ ., data = hub_data, ntree = 1000, importance = TRUE)
  eror <- importance(RF) %>% as.data.frame()
  imp <- data.frame(ID = rownames(eror), imp = eror$MeanDecreaseAccuracy)
  imp$relative_imp <- (imp$imp - min(imp$imp)) / (max(imp$imp) - min(imp$imp))
  RF_genes <- imp$ID[imp$relative_imp > 0.5]
  rf_data<- imp[order(imp$relative_imp,decreasing = T),c(1,3)]
  write.table(RF_genes, file = file.path(output_dir, "RF_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$RandomForest <- RF_genes
  cat("Random Forest Selected Genes:\n", paste(RF_genes, collapse = ", "), "\n\n")
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Output directory
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (!exists("all_genes") || !is.list(all_genes)) {
    all_genes <- list()
  }
  
  # Prepare RF importance data
  
  importance_cutoff <- 0.5
  
  rf_data_plot <- rf_data %>%
    mutate(
      ID = as.character(ID),
      relative_imp = as.numeric(relative_imp)
    ) %>%
    arrange(relative_imp) %>%
    mutate(
      ID = factor(ID, levels = ID)
    )
  
  RF_genes <- rf_data_plot %>%
    filter(relative_imp >= importance_cutoff) %>%
    pull(ID) %>%
    as.character()
  
  all_genes$RF <- RF_genes
  
  cat("RF Selected Genes:\n", paste(RF_genes, collapse = ", "), "\n\n")
  
  write.table(
    RF_genes,
    file = file.path(output_dir, "RF_Genes.txt"),
    row.names = FALSE,
    col.names = "Genes",
    quote = FALSE,
    sep = "\t"
  )
  
  write.table(
    rf_data_plot,
    file = file.path(output_dir, "RF_Relative_Importance.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  # RF relative importance plot
  
  rf_height <- max(6, nrow(rf_data_plot) * 0.32)
  
  p_rf_imp <- ggplot(
    rf_data_plot,
    aes(x = ID, y = relative_imp)
  ) +
    geom_col(
      width = 0.65,
      fill = "#1E90FF"
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    geom_hline(
      yintercept = importance_cutoff,
      linetype = "dashed",
      color = "red",
      size = 0.7
    ) +
    coord_flip() +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "Relative importance",
      title = "Random forest feature importance"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 12,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_RF_Feature_Importance.pdf"),
    plot = p_rf_imp,
    width = 10,
    height = rf_height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_RF_Feature_Importance.png"),
    plot = p_rf_imp,
    width = 10,
    height = rf_height,
    dpi = 600,
    limitsize = FALSE
  )
  
  # Prepare RF OOB error data
  
  oob_col <- if ("OOB" %in% colnames(RF$err.rate)) {
    "OOB"
  } else {
    colnames(RF$err.rate)[1]
  }
  
  tmp <- data.frame(
    trees = seq_len(nrow(RF$err.rate)),
    error = RF$err.rate[, oob_col],
    stringsAsFactors = FALSE
  )
  
  best_tree <- tmp$trees[which.min(tmp$error)]
  best_error <- min(tmp$error, na.rm = TRUE)
  
  cat("Best number of trees:\n", best_tree, "\n\n")
  cat("Minimum OOB error:\n", best_error, "\n\n")
  
  write.table(
    tmp,
    file = file.path(output_dir, "RF_OOB_Error.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
  )
  
  # RF OOB error curve
  
  p_rf_oob <- ggplot(
    tmp,
    aes(x = trees, y = error)
  ) +
    geom_line(
      size = 0.7,
      color = "gray35"
    ) +
    geom_vline(
      xintercept = best_tree,
      linetype = "dashed",
      color = "red",
      size = 0.8
    ) +
    geom_point(
      data = tmp %>% filter(trees == best_tree),
      aes(x = trees, y = error),
      color = "red",
      size = 3
    ) +
    annotate(
      geom = "text",
      x = best_tree,
      y = best_error,
      label = paste0("Best tree = ", best_tree),
      color = "red",
      size = 4,
      hjust = -0.05,
      vjust = -0.8
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = "Number of trees",
      y = "OOB error",
      title = "Random forest OOB error curve"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text = element_text(
        size = 13,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_RF_OOB_Error_Curve.pdf"),
    plot = p_rf_oob,
    width = 9,
    height = 7,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_RF_OOB_Error_Curve.png"),
    plot = p_rf_oob,
    width = 9,
    height = 7,
    dpi = 600,
    limitsize = FALSE
  )
  
  # Optional: RF full error-rate curve
  
  rf_error_df <- as.data.frame(RF$err.rate) %>%
    mutate(
      trees = seq_len(nrow(RF$err.rate))
    ) %>%
    pivot_longer(
      cols = -trees,
      names_to = "Error_type",
      values_to = "Error"
    )
  
  p_rf_error_all <- ggplot(
    rf_error_df,
    aes(x = trees, y = Error, color = Error_type)
  ) +
    geom_line(
      size = 0.7
    ) +
    geom_vline(
      xintercept = best_tree,
      linetype = "dashed",
      color = "red",
      size = 0.8
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = "Number of trees",
      y = "Error rate",
      title = "Random forest error-rate curve"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text = element_text(
        size = 13,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 12,
        color = "black"
      ),
      legend.position = "top",
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_RF_All_Error_Rates.pdf"),
    plot = p_rf_error_all,
    width = 9,
    height = 7,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_RF_All_Error_Rates.png"),
    plot = p_rf_error_all,
    width = 9,
    height = 7,
    dpi = 600,
    limitsize = FALSE
  )
  
  
  
  
#svm Model----
  cat("Svm model...\n")
  group
  control <- rfeControl(functions = caretFuncs,method = "cv", number = 10)  
  
  # SVM-RFE
  set.seed(123)
  n <- dim(hub_data)[1]
  results <- rfe(hub_data,            #
                 as.factor(group),   #
                 sizes = c(n:1), 
                 rfeControl = control,
                 method = "svmLinear")
  
  #rfe5
  #https://github.com/cran/caret/tree/master/R 
  top = 5 #
  cat("The top ",
      min(top, results$bestSubset),
      " variables (out of ",
      results$bestSubset,
      "):\n   ",
      paste(results$optVariables[1:min(top, results$bestSubset)], collapse = ", "),
      "\n\n",
      sep = "")

  
  write.table(Svm_genes, file = file.path(output_dir, "Svm_Genes.txt"), row.names = FALSE, col.names = "Genes")
  all_genes$Svm <- Svm_genes
  cat("Svm Selected Genes:\n", paste(Svm_genes, collapse = ", "), "\n\n") 
  

  library(ggplot2)
  library(dplyr)
  
  # Prepare SVM-RFE plotting data
  
  svm_rfe_df <- results$results
  
  if (!"Variables" %in% colnames(svm_rfe_df)) {
    stop("The column 'Variables' was not found in results$results.")
  }
  
  metric_name <- if ("Accuracy" %in% colnames(svm_rfe_df)) {
    "Accuracy"
  } else if ("Kappa" %in% colnames(svm_rfe_df)) {
    "Kappa"
  } else if (!is.null(results$metric) && results$metric %in% colnames(svm_rfe_df)) {
    results$metric
  } else {
    stop("No valid performance metric was found in results$results.")
  }
  
  metric_sd_name <- paste0(metric_name, "SD")
  
  svm_rfe_df <- svm_rfe_df %>%
    mutate(
      Variables = as.numeric(Variables),
      Performance = .data[[metric_name]],
      PerformanceSD = if (metric_sd_name %in% colnames(.)) {
        .data[[metric_sd_name]]
      } else {
        NA_real_
      }
    ) %>%
    arrange(Variables)
  
  best_subset <- results$bestSubset
  
  best_perf <- svm_rfe_df %>%
    filter(Variables == best_subset) %>%
    pull(Performance)
  
  if (length(best_perf) == 0) {
    best_perf <- max(svm_rfe_df$Performance, na.rm = TRUE)
  }
  
  # Optimized SVM-RFE performance curve
  
  p_svm_rfe <- ggplot(
    svm_rfe_df,
    aes(x = Variables, y = Performance)
  ) +
    geom_line(
      color = "grey35",
      size = 0.8
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    {
      if (!all(is.na(svm_rfe_df$PerformanceSD))) {
        geom_errorbar(
          aes(
            ymin = Performance - PerformanceSD,
            ymax = Performance + PerformanceSD
          ),
          width = 0.15,
          color = "grey55",
          size = 0.5
        )
      }
    } +
    geom_vline(
      xintercept = best_subset,
      linetype = "dashed",
      color = "#E41A1C",
      size = 0.8
    ) +
    geom_point(
      data = svm_rfe_df %>% filter(Variables == best_subset),
      aes(x = Variables, y = Performance),
      color = "#E41A1C",
      size = 4
    ) +
    annotate(
      geom = "text",
      x = best_subset,
      y = best_perf,
      label = paste0("Optimal features = ", best_subset),
      color = "#E41A1C",
      size = 4.2,
      hjust = -0.05,
      vjust = -0.8
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = "Number of selected features",
      y = metric_name,
      title = "SVM-RFE cross-validation performance"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.text = element_text(
        size = 13,
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 25, 10, 10,
        unit = "pt"
      )
    ) +
    scale_x_continuous(
      breaks = pretty(svm_rfe_df$Variables, n = 8)
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_SVM_RFE_Performance_optimized.pdf"),
    plot = p_svm_rfe,
    width = 9,
    height = 7,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_SVM_RFE_Performance_optimized.png"),
    plot = p_svm_rfe,
    width = 9,
    height = 7,
    dpi = 600,
    limitsize = FALSE
  )
  
  # Selected SVM-RFE gene ranking plot
  
  top <- 5
  
  Svm_genes <- results$optVariables[
    1:min(top, results$bestSubset)
  ]
  
  svm_gene_rank_df <- data.frame(
    Gene = Svm_genes,
    Rank = seq_along(Svm_genes),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Gene = factor(Gene, levels = rev(Gene))
    )
  
  p_svm_genes <- ggplot(
    svm_gene_rank_df,
    aes(x = Gene, y = Rank)
  ) +
    geom_col(
      width = 0.65,
      fill = "#1E90FF"
    ) +
    geom_point(
      color = "#B22222",
      size = 3
    ) +
    coord_flip() +
    scale_y_reverse(
      breaks = svm_gene_rank_df$Rank
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "SVM-RFE rank",
      title = paste0("Top ", length(Svm_genes), " SVM-RFE-selected genes")
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_SVM_RFE_Selected_Genes_optimized.pdf"),
    plot = p_svm_genes,
    width = 8,
    height = max(5, length(Svm_genes) * 0.45),
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_SVM_RFE_Selected_Genes_optimized.png"),
    plot = p_svm_genes,
    width = 8,
    height = max(5, length(Svm_genes) * 0.45),
    dpi = 600,
    limitsize = FALSE
  ) 
  
  
  
  
  
  
#XGBoost Model-----

  cat("Running XGBoost model...\n")
  library(xgboost)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  
  set.seed(123)
  
  xgboost_hub_data <- hub_data

  positive_class <- "CAA"
  label <- ifelse(as.character(group) == positive_class, 1, 0)
  label <- as.numeric(label)
  cat("XGBoost labels:\n")
  print(table(label))
  
  xgboost_hub_data <- as.data.frame(xgboost_hub_data, check.names = FALSE)
  
  xgboost_hub_data[] <- lapply(xgboost_hub_data, function(x) {
    if (is.factor(x)) {
      as.numeric(as.character(x))
    } else {
      as.numeric(x)
    }
  })
  #  sparse matrix
  train_matrix <- Matrix::Matrix(
    as.matrix(xgboost_hub_data),
    sparse = TRUE
  )
  
  colnames(train_matrix) <- colnames(xgboost_hub_data)
  
  #  DMatrix
  dtrain <- xgb.DMatrix(
    data = train_matrix,
    label = label
  )
  
  params <- list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 5,
    learning_rate = 0.5
  )
  
  res.xgb <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 25,
    verbose = 1
  )
  cat("XGBoost model finished.\n")
  xgb_importance <- xgb.importance(
    feature_names = colnames(train_matrix),
    model = res.xgb
  )
  xgb_genes <- xgb_importance$Feature
  write.table(
    data.frame(Genes = xgb_genes),
    file = file.path(output_dir, "XGBoost_Selected_Genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  cat("XGBoost Selected Genes:\n")
  cat(paste(xgb_genes, collapse = ", "), "\n\n")
  xgboost_number <- nrow(xgb_importance)
  cat("XGBoost feature importance plot saved.\n")
  all_genes$XGBoost <- xgb_genes

  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  
  # Prepare XGBoost plotting data
  
  if (!all(c("Feature", "Gain") %in% colnames(xgb_importance))) {
    stop("xgb_importance must contain at least two columns: Feature and Gain.")
  }
  
  xgb_top_n <- min(20, nrow(xgb_importance))
  
  xgb_importance_plot <- xgb_importance %>%
    as.data.frame() %>%
    mutate(
      Feature = as.character(Feature),
      Gain = as.numeric(Gain)
    ) %>%
    arrange(desc(Gain)) %>%
    slice_head(n = xgb_top_n) %>%
    mutate(
      Feature = factor(Feature, levels = rev(Feature))
    )
  
  xgb_plot_height <- max(6, nrow(xgb_importance_plot) * 0.35)
  
  # XGBoost Gain importance plot
  
  p_xgb_gain <- ggplot(
    xgb_importance_plot,
    aes(x = Feature, y = Gain)
  ) +
    geom_col(
      width = 0.72,
      fill = "#1E90FF",
      color = "black",
      size = 0.25
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "Gain",
      title = "XGBoost feature importance"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_XGBoost_Gain_Importance_optimized.pdf"),
    plot = p_xgb_gain,
    width = 10,
    height = xgb_plot_height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_XGBoost_Gain_Importance_optimized.png"),
    plot = p_xgb_gain,
    width = 10,
    height = xgb_plot_height,
    dpi = 600,
    limitsize = FALSE
  )
  
  
#lightgbm+shap------------
cat("Running LightGBM + SHAP model...\n")
library(lightgbm)
library(dplyr)
library(ggplot2)
library(tidyr)
set.seed(123)
lightgbm_hub_data <- hub_data
#  data.frame
lightgbm_hub_data <- as.data.frame(lightgbm_hub_data, check.names = FALSE)
#  numeric
lightgbm_hub_data[] <- lapply(lightgbm_hub_data, function(x) {
    if (is.factor(x)) {
      as.numeric(as.character(x))
    } else {
      as.numeric(x)
    }
  })

positive_class <- "CAA"
label <- ifelse(as.character(group) == positive_class, 1, 0)
label <- as.numeric(label)
cat("LightGBM labels:\n")
print(table(label))
train_matrix <- as.matrix(lightgbm_hub_data)
storage.mode(train_matrix) <- "numeric"
feature_names <- colnames(train_matrix)
lgb_train <- lgb.Dataset(
    data = train_matrix,
    label = label,
    free_raw_data = FALSE
  )
params <- list(
    objective = "binary",
    metric = "binary_logloss",
    boosting = "gbdt",
    learning_rate = 0.03,
    num_leaves = 4,
    max_depth = 2,
    min_data_in_leaf = 3,
    feature_fraction = 0.8,
    bagging_fraction = 0.8,
    bagging_freq = 1,
    lambda_l1 = 0.1,
    lambda_l2 = 0.1,
    verbosity = -1
  )
  
lgb_model <- lgb.train(
    params = params,
    data = lgb_train,
    nrounds = 200,
    verbose = -1
  )
  
cat("LightGBM model finished.\n")
lgb_importance <- lgb.importance(
    model = lgb_model,
    percentage = TRUE
  )

write.csv(
    lgb_importance,
    file = file.path(output_dir, "LightGBM_feature_importance.csv"),
    row.names = FALSE
  )
  
cat("LightGBM feature importance:\n")
print(lgb_importance)
###shap
shap_values <- predict(
  object = lgb_model,
  newdata = train_matrix,
  type = "contrib"
)

shap_values <- as.data.frame(shap_values, check.names = FALSE)

# LightGBM type = "contrib" usually returns:
# n_features columns + 1 final base-value column for binary classification.
# The final base-value column is not a gene and should be removed.

if (ncol(shap_values) == length(feature_names) + 1) {
  
  shap_values_no_bias <- shap_values[, 1:length(feature_names), drop = FALSE]
  
} else {
  
  # Fallback: remove BIAS / base column if it has an identifiable name
  bias_col <- grep(
    "BIAS|base|intercept",
    colnames(shap_values),
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(bias_col) > 0) {
    shap_values_no_bias <- shap_values[, !colnames(shap_values) %in% bias_col, drop = FALSE]
  } else {
    shap_values_no_bias <- shap_values
  }
}



colnames(shap_values_no_bias) <- feature_names

# Calculate mean absolute SHAP value for each gene
shap_importance <- data.frame(
  gene = colnames(shap_values_no_bias),
  mean_abs_SHAP = colMeans(abs(shap_values_no_bias), na.rm = TRUE),
  stringsAsFactors = FALSE
) %>%
  dplyr::arrange(desc(mean_abs_SHAP))

write.csv(
  shap_importance,
  file = file.path(output_dir, "LightGBM_SHAP_importance.csv"),
  row.names = FALSE
)

cat("LightGBM SHAP importance:\n")
print(shap_importance)


lgb_gain_df <- lgb_importance %>%
    dplyr::select(Feature, Gain) %>%
    dplyr::rename(
      gene = Feature,
      Gain_importance = Gain
    )
  
combined_lgb_shap <- shap_importance %>%
    left_join(lgb_gain_df, by = "gene") %>%
    mutate(
      Gain_importance = ifelse(is.na(Gain_importance), 0, Gain_importance),
      SHAP_rank = rank(-mean_abs_SHAP, ties.method = "min"),
      Gain_rank = rank(-Gain_importance, ties.method = "min"),
      mean_rank = rowMeans(cbind(SHAP_rank, Gain_rank), na.rm = TRUE)
    ) %>%
    arrange(mean_rank)
  
write.csv(
    combined_lgb_shap,
    file = file.path(output_dir, "LightGBM_SHAP_combined_importance.csv"),
    row.names = FALSE
  )
  
# 8. Select LightGBM + SHAP genes

  top_n <- min(20, nrow(shap_importance))
  
  lgb_shap_genes_topN <- shap_importance$gene[1:top_n]
  
  #  mean_abs_SHAP 
  lgb_shap_genes_above_mean <- shap_importance %>%
    filter(mean_abs_SHAP > mean(mean_abs_SHAP, na.rm = TRUE)) %>%
    pull(gene)
  
  # SHAP top20  Gain top20 
  gain_topN <- lgb_importance$Feature[1:min(20, nrow(lgb_importance))]
  
  lgb_shap_intersect_genes <- intersect(
    lgb_shap_genes_topN,
    gain_topN
  )
  
  #  SHAP topN
  lgb_shap_selected_genes <- lgb_shap_genes_topN
  
  write.table(
    data.frame(Genes = lgb_shap_selected_genes),
    file = file.path(output_dir, "LightGBM_SHAP_Selected_Genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  
  write.table(
    data.frame(Genes = lgb_shap_intersect_genes),
    file = file.path(output_dir, "LightGBM_SHAP_Gain_Intersect_Genes.txt"),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  
cat("LightGBM + SHAP selected genes:\n")
cat(paste(lgb_shap_selected_genes, collapse = ", "), "\n\n")
cat("LightGBM SHAP-Gain intersect genes:\n")
cat(paste(lgb_shap_intersect_genes, collapse = ", "), "\n\n")
all_genes$LightGBM_SHAP = lgb_shap_selected_genes
lgb_plot_df <- lgb_importance %>%
    arrange(Gain) %>%
    mutate(
      Feature = factor(Feature, levels = Feature)
    )
p_lgb_gain <- ggplot(
    lgb_plot_df,
    aes(x = Feature, y = Gain)
  ) +
    geom_col(
      width = 0.62,
      fill = "#2E86DE",
      color = "black",
      linewidth = 0.35
    ) +
    geom_point(
      shape = 21,
      fill = "#F39C12",
      color = "black",
      size = 3.1,
      stroke = 0.4
    ) +
    coord_flip() +
    labs(
      x = "Gene",
      y = "Gain importance",
      title = "LightGBM Feature Importance"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 13, face = "bold"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.line = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  p_lgb_gain
  
pdf(
    file = file.path(output_dir, "Model_LightGBM_Gain.pdf"),
    width = 7,
    height = 5.5
  )
print(p_lgb_gain)
dev.off()
  
#SHAP
shap_plot_df <- shap_importance %>%
    arrange(mean_abs_SHAP) %>%
    mutate(
      gene = factor(gene, levels = gene),
      selected = ifelse(gene %in% lgb_shap_selected_genes, "Selected", "Other")
    )
  
p_lgb_shap <- ggplot(
    shap_plot_df,
    aes(x = gene, y = mean_abs_SHAP, fill = selected)
  ) +
    geom_col(
      width = 0.62,
      color = "black",
      linewidth = 0.35
    ) +
    geom_point(
      shape = 21,
      fill = "#F39C12",
      color = "black",
      size = 3.1,
      stroke = 0.4
    ) +
    scale_fill_manual(
      values = c(
        "Selected" = "#D55E00",
        "Other" = "grey75"
      )
    ) +
    coord_flip() +
    labs(
      x = "Gene",
      y = "Mean absolute SHAP value",
      fill = "",
      title = "LightGBM SHAP Feature Importance"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 13, face = "bold"),
      legend.position = "right",
      legend.text = element_text(color = "black", size = 11),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.line = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  p_lgb_shap
  
  pdf(
    file = file.path(output_dir, "Model_LightGBM_SHAP.pdf"),
    width = 7,
    height = 5.5
  )
  print(p_lgb_shap)
  dev.off()
  
  ggsave(
    filename = file.path(output_dir, "Model_LightGBM_SHAP.png"),
    plot = p_lgb_shap,
    width = 7,
    height = 5.5,
    dpi = 600
  )
  

  
combined_plot_df <- combined_lgb_shap %>%
    arrange(mean_rank) %>%
    mutate(
      gene = factor(gene, levels = rev(gene)),
      selected = ifelse(gene %in% lgb_shap_selected_genes, "Selected", "Other")
    )
  
  p_lgb_combined <- ggplot(
    combined_plot_df,
    aes(x = gene, y = mean_rank, fill = selected)
  ) +
    geom_col(
      width = 0.62,
      color = "black",
      linewidth = 0.35
    ) +
    coord_flip() +
    scale_y_reverse() +
    scale_fill_manual(
      values = c(
        "Selected" = "#D55E00",
        "Other" = "grey75"
      )
    ) +
    labs(
      x = "Gene",
      y = "Combined rank",
      fill = "",
      title = "Integrated LightGBM-SHAP Gene Ranking"
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),
      axis.text = element_text(color = "black", size = 11),
      axis.title = element_text(color = "black", size = 13, face = "bold"),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
      axis.line = element_blank(),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  p_lgb_combined
  
pdf(
    file = file.path(output_dir, "Model_LightGBM_SHAP_Combined_Rank.pdf"),
    width = 7,
    height = 5.5
  )
print(p_lgb_combined)
dev.off()
  
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)


shap_mat <- as.data.frame(shap_values_no_bias, check.names = FALSE)
expr_mat <- as.data.frame(lightgbm_hub_data, check.names = FALSE)

stopifnot(all(colnames(shap_mat) == colnames(expr_mat)))
stopifnot(nrow(shap_mat) == nrow(expr_mat))

sample_group <- ifelse(label == 1, "CAA", "HC")

# Select top genes based on mean absolute SHAP
top_n <- min(20, nrow(shap_importance))
top_genes <- shap_importance$gene[1:top_n]

# Long-format SHAP data
shap_long <- shap_mat %>%
  dplyr::select(dplyr::all_of(top_genes)) %>%
  dplyr::mutate(sample_id = rownames(shap_mat)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(top_genes),
    names_to = "gene",
    values_to = "SHAP_value"
  )

expr_long <- expr_mat %>%
  dplyr::select(dplyr::all_of(top_genes)) %>%
  dplyr::mutate(sample_id = rownames(expr_mat)) %>%
  tidyr::pivot_longer(
    cols = dplyr::all_of(top_genes),
    names_to = "gene",
    values_to = "expression"
  )

shap_plot_df <- shap_long %>%
  dplyr::left_join(expr_long, by = c("sample_id", "gene")) %>%
  dplyr::mutate(
    group = rep(sample_group, times = length(top_genes)),
    gene = factor(gene, levels = rev(top_genes))
  )

# Scale expression within each gene for color mapping
shap_plot_df <- shap_plot_df %>%
  dplyr::group_by(gene) %>%
  dplyr::mutate(
    expression_scaled = as.numeric(scale(expression))
  ) %>%
  dplyr::ungroup()

shap_plot_df$expression_scaled[is.na(shap_plot_df$expression_scaled)] <- 0
  
p_shap_beeswarm <- ggplot(
  shap_plot_df,
  aes(
    x = SHAP_value,
    y = gene,
    color = expression_scaled
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.5
  ) +
  geom_jitter(
    height = 0.22,
    width = 0,
    size = 2.4,
    alpha = 0.85
  ) +
  scale_color_gradient2(
    low = "#2C7BB6",
    mid = "grey90",
    high = "#D7191C",
    midpoint = 0,
    name = "Scaled\nexpression"
  ) +
  labs(
    x = "SHAP value",
    y = "Gene",
    title = "LightGBM SHAP summary plot",
    subtitle = "Each point represents one sample; color indicates scaled gene expression"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(color = "black", size = 12, hjust = 0.5),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    legend.title = element_text(color = "black", size = 11, face = "bold"),
    legend.text = element_text(color = "black", size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.line = element_blank()
  )

p_shap_beeswarm

ggsave(
  file.path(output_dir, "LightGBM_SHAP_summary_beeswarm.pdf"),
  p_shap_beeswarm,
  width = 7.5,
  height = 6
)

shap_bar_df <- shap_importance %>%
  dplyr::slice_head(n = top_n) %>%
  dplyr::arrange(mean_abs_SHAP) %>%
  dplyr::mutate(
    gene = factor(gene, levels = gene)
  )

p_shap_bar <- ggplot(
  shap_bar_df,
  aes(x = gene, y = mean_abs_SHAP)
) +
  geom_col(
    width = 0.65,
    fill = "#2E86DE",
    color = "black",
    linewidth = 0.35
  ) +
  geom_point(
    shape = 21,
    fill = "#F39C12",
    color = "black",
    size = 3,
    stroke = 0.4
  ) +
  coord_flip() +
  labs(
    x = "Gene",
    y = "Mean absolute SHAP value",
    title = "Global SHAP importance"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(color = "black", size = 16, face = "bold", hjust = 0.5),
    axis.text = element_text(color = "black", size = 11),
    axis.title = element_text(color = "black", size = 14, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.7),
    axis.line = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank()
  )

p_shap_bar

ggsave(
  file.path(output_dir, "LightGBM_SHAP_summary_bar.pdf"),
  p_shap_bar,
  width = 7,
  height = 5.8
)


  
  
  
  
  
  


#CatBoost----
library(catboost)
labels <- ifelse(group == "CAA",1,0)
train <- as.data.frame(x)
train_pool <- catboost.load_pool(data = train,label = labels)
train_pool
model <- catboost.train(train_pool,  NULL,
                        params = list(loss_function = 'Logloss', # 
                                      iterations = 100, # 100
                                      metric_period=10 # 101
                                      #prediction_type=c("Class","Probability")
                        ) 
)

importance <- catboost.get_feature_importance(model,type='FeatureImportance')
importance <- data.frame(importance)
importance <- arrange(importance, desc(importance))
importance <- importance[importance$importance>mean(importance$importance),,drop = FALSE]
importance$gene <-rownames(importance) 
CatBoost_genes <- rownames(importance)
all_genes$CatBoost <- rownames(importance)
cat("Catboost Selected Genes:\n", paste(CatBoost_genes, collapse = ", "), "\n\n")

library(ggplot2)
library(dplyr)

# Prepare CatBoost plotting data

if (!"importance" %in% colnames(importance)) {
  stop("The column 'importance' was not found in the CatBoost importance table.")
}

if (!"gene" %in% colnames(importance)) {
  importance$gene <- rownames(importance)
}

importance_plot <- importance %>%
  as.data.frame() %>%
  mutate(
    gene = as.character(gene),
    importance = as.numeric(importance)
  ) %>%
  arrange(desc(importance)) %>%
  mutate(
    gene = factor(gene, levels = rev(gene))
  )

catboost_plot_height <- max(6, nrow(importance_plot) * 0.35)

# CatBoost feature importance plot

p_catboost <- ggplot(
  importance_plot,
  aes(x = gene, y = importance)
) +
  geom_col(
    width = 0.72,
    fill = "#1E90FF",
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Importance",
    title = "CatBoost feature importance"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_CatBoost_Feature_Importance_optimized.pdf"),
  plot = p_catboost,
  width = 10,
  height = catboost_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_CatBoost_Feature_Importance_optimized.png"),
  plot = p_catboost,
  width = 10,
  height = catboost_plot_height,
  dpi = 600,
  limitsize = FALSE
)











#----
positive_class <- "CAA"
x_df <- as.data.frame(hub_data, check.names = FALSE)
x_df[] <- lapply(x_df, function(x) {
  if (is.factor(x)) {
    as.numeric(as.character(x))
  } else {
    as.numeric(x)
  }
})
y_factor <- factor(group)
y01 <- ifelse(as.character(group) == positive_class, 1, 0)
y01 <- as.numeric(y01)
cat("Group distribution:\n")
print(table(y_factor))
add_gene_result <- function(all_genes, method_name, genes, output_dir) {
  genes <- unique(na.omit(as.character(genes)))
  all_genes[[method_name]] <- genes
  
  write.table(
    data.frame(Genes = genes),
    file = file.path(output_dir, paste0(method_name, "_Genes.txt")),
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    sep = "\t"
  )
  
  cat(method_name, "selected genes:\n")
  cat(paste(genes, collapse = ", "), "\n\n")
  
  return(all_genes)
}

#mRMR--------
library(mRMRe)
mrmr_df <- data.frame(
  target = y01,
  x_df,
  check.names = FALSE
)

mrmr_df[] <- lapply(mrmr_df, function(z) as.numeric(z))

mrmr_data <- mRMRe::mRMR.data(data = mrmr_df)

feature_count <- min(20, ncol(x_df))

mrmr_model <- mRMRe::mRMR.classic(
  data = mrmr_data,
  target_indices = 1,
  feature_count = feature_count
)

mrmr_indices <- mRMRe::solutions(mrmr_model)[[1]]
mrmr_feature_names <- mRMRe::featureNames(mrmr_data)

mrmr_genes <- mrmr_feature_names[mrmr_indices]
mrmr_genes <- setdiff(mrmr_genes, "target")

write.table(
  data.frame(Genes = mrmr_genes),
  file = file.path(output_dir, "mRMR_Genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

all_genes <- add_gene_result(
  all_genes,
  "mRMR",
  mrmr_genes,
  output_dir
)



library(ggplot2)
library(dplyr)

# Prepare mRMR selection-order plot

mrmr_rank_df <- data.frame(
  Gene = mrmr_genes,
  Rank = seq_along(mrmr_genes),
  stringsAsFactors = FALSE
) %>%
  mutate(
    Gene = factor(Gene, levels = rev(Gene))
  )

mrmr_plot_height <- max(6, nrow(mrmr_rank_df) * 0.35)

p_mrmr_rank <- ggplot(
  mrmr_rank_df,
  aes(x = Gene, y = Rank)
) +
  geom_col(
    width = 0.72,
    fill = "#1E90FF",
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  coord_flip() +
  scale_y_reverse(
    breaks = seq_len(max(mrmr_rank_df$Rank)),
    expand = expansion(mult = c(0.03, 0.05))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "mRMR selection rank",
    title = "mRMR-selected feature ranking"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_mRMR_Selection_Rank_optimized.pdf"),
  plot = p_mrmr_rank,
  width = 10,
  height = mrmr_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_mRMR_Selection_Rank_optimized.png"),
  plot = p_mrmr_rank,
  width = 10,
  height = mrmr_plot_height,
  dpi = 600,
  limitsize = FALSE
)




#sparse PLS-DA-----
library(mixOmics)
X <- as.matrix(x_df)
Y <- y_factor

keepX_value <- min(20, ncol(X))

set.seed(123)

splsda_fit <- mixOmics::splsda(
  X = X,
  Y = Y,
  ncomp = 1,
  keepX = keepX_value
)

splsda_var <- mixOmics::selectVar(
  splsda_fit,
  comp = 1
)

splsda_genes <- splsda_var$name

write.table(
  data.frame(Genes = splsda_genes),
  file = file.path(output_dir, "sPLSDA_Genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

all_genes <- add_gene_result(
  all_genes,
  "sPLSDA",
  splsda_genes,
  output_dir
)

library(ggplot2)
library(dplyr)
library(tibble)

# Prepare output directory

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Prepare sPLS-DA loading data

loading_vec <- splsda_fit$loadings$X[, 1]

splsda_loading_df <- data.frame(
  Gene = names(loading_vec),
  Loading = as.numeric(loading_vec),
  stringsAsFactors = FALSE
) %>%
  filter(Gene %in% splsda_genes) %>%
  mutate(
    AbsLoading = abs(Loading),
    Direction = ifelse(Loading >= 0, "Positive loading", "Negative loading")
  ) %>%
  arrange(desc(AbsLoading)) %>%
  mutate(
    Gene = factor(Gene, levels = rev(Gene))
  )

write.table(
  splsda_loading_df,
  file = file.path(output_dir, "sPLSDA_Selected_Gene_Loadings.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

# sPLS-DA selected gene loading plot

splsda_plot_height <- max(6, nrow(splsda_loading_df) * 0.35)

p_splsda_loading <- ggplot(
  splsda_loading_df,
  aes(x = Gene, y = AbsLoading)
) +
  geom_col(
    width = 0.72,
    fill = "#1E90FF",
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Absolute loading",
    title = "sPLS-DA selected gene loadings"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_sPLSDA_Selected_Gene_Loadings.pdf"),
  plot = p_splsda_loading,
  width = 10,
  height = splsda_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_sPLSDA_Selected_Gene_Loadings.png"),
  plot = p_splsda_loading,
  width = 10,
  height = splsda_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# sPLS-DA signed loading plot

p_splsda_signed <- ggplot(
  splsda_loading_df,
  aes(x = Gene, y = Loading)
) +
  geom_col(
    aes(fill = Direction),
    width = 0.72,
    color = "black",
    size = 0.25
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Positive loading" = "#D55E00",
      "Negative loading" = "#0072B2"
    )
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "sPLS-DA loading",
    title = "sPLS-DA signed gene loadings"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "top",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_sPLSDA_Signed_Gene_Loadings.pdf"),
  plot = p_splsda_signed,
  width = 10,
  height = splsda_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_sPLSDA_Signed_Gene_Loadings.png"),
  plot = p_splsda_signed,
  width = 10,
  height = splsda_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# Prepare sample score data

splsda_score_df <- data.frame(
  Sample = rownames(X),
  Component1 = as.numeric(splsda_fit$variates$X[, 1]),
  Group = as.factor(Y),
  stringsAsFactors = FALSE
)

write.table(
  splsda_score_df,
  file = file.path(output_dir, "sPLSDA_Sample_Component1_Scores.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

# sPLS-DA component 1 sample distribution plot

p_splsda_score <- ggplot(
  splsda_score_df,
  aes(x = Group, y = Component1, fill = Group)
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    color = "black",
    alpha = 0.75
  ) +
  geom_jitter(
    aes(color = Group),
    width = 0.12,
    size = 3,
    alpha = 0.9
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "sPLS-DA component 1 score",
    title = "sPLS-DA sample distribution"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_sPLSDA_Component1_Sample_Distribution.pdf"),
  plot = p_splsda_score,
  width = 7,
  height = 6,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_sPLSDA_Component1_Sample_Distribution.png"),
  plot = p_splsda_score,
  width = 7,
  height = 6,
  dpi = 600,
  limitsize = FALSE
)


#VSURF-----
if (!requireNamespace("VSURF", quietly = TRUE)) {
  stop("Package 'VSURF' is required.", call. = FALSE)
}
library(VSURF)

set.seed(123)
vsurf_fit <- VSURF::VSURF(
  x = x_df,
  y = y_factor,
  ntree = 1000,
  nfor.thres = 100,
  nfor.interp = 100,
  nfor.pred = 100
)

vsurf_genes_interp <- colnames(x_df)[vsurf_fit$varselect.interp]

write.table(
  data.frame(Genes = vsurf_genes_interp),
  file = file.path(output_dir, "VSURF_interpretation_Genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)


all_genes <- add_gene_result(
  all_genes,
  "VSURF_interpretation",
  vsurf_genes_interp,
  output_dir
)

# VSURF model

if (!requireNamespace("VSURF", quietly = TRUE)) {
  stop("Package 'VSURF' is required.", call. = FALSE)
}

library(VSURF)
library(ggplot2)
library(dplyr)

# Output directory

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Run VSURF

set.seed(123)

vsurf_fit <- VSURF::VSURF(
  x = x_df,
  y = y_factor,
  ntree = 1000,
  nfor.thres = 100,
  nfor.interp = 100,
  nfor.pred = 100
)

# Extract interpretation genes only

vsurf_genes_interp <- colnames(x_df)[vsurf_fit$varselect.interp]

write.table(
  data.frame(Genes = vsurf_genes_interp),
  file = file.path(output_dir, "VSURF_interpretation_Genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

all_genes <- add_gene_result(
  all_genes,
  "VSURF_interpretation",
  vsurf_genes_interp,
  output_dir
)

cat("VSURF interpretation genes:\n", paste(vsurf_genes_interp, collapse = ", "), "\n\n")

# Prepare VSURF interpretation plotting data

if (!is.null(vsurf_fit$imp.mean.dec)) {
  
  vsurf_imp <- as.numeric(vsurf_fit$imp.mean.dec)
  
  if (!is.null(names(vsurf_fit$imp.mean.dec))) {
    vsurf_imp_genes <- names(vsurf_fit$imp.mean.dec)
  } else {
    vsurf_imp_genes <- colnames(x_df)[seq_along(vsurf_imp)]
  }
  
  vsurf_interp_imp_df <- data.frame(
    Gene = vsurf_imp_genes,
    Importance = vsurf_imp,
    stringsAsFactors = FALSE
  ) %>%
    filter(Gene %in% vsurf_genes_interp) %>%
    filter(!is.na(Importance)) %>%
    arrange(desc(Importance)) %>%
    mutate(
      Gene = factor(Gene, levels = rev(Gene))
    )
  
} else {
  
  vsurf_interp_imp_df <- data.frame(
    Gene = vsurf_genes_interp,
    Importance = rev(seq_along(vsurf_genes_interp)),
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      Gene = factor(Gene, levels = rev(Gene))
    )
}

write.table(
  vsurf_interp_imp_df,
  file = file.path(output_dir, "VSURF_interpretation_Feature_Importance.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

# VSURF interpretation feature importance plot

vsurf_plot_height <- max(6, nrow(vsurf_interp_imp_df) * 0.35)

p_vsurf_interp <- ggplot(
  vsurf_interp_imp_df,
  aes(x = Gene, y = Importance)
) +
  geom_col(
    width = 0.72,
    fill = "#1E90FF",
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  coord_flip() +
  scale_y_continuous(
    expand = expansion(mult = c(0.01, 0.05))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Variable importance",
    title = "VSURF interpretation-selected features"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_VSURF_interpretation_Feature_Importance.pdf"),
  plot = p_vsurf_interp,
  width = 10,
  height = vsurf_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_VSURF_interpretation_Feature_Importance.png"),
  plot = p_vsurf_interp,
  width = 10,
  height = vsurf_plot_height,
  dpi = 600,
  limitsize = FALSE
)








#ReliefF feature selection----
library(CORElearn)
relief_data <- data.frame(x_df, group = y_factor, check.names = FALSE)
relief_score <- CORElearn::attrEval(
  group ~ .,
  data = relief_data,
  estimator = "ReliefFequalK"
)

relief_res <- data.frame(
  gene = names(relief_score),
  ReliefF_score = as.numeric(relief_score),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(ReliefF_score))

write.csv(
  relief_res,
  file = file.path(output_dir, "ReliefF_gene_importance.csv"),
  row.names = FALSE
)

relief_genes <- relief_res %>%
  filter(ReliefF_score > mean(ReliefF_score, na.rm = TRUE)) %>%
  pull(gene)

if (length(relief_genes) == 0) {
  relief_genes <- relief_res$gene[1:min(20, nrow(relief_res))]
}

all_genes <- add_gene_result(
  all_genes,
  "ReliefF",
  relief_genes,
  output_dir
)

library(ggplot2)
library(dplyr)

# Prepare ReliefF plotting data

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

relief_threshold <- mean(relief_res$ReliefF_score, na.rm = TRUE)

relief_top_n <- min(20, nrow(relief_res))

relief_plot_df <- relief_res %>%
  mutate(
    gene = as.character(gene),
    ReliefF_score = as.numeric(ReliefF_score),
    Selected = ifelse(gene %in% relief_genes, "Selected", "Not selected")
  ) %>%
  arrange(desc(ReliefF_score)) %>%
  slice_head(n = relief_top_n) %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),
    Selected = factor(Selected, levels = c("Selected", "Not selected"))
  )

relief_plot_height <- max(6, nrow(relief_plot_df) * 0.35)

# ReliefF top-ranked feature importance plot

p_relief <- ggplot(
  relief_plot_df,
  aes(x = gene, y = ReliefF_score)
) +
  geom_col(
    aes(fill = Selected),
    width = 0.72,
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  geom_hline(
    yintercept = relief_threshold,
    linetype = "dashed",
    color = "grey35",
    size = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Selected" = "#1E90FF",
      "Not selected" = "grey75"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "ReliefF score",
    title = "ReliefF feature importance"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "top",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_ReliefF_Feature_Importance.pdf"),
  plot = p_relief,
  width = 10,
  height = relief_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_ReliefF_Feature_Importance.png"),
  plot = p_relief,
  width = 10,
  height = relief_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# ReliefF selected genes only

relief_selected_plot_df <- relief_res %>%
  filter(gene %in% relief_genes) %>%
  mutate(
    gene = as.character(gene),
    ReliefF_score = as.numeric(ReliefF_score)
  ) %>%
  arrange(desc(ReliefF_score)) %>%
  mutate(
    gene = factor(gene, levels = rev(gene))
  )

if (nrow(relief_selected_plot_df) > 0) {
  
  relief_selected_height <- max(6, nrow(relief_selected_plot_df) * 0.35)
  
  p_relief_selected <- ggplot(
    relief_selected_plot_df,
    aes(x = gene, y = ReliefF_score)
  ) +
    geom_col(
      width = 0.72,
      fill = "#1E90FF",
      color = "black",
      size = 0.25
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "ReliefF score",
      title = "ReliefF-selected features"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_ReliefF_Selected_Features.pdf"),
    plot = p_relief_selected,
    width = 10,
    height = relief_selected_height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_ReliefF_Selected_Features.png"),
    plot = p_relief_selected,
    width = 10,
    height = relief_selected_height,
    dpi = 600,
    limitsize = FALSE
  )
}

#ranger permutation importance-----
library(ranger)
ranger_data <- data.frame(x_df, group = y_factor, check.names = FALSE)
set.seed(123)
ranger_fit <- ranger::ranger(
  group ~ .,
  data = ranger_data,
  num.trees = 2000,
  importance = "permutation",
  probability = TRUE,
  seed = 123
)

ranger_imp <- data.frame(
  gene = names(ranger_fit$variable.importance),
  importance = as.numeric(ranger_fit$variable.importance),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(importance))

write.csv(
  ranger_imp,
  file = file.path(output_dir, "Ranger_permutation_importance.csv"),
  row.names = FALSE
)

positive_imp <- ranger_imp$importance[ranger_imp$importance > 0]

ranger_genes <- ranger_imp %>%
  filter(importance > mean(positive_imp, na.rm = TRUE)) %>%
  pull(gene)

if (length(ranger_genes) == 0) {
  ranger_genes <- ranger_imp$gene[1:min(20, nrow(ranger_imp))]
}

all_genes <- add_gene_result(
  all_genes,
  "Ranger_RF_Permutation",
  ranger_genes,
  output_dir
)

library(ggplot2)
library(dplyr)

# Prepare Ranger plotting data

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

positive_imp <- ranger_imp$importance[ranger_imp$importance > 0]
ranger_threshold <- mean(positive_imp, na.rm = TRUE)

if (!is.finite(ranger_threshold)) {
  ranger_threshold <- mean(ranger_imp$importance, na.rm = TRUE)
}

ranger_top_n <- min(20, nrow(ranger_imp))

ranger_plot_df <- ranger_imp %>%
  mutate(
    gene = as.character(gene),
    importance = as.numeric(importance),
    Selected = ifelse(gene %in% ranger_genes, "Selected", "Not selected")
  ) %>%
  arrange(desc(importance)) %>%
  slice_head(n = ranger_top_n) %>%
  mutate(
    gene = factor(gene, levels = rev(gene)),
    Selected = factor(Selected, levels = c("Selected", "Not selected"))
  )

ranger_plot_height <- max(6, nrow(ranger_plot_df) * 0.35)

# Ranger permutation importance plot

p_ranger <- ggplot(
  ranger_plot_df,
  aes(x = gene, y = importance)
) +
  geom_col(
    aes(fill = Selected),
    width = 0.72,
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  {
    if (is.finite(ranger_threshold)) {
      geom_hline(
        yintercept = ranger_threshold,
        linetype = "dashed",
        color = "grey35",
        size = 0.7
      )
    }
  } +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Selected" = "#1E90FF",
      "Not selected" = "grey75"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.05, 0.08))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Permutation importance",
    title = "Ranger random forest feature importance"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "top",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_Ranger_RF_Permutation_Importance.pdf"),
  plot = p_ranger,
  width = 10,
  height = ranger_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_Ranger_RF_Permutation_Importance.png"),
  plot = p_ranger,
  width = 10,
  height = ranger_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# Ranger selected genes only plot

ranger_selected_plot_df <- ranger_imp %>%
  filter(gene %in% ranger_genes) %>%
  mutate(
    gene = as.character(gene),
    importance = as.numeric(importance)
  ) %>%
  arrange(desc(importance)) %>%
  mutate(
    gene = factor(gene, levels = rev(gene))
  )

if (nrow(ranger_selected_plot_df) > 0) {
  
  ranger_selected_height <- max(6, nrow(ranger_selected_plot_df) * 0.35)
  
  p_ranger_selected <- ggplot(
    ranger_selected_plot_df,
    aes(x = gene, y = importance)
  ) +
    geom_col(
      width = 0.72,
      fill = "#1E90FF",
      color = "black",
      size = 0.25
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.08))
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "Permutation importance",
      title = "Ranger-selected features"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_Ranger_RF_Selected_Features.pdf"),
    plot = p_ranger_selected,
    width = 10,
    height = ranger_selected_height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_Ranger_RF_Selected_Features.png"),
    plot = p_ranger_selected,
    width = 10,
    height = ranger_selected_height,
    dpi = 600,
    limitsize = FALSE
  )
}



#ROC-AUC filter----
library(pROC)
roc_res <- lapply(colnames(x_df), function(gene) {
  v <- x_df[[gene]]
auc_value <- tryCatch({
    roc_obj <- pROC::roc(
      response = y01,
      predictor = v,
      quiet = TRUE,
      direction = "auto"
    )
    as.numeric(pROC::auc(roc_obj))
  }, error = function(e) NA_real_)
data.frame(
    gene = gene,
    AUC = auc_value,
    AUC_direction_free = ifelse(is.na(auc_value), NA, max(auc_value, 1 - auc_value)),
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  arrange(desc(AUC_direction_free))

write.csv(
  roc_res,
  file = file.path(output_dir, "ROC_AUC_gene_statistics.csv"),
  row.names = FALSE
)

roc_genes <- roc_res %>%
  filter(AUC_direction_free >= 0.80) %>%
  pull(gene)

if (length(roc_genes) == 0) {
  roc_genes <- roc_res$gene[1:min(20, nrow(roc_res))]
}

all_genes <- add_gene_result(
  all_genes,
  "ROC_AUC",
  roc_genes,
  output_dir
)


library(ggplot2)
library(dplyr)
library(pROC)

# Prepare response variable

roc_response <- y01

if (is.factor(roc_response) || is.character(roc_response)) {
  roc_response <- ifelse(as.character(roc_response) == "CAA", 1, 0)
}

roc_response <- as.numeric(roc_response)

# Select top genes for ROC curve plotting

roc_curve_top_n <- min(6, nrow(roc_res))

roc_curve_genes <- roc_res %>%
  filter(!is.na(AUC_direction_free)) %>%
  arrange(desc(AUC_direction_free)) %>%
  slice_head(n = roc_curve_top_n) %>%
  pull(gene)

# Build corrected ROC curve data

roc_curve_list <- lapply(
  roc_curve_genes,
  function(gene) {
    
    predictor_raw <- as.numeric(x_df[[gene]])
    
    tmp_df <- data.frame(
      response = roc_response,
      predictor = predictor_raw
    ) %>%
      filter(
        !is.na(response),
        !is.na(predictor)
      )
    
    if (length(unique(tmp_df$response)) != 2) {
      return(NULL)
    }
    
    if (length(unique(tmp_df$predictor)) < 2) {
      return(NULL)
    }
    
    roc_obj <- tryCatch(
      {
        pROC::roc(
          response = tmp_df$response,
          predictor = tmp_df$predictor,
          levels = c(0, 1),
          direction = "auto",
          quiet = TRUE
        )
      },
      error = function(e) {
        NULL
      }
    )
    
    if (is.null(roc_obj)) {
      return(NULL)
    }
    
    auc_value <- as.numeric(pROC::auc(roc_obj))
    
    if (auc_value < 0.5) {
      roc_obj <- pROC::roc(
        response = tmp_df$response,
        predictor = -tmp_df$predictor,
        levels = c(0, 1),
        direction = "auto",
        quiet = TRUE
      )
      auc_value <- as.numeric(pROC::auc(roc_obj))
    }
    
    roc_coords <- as.data.frame(
      pROC::coords(
        roc_obj,
        x = "all",
        ret = c("specificity", "sensitivity"),
        transpose = FALSE
      )
    )
    
    roc_curve_df <- roc_coords %>%
      transmute(
        FPR = 1 - as.numeric(specificity),
        TPR = as.numeric(sensitivity)
      ) %>%
      filter(
        is.finite(FPR),
        is.finite(TPR)
      ) %>%
      bind_rows(
        data.frame(FPR = 0, TPR = 0),
        data.frame(FPR = 1, TPR = 1)
      ) %>%
      arrange(FPR, TPR) %>%
      mutate(
        gene = gene,
        AUC = auc_value,
        Curve_label = paste0(gene, " (AUC=", sprintf("%.3f", auc_value), ")")
      )
    
    roc_curve_df
  }
)

roc_curve_df <- bind_rows(roc_curve_list)

# Draw corrected ROC curves

if (nrow(roc_curve_df) > 0) {
  
  p_roc_curve_fixed <- ggplot(
    roc_curve_df,
    aes(
      x = FPR,
      y = TPR,
      color = Curve_label,
      group = Curve_label
    )
  ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "grey60",
      size = 0.7
    ) +
    geom_step(
      direction = "vh",
      size = 1.1
    ) +
    coord_equal() +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0.01, 0.01)
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = c(0.01, 0.01)
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = "1 - Specificity",
      y = "Sensitivity",
      title = "ROC curves of top-ranked genes"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text = element_text(
        size = 13,
        color = "black"
      ),
      axis.title = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      legend.title = element_blank(),
      legend.text = element_text(
        size = 11,
        color = "black"
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_ROC_Curves_Top_Genes_fixed.pdf"),
    plot = p_roc_curve_fixed,
    width = 8.5,
    height = 7,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_ROC_Curves_Top_Genes_fixed.png"),
    plot = p_roc_curve_fixed,
    width = 8.5,
    height = 7,
    dpi = 600,
    limitsize = FALSE
  )
}





#Wilcoxon filter----
wilcox_res <- lapply(colnames(x_df), function(gene) {
  v <- x_df[[gene]]
  
  p <- tryCatch(
    wilcox.test(v[y01 == 1], v[y01 == 0])$p.value,
    error = function(e) NA_real_
  )
  
  data.frame(
    gene = gene,
    mean_CAA = mean(v[y01 == 1], na.rm = TRUE),
    mean_HC  = mean(v[y01 == 0], na.rm = TRUE),
    delta_CAA_vs_HC = mean(v[y01 == 1], na.rm = TRUE) - mean(v[y01 == 0], na.rm = TRUE),
    p_value = p,
    stringsAsFactors = FALSE
  )
}) %>%
  bind_rows() %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    abs_delta = abs(delta_CAA_vs_HC)
  ) %>%
  arrange(p_adj, desc(abs_delta))

write.csv(
  wilcox_res,
  file = file.path(output_dir, "Wilcoxon_gene_statistics.csv"),
  row.names = FALSE
)

wilcox_genes <- wilcox_res %>%
  filter(p_value < 0.05) %>%
  pull(gene)


all_genes <- add_gene_result(
  all_genes,
  "Wilcoxon",
  wilcox_genes,
  output_dir
)

library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare output directory

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

if (!exists("all_genes") || !is.list(all_genes)) {
  all_genes <- list()
}

# Check input data

if (!exists("x_df")) {
  stop("Object 'x_df' was not found.")
}

if (!exists("y01")) {
  stop("Object 'y01' was not found. Please provide a binary vector: CAA = 1, HC = 0.")
}

x_df <- as.data.frame(x_df, check.names = FALSE)
y01 <- as.numeric(y01)

if (nrow(x_df) != length(y01)) {
  stop("The number of rows in x_df must be equal to the length of y01.")
}

if (!all(na.omit(unique(y01)) %in% c(0, 1))) {
  stop("y01 must be a binary vector with HC = 0 and CAA = 1.")
}

# Wilcoxon test for each gene

wilcox_res <- lapply(
  colnames(x_df),
  function(gene) {
    
    v <- as.numeric(x_df[[gene]])
    
    p <- tryCatch(
      {
        wilcox.test(
          v[y01 == 1],
          v[y01 == 0],
          exact = FALSE
        )$p.value
      },
      error = function(e) {
        NA_real_
      }
    )
    
    data.frame(
      gene = gene,
      mean_CAA = mean(v[y01 == 1], na.rm = TRUE),
      mean_HC = mean(v[y01 == 0], na.rm = TRUE),
      delta_CAA_vs_HC = mean(v[y01 == 1], na.rm = TRUE) -
        mean(v[y01 == 0], na.rm = TRUE),
      p_value = p,
      stringsAsFactors = FALSE
    )
  }
) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    abs_delta = abs(delta_CAA_vs_HC)
  ) %>%
  dplyr::arrange(p_value, dplyr::desc(abs_delta))

# Export Wilcoxon statistics

write.csv(
  wilcox_res,
  file = file.path(output_dir, "Wilcoxon_gene_statistics.csv"),
  row.names = FALSE
)

# Select Wilcoxon genes

wilcox_genes <- wilcox_res %>%
  dplyr::filter(p_value < 0.05) %>%
  dplyr::pull(gene)

write.table(
  data.frame(Genes = wilcox_genes),
  file = file.path(output_dir, "Wilcoxon_Genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

if (exists("add_gene_result")) {
  all_genes <- add_gene_result(
    all_genes,
    "Wilcoxon",
    wilcox_genes,
    output_dir
  )
} else {
  all_genes$Wilcoxon <- wilcox_genes
}

cat("Wilcoxon Selected Genes:\n", paste(wilcox_genes, collapse = ", "), "\n\n")

# Prepare Wilcoxon plotting data

wilcox_plot_top_n <- min(20, nrow(wilcox_res))

wilcox_plot_df <- wilcox_res %>%
  dplyr::mutate(
    gene = as.character(gene),
    p_value = as.numeric(p_value),
    p_adj = as.numeric(p_adj),
    delta_CAA_vs_HC = as.numeric(delta_CAA_vs_HC),
    abs_delta = abs(delta_CAA_vs_HC),
    neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin)),
    Direction = ifelse(delta_CAA_vs_HC >= 0, "Higher in CAA", "Lower in CAA"),
    Significant = ifelse(p_value < 0.05, "P < 0.05", "NS")
  ) %>%
  dplyr::arrange(p_value, dplyr::desc(abs_delta)) %>%
  dplyr::slice_head(n = wilcox_plot_top_n) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(gene)),
    Direction = factor(Direction, levels = c("Higher in CAA", "Lower in CAA")),
    Significant = factor(Significant, levels = c("P < 0.05", "NS"))
  )

wilcox_plot_height <- max(6, nrow(wilcox_plot_df) * 0.35)

# Wilcoxon delta ranking plot

p_wilcox_delta <- ggplot(
  wilcox_plot_df,
  aes(x = gene, y = delta_CAA_vs_HC)
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "grey45",
    size = 0.7
  ) +
  geom_col(
    aes(fill = Direction),
    width = 0.72,
    color = "black",
    size = 0.25
  ) +
  geom_point(
    aes(shape = Significant),
    color = "#B22222",
    size = 2.8
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Higher in CAA" = "#D55E00",
      "Lower in CAA" = "#0072B2"
    )
  ) +
  scale_shape_manual(
    values = c(
      "P < 0.05" = 16,
      "NS" = 1
    )
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Mean difference: CAA - HC",
    title = "Wilcoxon-ranked differential features"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "top",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_Wilcoxon_Delta_Ranking.pdf"),
  plot = p_wilcox_delta,
  width = 10,
  height = wilcox_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_Wilcoxon_Delta_Ranking.png"),
  plot = p_wilcox_delta,
  width = 10,
  height = wilcox_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# Wilcoxon -log10 P-value ranking plot

p_wilcox_p <- ggplot(
  wilcox_plot_df,
  aes(x = gene, y = neg_log10_p)
) +
  geom_col(
    aes(fill = Significant),
    width = 0.72,
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey35",
    size = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "P < 0.05" = "#1E90FF",
      "NS" = "grey75"
    )
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = expression(-log[10](P)),
    title = "Wilcoxon statistical significance ranking"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "top",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_Wilcoxon_Pvalue_Ranking.pdf"),
  plot = p_wilcox_p,
  width = 10,
  height = wilcox_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_Wilcoxon_Pvalue_Ranking.png"),
  plot = p_wilcox_p,
  width = 10,
  height = wilcox_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# Prepare boxplot genes

boxplot_gene_number <- 12

if (length(wilcox_genes) == 0) {
  
  boxplot_genes <- wilcox_res %>%
    dplyr::arrange(p_value, dplyr::desc(abs_delta)) %>%
    dplyr::slice_head(n = min(boxplot_gene_number, nrow(.))) %>%
    dplyr::pull(gene)
  
} else {
  
  boxplot_genes <- wilcox_res %>%
    dplyr::filter(gene %in% wilcox_genes) %>%
    dplyr::arrange(p_value, dplyr::desc(abs_delta)) %>%
    dplyr::slice_head(n = min(boxplot_gene_number, nrow(.))) %>%
    dplyr::pull(gene)
}

# Prepare boxplot data

wilcox_box_df <- x_df[, boxplot_genes, drop = FALSE] %>%
  as.data.frame(check.names = FALSE) %>%
  dplyr::mutate(
    Group = ifelse(y01 == 1, "CAA", "HC"),
    Sample = rownames(x_df)
  ) %>%
  tidyr::pivot_longer(
    cols = tidyselect::all_of(boxplot_genes),
    names_to = "gene",
    values_to = "Expression"
  ) %>%
  dplyr::mutate(
    Expression = as.numeric(Expression),
    Group = factor(Group, levels = c("HC", "CAA")),
    gene = factor(gene, levels = boxplot_genes)
  )

# Prepare P-value labels

wilcox_label_df <- wilcox_res %>%
  dplyr::filter(gene %in% boxplot_genes) %>%
  dplyr::select(gene, p_value) %>%
  dplyr::mutate(
    p_label = paste0("P = ", signif(p_value, 3))
  )

y_position_df <- wilcox_box_df %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    y_min = min(Expression, na.rm = TRUE),
    y_max = max(Expression, na.rm = TRUE),
    y_range = y_max - y_min,
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    y_pos = ifelse(
      y_range == 0,
      y_max + 0.1,
      y_max + 0.10 * y_range
    )
  )

wilcox_label_df <- wilcox_label_df %>%
  dplyr::left_join(y_position_df, by = "gene") %>%
  dplyr::mutate(
    gene = factor(gene, levels = boxplot_genes),
    x_pos = 1.5
  )

# Wilcoxon-selected gene boxplot

p_wilcox_box <- ggplot(
  wilcox_box_df,
  aes(x = Group, y = Expression, fill = Group)
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA,
    color = "black",
    alpha = 0.75
  ) +
  geom_jitter(
    aes(color = Group),
    width = 0.12,
    size = 2.5,
    alpha = 0.9
  ) +
  geom_text(
    data = wilcox_label_df,
    aes(x = x_pos, y = y_pos, label = p_label),
    inherit.aes = FALSE,
    size = 3.5,
    color = "black"
  ) +
  facet_wrap(
    ~gene,
    scales = "free_y",
    ncol = 4
  ) +
  scale_fill_manual(
    values = c(
      "HC" = "grey75",
      "CAA" = "#1E90FF"
    )
  ) +
  scale_color_manual(
    values = c(
      "HC" = "grey35",
      "CAA" = "#B22222"
    )
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = NULL,
    y = "Expression",
    title = "Expression of Wilcoxon-selected features"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text = element_text(
      size = 11,
      color = "black"
    ),
    axis.title.y = element_text(
      size = 15,
      face = "bold",
      color = "black"
    ),
    strip.text = element_text(
      size = 12,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.6,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.6,
      color = "black"
    ),
    legend.position = "none",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_Wilcoxon_Selected_Genes_Boxplot.pdf"),
  plot = p_wilcox_box,
  width = 12,
  height = max(6, ceiling(length(boxplot_genes) / 4) * 3.2),
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_Wilcoxon_Selected_Genes_Boxplot.png"),
  plot = p_wilcox_box,
  width = 12,
  height = max(6, ceiling(length(boxplot_genes) / 4) * 3.2),
  dpi = 600,
  limitsize = FALSE
)










#Stability Selection with glmnet----
library(stabs)
library(glmnet)
set.seed(123)
stab_fit <- stabs::stabsel(
  x = as.matrix(x_df),
  y = y01,
  fitfun = stabs::glmnet.lasso,
  args.fitfun = list(family = "binomial"),
  cutoff = 0.6,
  PFER = 1
)

stab_selected <- stab_fit$selected

if (is.numeric(stab_selected)) {
  stab_genes <- colnames(x_df)[stab_selected]
} else {
  stab_genes <- as.character(stab_selected)
}
stab_genes <- unique(na.omit(stab_genes))
write.table(
  data.frame(Genes = stab_genes),
  file = file.path(output_dir, "StabilitySelection_Genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

all_genes <- add_gene_result(
  all_genes,
  "StabilitySelection_glmnet",
  stab_genes,
  output_dir
)
cat("Stability Selection selected genes:\n")
cat(paste(stab_genes, collapse = ", "), "\n")

library(ggplot2)
library(dplyr)

# Prepare output directory

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Extract stability-selection probabilities

if (!is.null(stab_fit$max)) {
  
  stab_prob <- stab_fit$max
  
} else if (!is.null(stab_fit$phat)) {
  
  stab_prob <- stab_fit$phat
  
} else {
  
  stop("No selection probability information was found in stab_fit.")
}

stab_prob <- as.numeric(stab_prob)

if (!is.null(names(stab_fit$max)) && length(names(stab_fit$max)) == length(stab_prob)) {
  
  stab_gene_names <- names(stab_fit$max)
  
} else {
  
  stab_gene_names <- colnames(x_df)[seq_along(stab_prob)]
}

stab_cutoff <- if (!is.null(stab_fit$cutoff)) {
  stab_fit$cutoff
} else {
  0.6
}

stab_prob_df <- data.frame(
  gene = stab_gene_names,
  selection_probability = stab_prob,
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(!is.na(selection_probability)) %>%
  dplyr::mutate(
    Selected = ifelse(gene %in% stab_genes, "Selected", "Not selected")
  ) %>%
  dplyr::arrange(dplyr::desc(selection_probability))

write.csv(
  stab_prob_df,
  file = file.path(output_dir, "StabilitySelection_probability_statistics.csv"),
  row.names = FALSE
)

# Stability selection probability ranking plot

stab_top_n <- min(20, nrow(stab_prob_df))

stab_plot_df <- stab_prob_df %>%
  dplyr::slice_head(n = stab_top_n) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(gene)),
    Selected = factor(Selected, levels = c("Selected", "Not selected"))
  )

stab_plot_height <- max(6, nrow(stab_plot_df) * 0.35)

p_stab <- ggplot(
  stab_plot_df,
  aes(x = gene, y = selection_probability)
) +
  geom_col(
    aes(fill = Selected),
    width = 0.72,
    color = "black",
    size = 0.25
  ) +
  geom_point(
    color = "#B22222",
    size = 2.8
  ) +
  geom_hline(
    yintercept = stab_cutoff,
    linetype = "dashed",
    color = "grey35",
    size = 0.7
  ) +
  coord_flip() +
  scale_fill_manual(
    values = c(
      "Selected" = "#1E90FF",
      "Not selected" = "grey75"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Selection probability",
    title = "Stability selection feature ranking"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 12,
      color = "black"
    ),
    legend.position = "top",
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_StabilitySelection_Probability_Ranking.pdf"),
  plot = p_stab,
  width = 10,
  height = stab_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_StabilitySelection_Probability_Ranking.png"),
  plot = p_stab,
  width = 10,
  height = stab_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# Stability-selected genes only plot

stab_selected_plot_df <- stab_prob_df %>%
  dplyr::filter(gene %in% stab_genes) %>%
  dplyr::arrange(dplyr::desc(selection_probability)) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(gene))
  )

if (nrow(stab_selected_plot_df) > 0) {
  
  stab_selected_height <- max(6, nrow(stab_selected_plot_df) * 0.35)
  
  p_stab_selected <- ggplot(
    stab_selected_plot_df,
    aes(x = gene, y = selection_probability)
  ) +
    geom_col(
      width = 0.72,
      fill = "#1E90FF",
      color = "black",
      size = 0.25
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    geom_hline(
      yintercept = stab_cutoff,
      linetype = "dashed",
      color = "grey35",
      size = 0.7
    ) +
    coord_flip() +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2),
      expand = expansion(mult = c(0.01, 0.04))
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "Selection probability",
      title = "Stability-selected features"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_StabilitySelection_Selected_Features.pdf"),
    plot = p_stab_selected,
    width = 10,
    height = stab_selected_height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_StabilitySelection_Selected_Features.png"),
    plot = p_stab_selected,
    width = 10,
    height = stab_selected_height,
    dpi = 600,
    limitsize = FALSE
  )
}

# Optional default diagnostic plot

pdf(
  file = file.path(output_dir, "Model_StabilitySelection_Default_Diagnostic.pdf"),
  width = 10,
  height = 8
)

plot(stab_fit)

dev.off()


#PAM / nearest shrunken centroid-----

library(pamr)
library(dplyr)
library(ggplot2)

cat("Running PAM model...\n")

x_pam <- x_df[, apply(x_df, 2, sd, na.rm = TRUE) > 0, drop = FALSE]

cat("Number of genes used for PAM:", ncol(x_pam), "\n")
pam_data <- list(
  x = t(as.matrix(x_pam)),
  y = y_factor,
  geneid = colnames(x_pam),
  genenames = colnames(x_pam)
)

set.seed(123)
pam_fit <- pamr::pamr.train(pam_data)
pam_cv <- pamr::pamr.cv(
  fit = pam_fit,
  data = pam_data,
  nfold = min(10, min(table(y_factor)))
)

#  CV 
pam_cv_df <- data.frame(
  threshold = pam_cv$threshold,
  error = pam_cv$error,
  stringsAsFactors = FALSE
)

write.csv(
  pam_cv_df,
  file = file.path(output_dir, "PAM_CV_threshold_error.csv"),
  row.names = FALSE
)


pam_cv_df <- data.frame(
  threshold = as.numeric(unlist(pam_cv$threshold)),
  error = as.numeric(unlist(pam_cv$error)),
  stringsAsFactors = FALSE
)

# Remove NA / Inf
pam_cv_df <- pam_cv_df[
  is.finite(pam_cv_df$threshold) &
    is.finite(pam_cv_df$error),
  ,
  drop = FALSE
]

# Minimum CV error
min_error <- min(pam_cv_df$error, na.rm = TRUE)
candidate_idx <- which(pam_cv_df$error == min_error)

best_threshold <- max(
  pam_cv_df$threshold[candidate_idx],
  na.rm = TRUE
)

cat("Best PAM threshold:", best_threshold, "\n")
cat("Minimum CV error:", min_error, "\n")

pam_gene_table <- pamr::pamr.listgenes(
  fit = pam_fit,
  data = pam_data,
  threshold = best_threshold,
  genenames = TRUE
)

pam_gene_table <- as.data.frame(pam_gene_table, check.names = FALSE)

write.csv(
  pam_gene_table,
  file = file.path(output_dir, "PAM_selected_gene_table.csv"),
  row.names = FALSE
)

gene_col <- intersect(
  c("id", "geneid", "genename", "Gene", "gene", "Name"),
  colnames(pam_gene_table)
)[1]

if (!is.na(gene_col)) {
  pam_genes <- unique(as.character(pam_gene_table[[gene_col]]))
} else {
  pam_genes <- unique(rownames(pam_gene_table))
}

pam_genes <- pam_genes[!is.na(pam_genes)]
pam_genes <- pam_genes[pam_genes != ""]
pam_genes <- pam_genes[pam_genes %in% colnames(x_df)]

cat("Number of PAM selected genes:", length(pam_genes), "\n")
cat("PAM selected genes:\n")
cat(paste(pam_genes, collapse = ", "), "\n")

all_genes <- add_gene_result(
  all_genes,
  "PAM",
  pam_genes,
  output_dir
)



library(ggplot2)
library(dplyr)
library(tidyr)

# Prepare output directory

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Prepare PAM cross-validation data

pam_threshold_vec <- as.numeric(unlist(pam_cv$threshold))
pam_error_vec <- as.numeric(unlist(pam_cv$error))

pam_len <- min(length(pam_threshold_vec), length(pam_error_vec))

pam_cv_plot_df <- data.frame(
  threshold = pam_threshold_vec[seq_len(pam_len)],
  error = pam_error_vec[seq_len(pam_len)],
  stringsAsFactors = FALSE
) %>%
  dplyr::filter(
    is.finite(threshold),
    is.finite(error)
  ) %>%
  dplyr::arrange(threshold)

if (!exists("best_threshold") || !is.finite(best_threshold)) {
  
  min_error <- min(pam_cv_plot_df$error, na.rm = TRUE)
  
  best_threshold <- pam_cv_plot_df %>%
    dplyr::filter(error == min_error) %>%
    dplyr::summarise(best_threshold = max(threshold, na.rm = TRUE)) %>%
    dplyr::pull(best_threshold)
}

if (!exists("min_error") || !is.finite(min_error)) {
  min_error <- min(pam_cv_plot_df$error, na.rm = TRUE)
}

best_cv_df <- pam_cv_plot_df %>%
  dplyr::filter(threshold == best_threshold) %>%
  dplyr::slice_head(n = 1)

# PAM CV threshold-error plot

p_pam_cv <- ggplot(
  pam_cv_plot_df,
  aes(x = threshold, y = error)
) +
  geom_line(
    color = "grey35",
    size = 0.8
  ) +
  geom_point(
    color = "#B22222",
    size = 2.6
  ) +
  geom_vline(
    xintercept = best_threshold,
    linetype = "dashed",
    color = "#E41A1C",
    size = 0.8
  ) +
  geom_point(
    data = best_cv_df,
    aes(x = threshold, y = error),
    color = "#E41A1C",
    size = 3.8
  ) +
  annotate(
    geom = "text",
    x = best_threshold,
    y = min_error,
    label = paste0("Best threshold = ", signif(best_threshold, 4)),
    color = "#E41A1C",
    size = 4.2,
    hjust = -0.05,
    vjust = -0.8
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Threshold",
    y = "Cross-validation error",
    title = "PAM cross-validation threshold selection"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text = element_text(
      size = 13,
      color = "black"
    ),
    axis.title = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      size = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      size = 0.7,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      10, 25, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_PAM_CV_Threshold_Error.pdf"),
  plot = p_pam_cv,
  width = 9,
  height = 7,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_PAM_CV_Threshold_Error.png"),
  plot = p_pam_cv,
  width = 9,
  height = 7,
  dpi = 600,
  limitsize = FALSE
)

# Prepare PAM selected-gene score data

pam_gene_table_plot <- pam_gene_table %>%
  as.data.frame(check.names = FALSE)

gene_col <- intersect(
  c("id", "geneid", "genename", "Gene", "gene", "Name", "name"),
  colnames(pam_gene_table_plot)
)[1]

if (!is.na(gene_col)) {
  pam_gene_table_plot$Gene <- as.character(pam_gene_table_plot[[gene_col]])
} else {
  pam_gene_table_plot$Gene <- rownames(pam_gene_table_plot)
}

pam_gene_table_plot <- pam_gene_table_plot %>%
  dplyr::filter(
    !is.na(Gene),
    Gene != "",
    Gene %in% pam_genes
  )

numeric_cols <- colnames(pam_gene_table_plot)[
  sapply(
    pam_gene_table_plot,
    function(z) {
      suppressWarnings(
        all(is.na(z) | is.finite(as.numeric(z)))
      )
    }
  )
]

numeric_cols <- setdiff(numeric_cols, c("Gene", gene_col))

score_col_candidates <- grep(
  "score|d|centroid|overall",
  numeric_cols,
  ignore.case = TRUE,
  value = TRUE
)

if (length(score_col_candidates) > 0) {
  score_col <- score_col_candidates[1]
} else if (length(numeric_cols) > 0) {
  score_col <- numeric_cols[1]
} else {
  score_col <- NA
}

if (!is.na(score_col)) {
  
  pam_score_df <- pam_gene_table_plot %>%
    dplyr::mutate(
      PAM_score = abs(as.numeric(.data[[score_col]]))
    ) %>%
    dplyr::select(Gene, PAM_score) %>%
    dplyr::filter(
      !is.na(PAM_score),
      is.finite(PAM_score)
    ) %>%
    dplyr::arrange(dplyr::desc(PAM_score))
  
} else {
  
  pam_score_df <- data.frame(
    Gene = pam_genes,
    PAM_score = rev(seq_along(pam_genes)),
    stringsAsFactors = FALSE
  )
}

pam_score_df <- pam_score_df %>%
  dplyr::distinct(Gene, .keep_all = TRUE) %>%
  dplyr::arrange(dplyr::desc(PAM_score)) %>%
  dplyr::mutate(
    Gene = factor(Gene, levels = rev(Gene))
  )

write.csv(
  pam_score_df,
  file = file.path(output_dir, "PAM_selected_gene_score_for_plot.csv"),
  row.names = FALSE
)

# PAM selected-gene score plot

if (nrow(pam_score_df) > 0) {
  
  pam_score_height <- max(6, nrow(pam_score_df) * 0.35)
  
  p_pam_score <- ggplot(
    pam_score_df,
    aes(x = Gene, y = PAM_score)
  ) +
    geom_col(
      width = 0.72,
      fill = "#1E90FF",
      color = "black",
      size = 0.25
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = "PAM score",
      title = "PAM-selected feature ranking"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_PAM_Selected_Genes_Score.pdf"),
    plot = p_pam_score,
    width = 10,
    height = pam_score_height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_PAM_Selected_Genes_Score.png"),
    plot = p_pam_score,
    width = 10,
    height = pam_score_height,
    dpi = 600,
    limitsize = FALSE
  )
}

# Prepare PAM selected-gene boxplot data

pam_box_gene_number <- min(12, length(pam_genes))

if (pam_box_gene_number > 0) {
  
  pam_box_genes <- pam_score_df %>%
    dplyr::arrange(dplyr::desc(PAM_score)) %>%
    dplyr::slice_head(n = pam_box_gene_number) %>%
    dplyr::pull(Gene) %>%
    as.character()
  
  pam_group <- as.factor(y_factor)
  pam_group_levels <- levels(pam_group)
  
  pam_fill_values <- setNames(
    c("grey75", "#1E90FF")[seq_along(pam_group_levels)],
    pam_group_levels
  )
  
  pam_color_values <- setNames(
    c("grey35", "#B22222")[seq_along(pam_group_levels)],
    pam_group_levels
  )
  
  pam_box_df <- x_df[, pam_box_genes, drop = FALSE] %>%
    as.data.frame(check.names = FALSE) %>%
    dplyr::mutate(
      Group = pam_group,
      Sample = rownames(x_df)
    ) %>%
    tidyr::pivot_longer(
      cols = tidyselect::all_of(pam_box_genes),
      names_to = "Gene",
      values_to = "Expression"
    ) %>%
    dplyr::mutate(
      Expression = as.numeric(Expression),
      Gene = factor(Gene, levels = pam_box_genes),
      Group = factor(Group, levels = pam_group_levels)
    )
  
  p_pam_box <- ggplot(
    pam_box_df,
    aes(x = Group, y = Expression, fill = Group)
  ) +
    geom_boxplot(
      width = 0.5,
      outlier.shape = NA,
      color = "black",
      alpha = 0.75
    ) +
    geom_jitter(
      aes(color = Group),
      width = 0.12,
      size = 2.5,
      alpha = 0.9
    ) +
    facet_wrap(
      ~Gene,
      scales = "free_y",
      ncol = 4
    ) +
    scale_fill_manual(
      values = pam_fill_values
    ) +
    scale_color_manual(
      values = pam_color_values
    ) +
    theme_classic(base_size = 15) +
    labs(
      x = NULL,
      y = "Expression",
      title = "Expression of PAM-selected features"
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text = element_text(
        size = 11,
        color = "black"
      ),
      axis.title.y = element_text(
        size = 15,
        face = "bold",
        color = "black"
      ),
      strip.text = element_text(
        size = 12,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.6,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.6,
        color = "black"
      ),
      legend.position = "none",
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_PAM_Selected_Genes_Boxplot.pdf"),
    plot = p_pam_box,
    width = 12,
    height = max(6, ceiling(length(pam_box_genes) / 4) * 3.2),
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_PAM_Selected_Genes_Boxplot.png"),
    plot = p_pam_box,
    width = 12,
    height = max(6, ceiling(length(pam_box_genes) / 4) * 3.2),
    dpi = 600,
    limitsize = FALSE
  )
}

# Optional PAM default diagnostic plots

pdf(
  file = file.path(output_dir, "Model_PAM_Default_CV_Plot.pdf"),
  width = 10,
  height = 8
)

pamr::pamr.plotcv(pam_cv)

dev.off()




#Information-theoretic feature selection----
if (!requireNamespace("FSelectorRcpp", quietly = TRUE)) {
  stop("Package 'FSelectorRcpp' is required.", call. = FALSE)
}

library(FSelectorRcpp)
library(dplyr)
library(ggplot2)

cat("Running information-theoretic feature selection...\n")

# 1. Prepare data

info_data <- data.frame(
  x_df,
  group = y_factor,
  check.names = FALSE
)

info_data$group <- as.factor(info_data$group)

# 2. Information Gain

ig_res <- FSelectorRcpp::information_gain(
  group ~ .,
  data = info_data,
  type = "infogain"
) %>%
  dplyr::arrange(dplyr::desc(importance))

write.csv(
  ig_res,
  file = file.path(output_dir, "InformationGain_gene_importance.csv"),
  row.names = FALSE
)

ig_genes <- ig_res %>%
  dplyr::filter(importance > mean(importance, na.rm = TRUE)) %>%
  dplyr::pull(attributes)

if (length(ig_genes) == 0) {
  ig_genes <- ig_res$attributes[1:min(20, nrow(ig_res))]
}

all_genes <- add_gene_result(
  all_genes,
  "InformationGain",
  ig_genes,
  output_dir
)

cat("Information Gain selected genes:\n")
cat(paste(ig_genes, collapse = ", "), "\n\n")


# 3. Gain Ratio

gr_res <- FSelectorRcpp::information_gain(
  group ~ .,
  data = info_data,
  type = "gainratio"
) %>%
  dplyr::arrange(dplyr::desc(importance))

write.csv(
  gr_res,
  file = file.path(output_dir, "GainRatio_gene_importance.csv"),
  row.names = FALSE
)

gr_genes <- gr_res %>%
  dplyr::filter(importance > mean(importance, na.rm = TRUE)) %>%
  dplyr::pull(attributes)

if (length(gr_genes) == 0) {
  gr_genes <- gr_res$attributes[1:min(20, nrow(gr_res))]
}

all_genes <- add_gene_result(
  all_genes,
  "GainRatio",
  gr_genes,
  output_dir
)

cat("Gain Ratio selected genes:\n")
cat(paste(gr_genes, collapse = ", "), "\n\n")

# Symmetrical Uncertainty

su_res <- FSelectorRcpp::information_gain(
  group ~ .,
  data = info_data,
  type = "symuncert"
) %>%
  dplyr::arrange(dplyr::desc(importance))

write.csv(
  su_res,
  file = file.path(output_dir, "SymmetricalUncertainty_gene_importance.csv"),
  row.names = FALSE
)

su_genes <- su_res %>%
  dplyr::filter(importance > mean(importance, na.rm = TRUE)) %>%
  dplyr::pull(attributes)

if (length(su_genes) == 0) {
  su_genes <- su_res$attributes[1:min(20, nrow(su_res))]
}

all_genes <- add_gene_result(
  all_genes,
  "SymmetricalUncertainty",
  su_genes,
  output_dir
)

cat("Symmetrical Uncertainty selected genes:\n")
cat(paste(su_genes, collapse = ", "), "\n\n")

##
plot_info_filter <- function(
    res,
    method_name,
    output_dir,
    top_n = 20
) {
  
  library(dplyr)
  library(ggplot2)
  
 
  res <- as.data.frame(res, check.names = FALSE)
  
  if (!all(c("attributes", "importance") %in% colnames(res))) {
    stop("Input result must contain columns: attributes and importance.")
  }
  
  plot_df <- res %>%
    dplyr::select(
      gene = attributes,
      importance = importance
    ) %>%
    dplyr::filter(
      !is.na(gene),
      gene != "",
      !is.na(importance),
      is.finite(importance)
    ) %>%
    dplyr::arrange(dplyr::desc(importance))
  
  if (nrow(plot_df) == 0) {
    stop("No valid features available for plotting.")
  }
  

  top_k <- min(top_n, nrow(plot_df))
  cutoff_value <- mean(plot_df$importance, na.rm = TRUE)
  
  plot_df <- plot_df[seq_len(top_k), , drop = FALSE]
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      selected = ifelse(
        importance > cutoff_value,
        "Selected",
        "Other"
      )
    ) %>%
    dplyr::arrange(importance) %>%
    dplyr::mutate(
      gene = factor(gene, levels = gene)
    )
  

  p <- ggplot(
    plot_df,
    aes(x = gene, y = importance, fill = selected)
  ) +
    geom_col(
      width = 0.65,
      color = "black",
      linewidth = 0.35
    ) +
    geom_point(
      shape = 21,
      fill = "#F39C12",
      color = "black",
      size = 3,
      stroke = 0.4
    ) +
    coord_flip() +
    scale_fill_manual(
      values = c(
        "Selected" = "#D55E00",
        "Other" = "grey75"
      )
    ) +
    labs(
      x = "Gene",
      y = "Importance score",
      fill = "",
      title = paste0(method_name, " feature ranking")
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(
        color = "black",
        size = 16,
        face = "bold",
        hjust = 0.5
      ),
      axis.text = element_text(
        color = "black",
        size = 11
      ),
      axis.title = element_text(
        color = "black",
        size = 14,
        face = "bold"
      ),
      legend.position = "right",
      legend.text = element_text(
        color = "black",
        size = 11
      ),
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.7
      ),
      axis.line = element_blank(),
      panel.grid.major.x = element_line(
        color = "grey90",
        linewidth = 0.3
      ),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
  

  ggsave(
    filename = file.path(output_dir, paste0(method_name, "_ranking_plot.pdf")),
    plot = p,
    width = 7,
    height = 5.8
  )
  
  ggsave(
    filename = file.path(output_dir, paste0(method_name, "_ranking_plot.png")),
    plot = p,
    width = 7,
    height = 5.8,
    dpi = 600
  )
  
  return(p)
}


p_ig <- plot_info_filter(
  res = ig_res,
  method_name = "InformationGain",
  output_dir = output_dir
)

p_gr <- plot_info_filter(
  res = gr_res,
  method_name = "GainRatio",
  output_dir = output_dir
)

p_su <- plot_info_filter(
  res = su_res,
  method_name = "SymmetricalUncertainty",
  output_dir = output_dir
)




#Genetic Algorithm Feature Selection-------
library(caret)
set.seed(123)

ga_ctrl <- gafsControl(
  functions = caretGA,
  method = "cv",
  number = min(5, min(table(y_factor)))
)

ga_fit <- gafs(
  x = x_df,
  y = y_factor,
  iters = 10,
  popSize = 10,
  gafsControl = ga_ctrl,
  method = "rf"
)

ga_genes <- ga_fit$optVariables
ga_genes <- unique(na.omit(ga_genes))

write.table(
  data.frame(Genes = ga_genes),
  file = file.path(output_dir, "GA_FeatureSelection_Genes.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

all_genes <- add_gene_result(
  all_genes,
  "GeneticAlgorithm",
  ga_genes,
  output_dir
)

cat("Genetic Algorithm selected genes:\n")
cat(paste(ga_genes, collapse = ", "), "\n")


library(ggplot2)
library(dplyr)
library(caret)

# Prepare output directory

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Prepare GA-selected feature importance

ga_imp_df <- tryCatch(
  {
    vi <- caret::varImp(ga_fit$fit, scale = TRUE)
    vi_df <- as.data.frame(vi$importance)
    
    numeric_cols <- colnames(vi_df)[
      vapply(vi_df, is.numeric, logical(1))
    ]
    
    if ("Overall" %in% colnames(vi_df)) {
      score_vec <- vi_df$Overall
    } else if (length(numeric_cols) > 0) {
      score_vec <- rowMeans(vi_df[, numeric_cols, drop = FALSE], na.rm = TRUE)
    } else {
      stop("No numeric importance column was found.")
    }
    
    data.frame(
      gene = rownames(vi_df),
      importance = as.numeric(score_vec),
      stringsAsFactors = FALSE
    )
  },
  error = function(e) {
    NULL
  }
)

if (!is.null(ga_imp_df)) {
  
  ga_plot_df <- ga_imp_df %>%
    dplyr::filter(gene %in% ga_genes) %>%
    dplyr::filter(!is.na(importance)) %>%
    dplyr::arrange(dplyr::desc(importance))
  
  plot_y_label <- "RF importance in GA-selected model"
  plot_title <- "Genetic algorithm-selected feature importance"
  
} else {
  
  ga_plot_df <- data.frame(
    gene = ga_genes,
    importance = rev(seq_along(ga_genes)),
    stringsAsFactors = FALSE
  )
  
  plot_y_label <- "Selection order score"
  plot_title <- "Genetic algorithm-selected features"
}

ga_plot_df <- ga_plot_df %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(gene))
  )

write.csv(
  ga_plot_df,
  file = file.path(output_dir, "GA_Selected_Feature_Importance_for_plot.csv"),
  row.names = FALSE
)

# GA-selected feature plot

if (nrow(ga_plot_df) > 0) {
  
  ga_plot_height <- max(6, nrow(ga_plot_df) * 0.35)
  
  p_ga_features <- ggplot(
    ga_plot_df,
    aes(x = gene, y = importance)
  ) +
    geom_col(
      width = 0.72,
      fill = "#1E90FF",
      color = "black",
      size = 0.25
    ) +
    geom_point(
      color = "#B22222",
      size = 2.8
    ) +
    coord_flip() +
    scale_y_continuous(
      expand = expansion(mult = c(0.01, 0.05))
    ) +
    theme_classic(base_size = 16) +
    labs(
      x = NULL,
      y = plot_y_label,
      title = plot_title
    ) +
    theme(
      plot.title = element_text(
        size = 18,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      axis.text.y = element_text(
        size = 13,
        color = "black"
      ),
      axis.text.x = element_text(
        size = 13,
        color = "black"
      ),
      axis.title.x = element_text(
        size = 16,
        face = "bold",
        color = "black"
      ),
      axis.line = element_line(
        size = 0.7,
        color = "black"
      ),
      axis.ticks = element_line(
        size = 0.7,
        color = "black"
      ),
      plot.margin = ggplot2::margin(
        10, 20, 10, 10,
        unit = "pt"
      )
    )
  
  ggsave(
    filename = file.path(output_dir, "Model_GA_Selected_Features.pdf"),
    plot = p_ga_features,
    width = 10,
    height = ga_plot_height,
    limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(output_dir, "Model_GA_Selected_Features.png"),
    plot = p_ga_features,
    width = 10,
    height = ga_plot_height,
    dpi = 600,
    limitsize = FALSE
  )
}

# Prepare GA performance data if available

ga_perf_raw <- NULL

if (!is.null(ga_fit$external) && is.data.frame(ga_fit$external)) {
  ga_perf_raw <- ga_fit$external
} else if (!is.null(ga_fit$ga) && is.data.frame(ga_fit$ga)) {
  ga_perf_raw <- ga_fit$ga
} else if (!is.null(ga_fit$results) && is.data.frame(ga_fit$results)) {
  ga_perf_raw <- ga_fit$results
}

if (!is.null(ga_perf_raw)) {
  
  ga_perf_raw <- as.data.frame(ga_perf_raw)
  
  possible_iter_cols <- intersect(
    c("Iter", "Iteration", "iter", "generation", "Generation"),
    colnames(ga_perf_raw)
  )
  
  if (length(possible_iter_cols) > 0) {
    iter_col <- possible_iter_cols[1]
    ga_perf_raw$GA_iteration <- suppressWarnings(as.numeric(ga_perf_raw[[iter_col]]))
  } else {
    ga_perf_raw$GA_iteration <- seq_len(nrow(ga_perf_raw))
  }
  
  numeric_cols <- colnames(ga_perf_raw)[
    vapply(ga_perf_raw, is.numeric, logical(1))
  ]
  
  metric_candidates <- intersect(
    c("Accuracy", "ROC", "Kappa", "Fitness", "fitness", "Sens", "Spec"),
    numeric_cols
  )
  
  metric_candidates <- setdiff(metric_candidates, "GA_iteration")
  
  if (length(metric_candidates) > 0) {
    metric_col <- metric_candidates[1]
  } else {
    metric_col <- setdiff(numeric_cols, "GA_iteration")[1]
  }
  
  if (!is.na(metric_col)) {
    
    ga_perf_df <- ga_perf_raw %>%
      dplyr::mutate(
        GA_iteration = as.numeric(GA_iteration),
        Performance = as.numeric(.data[[metric_col]])
      ) %>%
      dplyr::filter(
        is.finite(GA_iteration),
        is.finite(Performance)
      ) %>%
      dplyr::group_by(GA_iteration) %>%
      dplyr::summarise(
        mean_performance = mean(Performance, na.rm = TRUE),
        sd_performance = sd(Performance, na.rm = TRUE),
        .groups = "drop"
      )
    
    write.csv(
      ga_perf_df,
      file = file.path(output_dir, "GA_Performance_Curve_for_plot.csv"),
      row.names = FALSE
    )
    
    if (nrow(ga_perf_df) > 1) {
      
      p_ga_perf <- ggplot(
        ga_perf_df,
        aes(x = GA_iteration, y = mean_performance)
      ) +
        geom_line(
          color = "grey35",
          size = 0.8
        ) +
        geom_point(
          color = "#B22222",
          size = 2.6
        ) +
        theme_classic(base_size = 16) +
        labs(
          x = "GA iteration",
          y = metric_col,
          title = "Genetic algorithm performance"
        ) +
        theme(
          plot.title = element_text(
            size = 18,
            face = "bold",
            hjust = 0.5,
            color = "black"
          ),
          axis.text = element_text(
            size = 13,
            color = "black"
          ),
          axis.title = element_text(
            size = 16,
            face = "bold",
            color = "black"
          ),
          axis.line = element_line(
            size = 0.7,
            color = "black"
          ),
          axis.ticks = element_line(
            size = 0.7,
            color = "black"
          ),
          plot.margin = ggplot2::margin(
            10, 20, 10, 10,
            unit = "pt"
          )
        )
      
      ggsave(
        filename = file.path(output_dir, "Model_GA_Performance_Curve.pdf"),
        plot = p_ga_perf,
        width = 9,
        height = 7,
        limitsize = FALSE
      )
      
      ggsave(
        filename = file.path(output_dir, "Model_GA_Performance_Curve.png"),
        plot = p_ga_perf,
        width = 9,
        height = 7,
        dpi = 600,
        limitsize = FALSE
      )
    }
  }
}

# Optional default GA diagnostic plot

pdf(
  file = file.path(output_dir, "Model_GA_Default_Diagnostic.pdf"),
  width = 10,
  height = 8
)

try(
  plot(ga_fit),
  silent = TRUE
)

dev.off()




#Simulated Annealing Feature Selection###----
#5020%
library(caret)
library(dplyr)
library(tidyr)
library(ggplot2)

n_repeats <- 50     # repeated SA runs
sa_iters  <- 5      # iterations inside each SA run
sa_result_list <- list()
cat("Running SA repeat", i, "...\n")
set.seed(1000 + i)
  sa_fit_tmp <- safs(
    x = x_df,
    y = y_factor,
    iters = sa_iters,
    safsControl = sa_ctrl,
    method = "rf"
  )
  
  sa_genes_tmp <- unique(na.omit(as.character(sa_fit_tmp$optVariables)))
  

  
  cat("Selected genes:", paste(sa_genes_tmp, collapse = ", "), "\n\n")


all_genes$SimulatedAnnealing_single <- sa_genes_tmp 

  library(dplyr)
library(ggplot2)

# Prepare output directory

if (!exists("output_dir")) {
  output_dir <- "ML_screening_results"
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Check existing SA results

if (!exists("sa_result_list")) {
  stop("Object 'sa_result_list' was not found. Please run SA feature selection first.")
}

sa_result_list <- sa_result_list[
  vapply(sa_result_list, function(z) length(z) > 0, logical(1))
]

if (length(sa_result_list) == 0) {
  stop("sa_result_list is empty. No SA-selected genes were found.")
}

n_trials <- length(sa_result_list)

cat("Number of available SA runs:", n_trials, "\n")

# Convert existing SA results into long format

sa_long <- lapply(seq_along(sa_result_list), function(i) {
  
  genes_i <- unique(na.omit(as.character(sa_result_list[[i]])))
  
  if (length(genes_i) == 0) {
    return(NULL)
  }
  
  data.frame(
    trial = i,
    gene = genes_i,
    selected = 1,
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

if (nrow(sa_long) == 0) {
  stop("No selected genes were found in sa_result_list.")
}

# Summarize selection frequency

sa_freq <- sa_long %>%
  dplyr::count(gene, name = "n_selected") %>%
  dplyr::mutate(
    selection_frequency = n_selected / n_trials
  ) %>%
  dplyr::arrange(
    dplyr::desc(selection_frequency),
    gene
  )

write.csv(
  sa_freq,
  file = file.path(output_dir, "SA_gene_selection_summary_existing_results.csv"),
  row.names = FALSE
)

sa_genes <- sa_freq$gene

write.table(
  data.frame(Genes = sa_genes),
  file = file.path(output_dir, "SA_Selected_Genes_existing_results.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

if (exists("add_gene_result")) {
  all_genes <- add_gene_result(
    all_genes,
    "SimulatedAnnealing",
    sa_genes,
    output_dir
  )
} else {
  if (!exists("all_genes") || !is.list(all_genes)) {
    all_genes <- list()
  }
  all_genes$SimulatedAnnealing <- sa_genes
}

cat("SA selected genes from existing results:\n")
cat(paste(sa_genes, collapse = ", "), "\n\n")

# Plot data

top_plot_n <- 30

sa_plot_df <- sa_freq %>%
  dplyr::slice_head(n = min(top_plot_n, nrow(sa_freq))) %>%
  dplyr::arrange(
    dplyr::desc(selection_frequency),
    gene
  ) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(gene))
  )

sa_plot_height <- max(5, nrow(sa_plot_df) * 0.35)

# Clean lollipop plot without special DUSP2 highlighting

p_sa_lollipop <- ggplot(
  sa_plot_df,
  aes(x = selection_frequency, y = gene)
) +
  geom_segment(
    aes(
      x = 0,
      xend = selection_frequency,
      y = gene,
      yend = gene
    ),
    color = "grey65",
    linewidth = 0.8
  ) +
  geom_point(
    shape = 21,
    fill = "#1E90FF",
    color = "black",
    size = 4,
    stroke = 0.35
  ) +
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "Selection frequency",
    y = NULL,
    title = "Simulated annealing-selected features"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.7,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_SA_Selected_Features_Lollipop_Normal.pdf"),
  plot = p_sa_lollipop,
  width = 8,
  height = sa_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_SA_Selected_Features_Lollipop_Normal.png"),
  plot = p_sa_lollipop,
  width = 8,
  height = sa_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# Clean barplot version without special DUSP2 highlighting

p_sa_bar <- ggplot(
  sa_plot_df,
  aes(x = gene, y = selection_frequency)
) +
  geom_col(
    width = 0.72,
    fill = "#1E90FF",
    color = "black",
    linewidth = 0.25
  ) +
  coord_flip() +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.2),
    expand = expansion(mult = c(0.01, 0.04))
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = NULL,
    y = "Selection frequency",
    title = "Simulated annealing-selected features"
  ) +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.text.y = element_text(
      size = 13,
      color = "black"
    ),
    axis.text.x = element_text(
      size = 13,
      color = "black"
    ),
    axis.title.x = element_text(
      size = 16,
      face = "bold",
      color = "black"
    ),
    axis.line = element_line(
      linewidth = 0.7,
      color = "black"
    ),
    axis.ticks = element_line(
      linewidth = 0.7,
      color = "black"
    ),
    plot.margin = ggplot2::margin(
      10, 20, 10, 10,
      unit = "pt"
    )
  )

ggsave(
  filename = file.path(output_dir, "Model_SA_Selected_Features_Barplot_Normal.pdf"),
  plot = p_sa_bar,
  width = 8,
  height = sa_plot_height,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Model_SA_Selected_Features_Barplot_Normal.png"),
  plot = p_sa_bar,
  width = 8,
  height = sa_plot_height,
  dpi = 600,
  limitsize = FALSE
)

# SimulatedAnnealing_single is already stored above from sa_genes_tmp.
##-----
# -----
####

# Final one-click workflow:
# Multi-algorithm consensus feature-selection visualization
# No UpSet; no target-gene highlighting
# Class blocks ordered by number of algorithms
# Large-font Nature-style output

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "stringr",
  "scales",
  "patchwork"
)

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.", call. = FALSE)
  }
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
  library(scales)
  library(patchwork)
  library(grid)
})



# 1. Basic settings


if (!exists("output_dir")) {
  output_dir <- "ML_screening_results"
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

top_n_vote <- 30
heatmap_top_n <- 40


# Prepare all_genes

if (!exists("all_genes")) {
  all_genes <- list()
}

if (is.data.frame(all_genes)) {
  all_genes <- as.list(all_genes)
}

if (!is.list(all_genes)) {
  stop("all_genes must be a list or data.frame.")
}


#3. Automatically add existing gene-selection objects


candidate_sources <- c(
  "Wilcoxon" = "wilcox_genes",
  "ROC_AUC" = "roc_genes",
  "RandomForest" = "rf_genes",
  "Ranger_RF_Permutation" = "ranger_genes",
  "Boruta" = "boruta_genes",
  "VSURF_interpretation" = "vsurf_genes_interp",
  "VSURF_prediction" = "vsurf_genes_pred",
  "XGBoost" = "xgb_genes",
  "LightGBM" = "lgb_genes",
  "LightGBM_SHAP" = "lgb_shap_genes_topN",
  "CatBoost" = "catboost_genes",
  "LASSO" = "lasso_genes",
  "Ridge" = "ridge_genes",
  "ElasticNet" = "enet_genes",
  "StabilitySelection_glmnet" = "stab_genes",
  "PAM" = "pam_genes",
  "SVM" = "svm_genes",
  "SVM_RFE" = "svm_rfe_genes",
  "mRMR" = "mrmr_genes",
  "InformationGain" = "ig_genes",
  "GainRatio" = "gr_genes",
  "SymmetricalUncertainty" = "su_genes",
  "ReliefF" = "relief_genes",
  "sPLSDA" = "splsda_genes",
  "GeneticAlgorithm" = "ga_genes",
  "SimulatedAnnealing" = "sa_genes"
)

for (method_name in names(candidate_sources)) {
  
  obj_name <- candidate_sources[[method_name]]
  
  if (exists(obj_name, inherits = TRUE)) {
    
    tmp_genes <- get(obj_name, inherits = TRUE)
    tmp_genes <- unique(na.omit(as.character(unlist(tmp_genes))))
    tmp_genes <- tmp_genes[tmp_genes != ""]
    tmp_genes <- tmp_genes[tmp_genes != "NA"]
    tmp_genes <- tmp_genes[tmp_genes != "NULL"]
    
    if (length(tmp_genes) > 0) {
      if (is.null(names(all_genes)) || !method_name %in% names(all_genes)) {
        all_genes[[method_name]] <- tmp_genes
      }
    }
  }
}


#  4. Robust gene extraction function

extract_genes_from_object <- function(x) {
  
  if (is.null(x)) {
    return(character(0))
  }
  
  if (is.data.frame(x)) {
    
    possible_cols <- c(
      "Genes", "genes", "Gene", "gene",
      "Feature", "feature", "features",
      "selected.fea", "selected_features",
      "ID", "id",
      "all", "All"
    )
    
    gene_col <- intersect(possible_cols, colnames(x))[1]
    
    if (!is.na(gene_col)) {
      genes <- x[[gene_col]]
    } else {
      genes <- x[[1]]
    }
    
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
  
  return(genes)
}


# 5. Convert all_genes to long-format voting table

method_names <- names(all_genes)

if (is.null(method_names)) {
  method_names <- paste0("Method_", seq_along(all_genes))
}

method_names <- make.unique(method_names)

vote_long <- lapply(seq_along(all_genes), function(i) {
  
  genes <- extract_genes_from_object(all_genes[[i]])
  
  if (length(genes) == 0) {
    return(NULL)
  }
  
  data.frame(
    method = method_names[i],
    gene = genes,
    stringsAsFactors = FALSE
  )
}) %>%
  dplyr::bind_rows()

vote_long <- as.data.frame(vote_long, stringsAsFactors = FALSE)

if (nrow(vote_long) == 0) {
  stop("No valid gene records were extracted from all_genes.")
}

vote_long$method <- as.character(vote_long$method)
vote_long$gene <- as.character(vote_long$gene)

vote_long <- vote_long %>%
  dplyr::filter(
    !is.na(method),
    method != "",
    !is.na(gene),
    gene != "",
    gene != "NA",
    gene != "NULL"
  ) %>%
  dplyr::distinct(method, gene)

n_methods_total <- dplyr::n_distinct(vote_long$method)

cat("Total algorithms included:", n_methods_total, "\n")
cat("Total unique genes:", dplyr::n_distinct(vote_long$gene), "\n\n")

write.csv(
  vote_long,
  file = file.path(output_dir, "All_methods_gene_long_table.csv"),
  row.names = FALSE
)


# 6. Gene-level consensus voting table

vote_table <- vote_long %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    n_methods = dplyr::n_distinct(method),
    support_rate = n_methods / n_methods_total,
    methods = paste(sort(unique(method)), collapse = "; "),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_methods), gene)

write.csv(
  vote_table,
  file = file.path(output_dir, "Consensus_gene_vote_table.csv"),
  row.names = FALSE
)


#7. Consensus threshold and consensus genes

min_support_methods <- min(
  n_methods_total,
  max(3, ceiling(0.30 * n_methods_total))
)

consensus_genes <- vote_table %>%
  dplyr::filter(n_methods >= min_support_methods) %>%
  dplyr::pull(gene)

write.table(
  data.frame(Genes = consensus_genes),
  file = file.path(output_dir, "Consensus_Genes_selected_by_voting.txt"),
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  sep = "\t"
)

cat("Consensus threshold: >=", min_support_methods, "algorithms\n")
cat("Consensus genes:\n")
cat(paste(consensus_genes, collapse = ", "), "\n\n")


#  8. Method class annotation

method_info <- vote_long %>%
  dplyr::distinct(method) %>%
  dplyr::mutate(
    method_class = dplyr::case_when(
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
  )

write.csv(
  method_info,
  file = file.path(output_dir, "Feature_selection_method_class_annotation.csv"),
  row.names = FALSE
)


#9. Nature-style theme, palettes and font settings

font_main_title   <- 28
font_panel_tag    <- 28

font_panel_title  <- 28
font_axis_title   <- 28
font_axis_text    <- 20

font_gene_A       <- 20
font_gene_B       <- 20
font_gene_C       <- 20

font_heatmap_x    <- 20
font_facet        <- 20
font_legend_title <- 20
font_legend_text  <- 20

theme_nature <- function(base_size = 15) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      text = ggplot2::element_text(color = "black"),
      
      plot.title = ggplot2::element_text(
        size = font_panel_title,
        face = "bold",
        hjust = 0.5,
        color = "black"
      ),
      
      axis.title = ggplot2::element_text(
        size = font_axis_title,
        face = "bold",
        color = "black"
      ),
      
      axis.text = ggplot2::element_text(
        size = font_axis_text,
        color = "black"
      ),
      
      axis.line = ggplot2::element_line(
        color = "black",
        linewidth = 0.55
      ),
      
      axis.ticks = ggplot2::element_line(
        color = "black",
        linewidth = 0.45
      ),
      
      panel.border = ggplot2::element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.65
      ),
      
      legend.title = ggplot2::element_text(
        size = font_legend_title,
        face = "bold",
        color = "black"
      ),
      
      legend.text = ggplot2::element_text(
        size = font_legend_text,
        color = "black"
      ),
      
      strip.background = ggplot2::element_rect(
        fill = "grey97",
        color = "grey40",
        linewidth = 0.45
      ),
      
      strip.text = ggplot2::element_text(
        size = font_facet,
        face = "bold",
        color = "black"
      ),
      
      plot.margin = ggplot2::margin(10, 10, 10, 10)
    )
}

pal_low    <- "#A6CEE3"
pal_high   <- "#1F78B4"
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


# 10. Panel A: Consensus vote lollipop plot

vote_plot_df <- vote_table %>%
  dplyr::arrange(dplyr::desc(n_methods), gene)

vote_plot_df <- vote_plot_df[
  seq_len(min(top_n_vote, nrow(vote_plot_df))),
  ,
  drop = FALSE
]

vote_plot_df <- vote_plot_df %>%
  dplyr::distinct(gene, .keep_all = TRUE) %>%
  dplyr::arrange(n_methods, gene) %>%
  dplyr::mutate(
    gene = factor(gene, levels = gene)
  )

x_break_by <- max(1, ceiling(n_methods_total / 8))

p_vote <- ggplot2::ggplot(
  vote_plot_df,
  ggplot2::aes(x = n_methods, y = gene)
) +
  ggplot2::geom_segment(
    ggplot2::aes(
      x = 0,
      xend = n_methods,
      y = gene,
      yend = gene,
      color = support_rate
    ),
    linewidth = 1.35,
    lineend = "round"
  ) +
  ggplot2::geom_point(
    ggplot2::aes(fill = support_rate),
    shape = 21,
    color = "black",
    size = 4.0,
    stroke = 0.45
  ) +
  ggplot2::geom_vline(
    xintercept = min_support_methods,
    linetype = "dashed",
    color = pal_cutoff,
    linewidth = 0.75
  ) +
  ggplot2::scale_color_gradient(
    low = pal_low,
    high = pal_high,
    labels = scales::percent_format(accuracy = 1)
  ) +
  ggplot2::scale_fill_gradient(
    low = pal_low,
    high = pal_high,
    labels = scales::percent_format(accuracy = 1)
  ) +
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
  theme_nature(base_size = 13) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = font_axis_text, color = "black"),
    axis.text.y = ggplot2::element_text(size = font_gene_A, color = "black"),
    legend.position = "right",
    legend.key.height = grid::unit(0.70, "cm"),
    legend.key.width  = grid::unit(0.45, "cm"),
    panel.grid.major.x = ggplot2::element_line(
      color = "grey90",
      linewidth = 0.30
    ),
    panel.grid.major.y = ggplot2::element_blank()
  ) +
  ggplot2::guides(
    color = "none",
    fill = ggplot2::guide_colorbar(
      title.position = "top",
      barwidth = grid::unit(0.45, "cm"),
      barheight = grid::unit(3.2, "cm")
    )
  )

p_vote


#11. Algorithm output size and class-size order
#     Key logic:
#     Class with more algorithms is placed further left.

method_count <- vote_long %>%
  dplyr::count(method, name = "n_selected_genes") %>%
  dplyr::left_join(method_info, by = "method")

write.csv(
  method_count,
  file = file.path(output_dir, "Algorithm_selected_gene_count_table.csv"),
  row.names = FALSE
)

class_order <- method_count %>%
  dplyr::group_by(method_class) %>%
  dplyr::summarise(
    n_algorithms = dplyr::n_distinct(method),
    total_selected_genes = sum(n_selected_genes, na.rm = TRUE),
    mean_selected_genes = mean(n_selected_genes, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_algorithms),
    dplyr::desc(total_selected_genes),
    dplyr::desc(mean_selected_genes),
    method_class
  ) %>%
  dplyr::pull(method_class)

method_count <- method_count %>%
  dplyr::mutate(
    method_class = factor(method_class, levels = class_order)
  ) %>%
  dplyr::arrange(
    method_class,
    dplyr::desc(n_selected_genes),
    method
  )

method_order <- method_count$method



# 12. Panel B: Method-gene consensus heatmap


heatmap_genes <- vote_table %>%
  dplyr::arrange(dplyr::desc(n_methods), gene) %>%
  dplyr::pull(gene)

heatmap_genes <- heatmap_genes[
  seq_len(min(heatmap_top_n, length(heatmap_genes)))
]

gene_order <- vote_table %>%
  dplyr::filter(gene %in% heatmap_genes) %>%
  dplyr::arrange(dplyr::desc(n_methods), gene) %>%
  dplyr::pull(gene)

heatmap_long <- expand.grid(
  gene = heatmap_genes,
  method = method_order,
  stringsAsFactors = FALSE
) %>%
  dplyr::left_join(
    vote_long %>% dplyr::mutate(selected = 1),
    by = c("gene", "method")
  ) %>%
  dplyr::mutate(
    selected = ifelse(is.na(selected), 0, selected)
  ) %>%
  dplyr::left_join(method_info, by = "method") %>%
  dplyr::mutate(
    method_class = factor(method_class, levels = class_order),
    tile_fill = ifelse(selected == 1, as.character(method_class), "Not selected"),
    gene = factor(gene, levels = rev(gene_order)),
    method = factor(method, levels = method_order)
  )

p_heatmap <- ggplot2::ggplot(
  heatmap_long,
  ggplot2::aes(x = method, y = gene, fill = tile_fill)
) +
  ggplot2::geom_tile(
    color = "white",
    linewidth = 0.35
  ) +
  ggplot2::facet_grid(
    . ~ method_class,
    scales = "free_x",
    space = "free_x"
  ) +
  ggplot2::scale_fill_manual(
    values = method_class_palette,
    breaks = class_order,
    name = "Algorithm class"
  ) +
  ggplot2::scale_x_discrete(
    labels = function(x) stringr::str_wrap(x, width = 12)
  ) +
  ggplot2::labs(
    x = "Feature-selection algorithms grouped by algorithm class",
    y = NULL,
    title = "Algorithm-level support for candidate genes"
  ) +
  theme_nature(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      size = font_panel_title,
      face = "bold",
      hjust = 0.5,
      color = "black"
    ),
    axis.title.x = ggplot2::element_text(
      size = font_axis_title,
      face = "bold",
      color = "black"
    ),
    axis.text.x = ggplot2::element_text(
      size = font_heatmap_x,
      angle = 50,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = ggplot2::element_text(
      size = font_gene_B,
      color = "black"
    ),
    
    # Remove facet titles and facet-strip background
    strip.text.x = ggplot2::element_blank(),
    strip.background = ggplot2::element_blank(),
    
    # Cleaner heatmap background
    panel.background = ggplot2::element_rect(fill = "white", color = NA),
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.border     = ggplot2::element_blank(),
    panel.grid       = ggplot2::element_blank(),
    
    # Reduce unnecessary visual elements
    axis.line.x  = ggplot2::element_blank(),
    axis.line.y  = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
    
    legend.position = "right",
    legend.title = ggplot2::element_text(
      size = font_legend_title,
      face = "bold",
      color = "black"
    ),
    legend.text = ggplot2::element_text(
      size = font_legend_text,
      color = "black"
    ),
    legend.key.height = grid::unit(0.60, "cm"),
    legend.key.width  = grid::unit(0.60, "cm"),
    
    # Keep only narrow spacing between algorithm-class blocks
    panel.spacing.x = grid::unit(0.08, "lines"),
    
    plot.margin = ggplot2::margin(8, 8, 8, 8)
  )

p_heatmap


# 13. Panel C: Algorithm output-size plot


method_count_plot <- method_count %>%
  dplyr::mutate(
    method = factor(method, levels = rev(method_order))
  )

p_method_count <- ggplot2::ggplot(
  method_count_plot,
  ggplot2::aes(x = n_selected_genes, y = method, fill = method_class)
) +
  ggplot2::geom_col(
    width = 0.70,
    color = "black",
    linewidth = 0.35
  ) +
  ggplot2::scale_fill_manual(
    values = method_class_palette,
    drop = FALSE,
    name = "Algorithm class"
  ) +
  ggplot2::scale_y_discrete(
    labels = function(x) stringr::str_wrap(x, width = 26)
  ) +
  ggplot2::labs(
    x = "Number of selected genes",
    y = NULL,
    title = "Feature output size of individual algorithms"
  ) +
  theme_nature(base_size = 13) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(size = font_axis_text, color = "black"),
    axis.text.y = ggplot2::element_text(size = font_gene_C, color = "black"),
    legend.position = "right",
    legend.key.height = grid::unit(0.60, "cm"),
    legend.key.width  = grid::unit(0.60, "cm"),
    panel.grid.major.x = ggplot2::element_line(
      color = "grey90",
      linewidth = 0.30
    ),
    panel.grid.major.y = ggplot2::element_blank()
  )

p_method_count



# 14. Dynamic figure size


n_methods_total <- dplyr::n_distinct(vote_long$method)
n_heatmap_genes <- length(heatmap_genes)

vote_height <- max(5.5, 0.30 * nrow(vote_plot_df) + 1.8)

heatmap_width <- max(24, min(46, 0.90 * n_methods_total + 12))
heatmap_height <- max(7.5, 0.30 * n_heatmap_genes + 3.0)

method_count_width <- 11
method_count_height <- max(6.5, 0.30 * nrow(method_count_plot) + 2.0)

combined_width <- max(24, min(48, 0.92 * n_methods_total + 13))
combined_height <- max(22, 0.38 * n_heatmap_genes + 12)



# 15. Save individual panels


ggplot2::ggsave(
  filename = file.path(output_dir, "NatureStyle_Consensus_gene_vote_lollipop_largefont.pdf"),
  plot = p_vote,
  width = 10,
  height = vote_height,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "NatureStyle_Consensus_gene_vote_lollipop_largefont.png"),
  plot = p_vote,
  width = 10,
  height = vote_height,
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "NatureStyle_Method_gene_consensus_heatmap_class_size_ordered_largefont.pdf"),
  plot = p_heatmap,
  width = heatmap_width,
  height = heatmap_height,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "NatureStyle_Method_gene_consensus_heatmap_class_size_ordered_largefont.png"),
  plot = p_heatmap,
  width = heatmap_width,
  height = heatmap_height,
  dpi = 600,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "NatureStyle_Algorithm_selected_gene_count_class_size_ordered_largefont.pdf"),
  plot = p_method_count,
  width = method_count_width,
  height = method_count_height,
  bg = "white"
)

ggplot2::ggsave(
  filename = file.path(output_dir, "NatureStyle_Algorithm_selected_gene_count_class_size_ordered_largefont.png"),
  plot = p_method_count,
  width = method_count_width,
  height = method_count_height,
  dpi = 600,
  bg = "white"
)



# 16. Final combined figure


p_combined_final <- p_vote / p_heatmap / p_method_count +
  patchwork::plot_layout(
    heights = c(1.00, 1.90, 1.5)
  ) +
  patchwork::plot_annotation(
    tag_levels = "A",
    title = "Multi-algorithm consensus feature-selection overview",
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(
        color = "black",
        size = font_main_title,
        face = "bold",
        hjust = 0.5
      ),
      plot.tag = ggplot2::element_text(
        color = "black",
        size = font_panel_tag,
        face = "bold"
      )
    )
  )

p_combined_final

ggplot2::ggsave(
  filename = file.path(output_dir, "NatureStyle_Consensus_feature_selection_overview_FINAL_largefont.pdf"),
  plot = p_combined_final,
  width = combined_width+12,
  height = combined_height+20,
  bg = "white"
)



cat("Final large-font Nature-style consensus figures have been saved in:\n")
cat(output_dir, "\n")
cat("Main output file:\n")
cat(file.path(output_dir, "NatureStyle_Consensus_feature_selection_overview_FINAL_largefont.pdf"), "\n")









}


