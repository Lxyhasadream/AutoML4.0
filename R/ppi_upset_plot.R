# Publication-style UpSet plot for PPI hub algorithm intersections.

auto_ppi_upset_plot <- function(ppi_result,
                                methods = ppi_methods(),
                                top_n = 10,
                                output_dir = NULL,
                                file_prefix = NULL,
                                max_intersections = 30,
                                min_intersection_size = 1,
                                width = 11,
                                height = 8,
                                dpi = 600,
                                title = NULL,
                                subtitle = NULL,
                                write_outputs = TRUE) {
  ppi_need_package("ggplot2")
  ppi_need_package("patchwork")

  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }
  max_intersections <- as.integer(max_intersections)
  if (is.na(max_intersections) || max_intersections < 1) {
    stop("max_intersections must be a positive integer.", call. = FALSE)
  }
  min_intersection_size <- as.integer(min_intersection_size)
  if (is.na(min_intersection_size) || min_intersection_size < 1) {
    stop("min_intersection_size must be a positive integer.", call. = FALSE)
  }

  ppi_obj <- ppi_extract_analysis_result(ppi_result)
  selected_sets <- ppi_extract_top_sets(ppi_obj, methods = methods, top_n = top_n)
  upset_data <- ppi_make_upset_data(
    selected_sets = selected_sets,
    max_intersections = max_intersections,
    min_intersection_size = min_intersection_size
  )

  if (is.null(title)) {
    title <- "PPI Hub Algorithm Intersection UpSet Plot"
  }
  if (is.null(subtitle)) {
    subtitle <- paste0("Top ", top_n, " genes from ", length(selected_sets), " selected algorithm(s)")
  }

  bar_plot <- ppi_upset_bar_plot(upset_data$intersection_summary, title, subtitle)
  matrix_plot <- ppi_upset_matrix_plot(upset_data$matrix_df)
  plot <- bar_plot / matrix_plot + patchwork::plot_layout(heights = c(0.55, 0.45))

  if (write_outputs) {
    if (is.null(output_dir)) {
      output_dir <- ppi_default_output_dir(ppi_result, "PPI_upset_results")
    }
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
    if (is.null(file_prefix)) {
      file_prefix <- paste0("PPI_upset_top", top_n, "_", length(selected_sets), "_methods")
    }
    utils::write.csv(
      upset_data$intersection_summary,
      file.path(output_dir, paste0(file_prefix, "_intersection_summary.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      upset_data$gene_membership,
      file.path(output_dir, paste0(file_prefix, "_gene_membership.csv")),
      row.names = FALSE
    )
    ggplot2::ggsave(
      file.path(output_dir, paste0(file_prefix, ".pdf")),
      plot,
      width = width,
      height = height,
      bg = "white",
      device = grDevices::cairo_pdf
    )
    ggplot2::ggsave(
      file.path(output_dir, paste0(file_prefix, ".png")),
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }

  structure(
    list(
      plot = plot,
      selected_sets = selected_sets,
      intersection_summary = upset_data$intersection_summary,
      gene_membership = upset_data$gene_membership,
      matrix_data = upset_data$matrix_df,
      output_dir = if (write_outputs) normalizePath(output_dir, mustWork = FALSE) else NA_character_
    ),
    class = "auto_ppi_upset_plot"
  )
}

ppi_extract_analysis_result <- function(ppi_result) {
  if (inherits(ppi_result, "auto_string_ppi_analysis")) {
    return(ppi_result$ppi)
  }
  if (inherits(ppi_result, "auto_ppi_analysis")) {
    return(ppi_result)
  }
  if (is.list(ppi_result) && !is.null(ppi_result$ppi) && inherits(ppi_result$ppi, "auto_ppi_analysis")) {
    return(ppi_result$ppi)
  }
  if (is.list(ppi_result) && (!is.null(ppi_result$scores) || !is.null(ppi_result$top_nodes))) {
    return(ppi_result)
  }
  stop("ppi_result must be an auto_ppi_analysis or auto_string_ppi_analysis result.", call. = FALSE)
}

ppi_extract_top_sets <- function(ppi_obj, methods, top_n) {
  methods <- unique(as.character(methods))
  methods <- methods[!is.na(methods) & nzchar(methods)]
  if (length(methods) < 2) {
    stop("At least two methods are required for an UpSet plot.", call. = FALSE)
  }

  if (!is.null(ppi_obj$scores)) {
    missing_methods <- setdiff(methods, colnames(ppi_obj$scores))
    if (length(missing_methods) > 0) {
      stop("Method(s) not found in PPI score table: ", paste(missing_methods, collapse = ", "), call. = FALSE)
    }
    sets <- lapply(methods, function(method) {
      score_df <- ppi_obj$scores
      ordered_index <- order(-score_df[[method]], score_df$Node)
      unique(as.character(utils::head(score_df$Node[ordered_index], top_n)))
    })
    names(sets) <- methods
    return(sets)
  }

  if (!is.null(ppi_obj$top_nodes) && is.list(ppi_obj$top_nodes)) {
    missing_methods <- setdiff(methods, names(ppi_obj$top_nodes))
    if (length(missing_methods) > 0) {
      stop("Method(s) not found in PPI top-node lists: ", paste(missing_methods, collapse = ", "), call. = FALSE)
    }
    sets <- lapply(methods, function(method) {
      unique(as.character(utils::head(ppi_obj$top_nodes[[method]], top_n)))
    })
    names(sets) <- methods
    return(sets)
  }

  stop("ppi_result must contain either scores or top_nodes.", call. = FALSE)
}

ppi_make_upset_data <- function(selected_sets, max_intersections, min_intersection_size) {
  all_genes <- sort(unique(unlist(selected_sets, use.names = FALSE)))
  if (length(all_genes) == 0) {
    stop("No genes were found in the selected top lists.", call. = FALSE)
  }

  membership <- data.frame(Gene = all_genes, stringsAsFactors = FALSE)
  for (method in names(selected_sets)) {
    membership[[method]] <- all_genes %in% selected_sets[[method]]
  }

  method_cols <- names(selected_sets)
  membership$Combination_ID <- apply(membership[, method_cols, drop = FALSE], 1, function(x) {
    present <- method_cols[as.logical(x)]
    if (length(present) == 0) {
      return(NA_character_)
    }
    paste(present, collapse = " + ")
  })
  membership$Method_Count <- rowSums(membership[, method_cols, drop = FALSE])
  membership <- membership[!is.na(membership$Combination_ID), , drop = FALSE]

  combo_rows <- split(membership, membership$Combination_ID)
  intersection_summary <- do.call(rbind, lapply(names(combo_rows), function(combo) {
    df <- combo_rows[[combo]]
    present_methods <- strsplit(combo, " + ", fixed = TRUE)[[1]]
    data.frame(
      Combination_ID = combo,
      Method_Count = length(present_methods),
      Intersection_Size = nrow(df),
      Genes = paste(sort(df$Gene), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }))
  intersection_summary <- intersection_summary[intersection_summary$Intersection_Size >= min_intersection_size, , drop = FALSE]
  intersection_summary <- intersection_summary[
    order(-intersection_summary$Intersection_Size, -intersection_summary$Method_Count, intersection_summary$Combination_ID),
    ,
    drop = FALSE
  ]
  intersection_summary <- utils::head(intersection_summary, max_intersections)
  if (nrow(intersection_summary) == 0) {
    stop("No intersections meet min_intersection_size.", call. = FALSE)
  }
  intersection_summary$Intersection <- paste0("I", seq_len(nrow(intersection_summary)))

  kept_combos <- intersection_summary$Combination_ID
  membership <- membership[membership$Combination_ID %in% kept_combos, , drop = FALSE]

  matrix_df <- do.call(rbind, lapply(seq_len(nrow(intersection_summary)), function(i) {
    combo <- intersection_summary$Combination_ID[i]
    present <- strsplit(combo, " + ", fixed = TRUE)[[1]]
    data.frame(
      Intersection = intersection_summary$Intersection[i],
      Combination_ID = combo,
      Method = method_cols,
      Present = method_cols %in% present,
      stringsAsFactors = FALSE
    )
  }))

  set_sizes <- data.frame(
    Method = names(selected_sets),
    Set_Size = as.integer(vapply(selected_sets, length, integer(1))),
    stringsAsFactors = FALSE
  )
  set_sizes$Method_Label <- paste0(set_sizes$Method, " (n=", set_sizes$Set_Size, ")")
  rownames(set_sizes) <- set_sizes$Method

  matrix_df$Intersection <- factor(matrix_df$Intersection, levels = intersection_summary$Intersection)
  matrix_df$Method <- factor(matrix_df$Method, levels = rev(method_cols))
  matrix_df$Method_Label <- set_sizes[as.character(matrix_df$Method), "Method_Label"]
  matrix_df$Method_Label <- factor(matrix_df$Method_Label, levels = rev(set_sizes[method_cols, "Method_Label"]))
  intersection_summary$Intersection <- factor(intersection_summary$Intersection, levels = intersection_summary$Intersection)

  list(
    intersection_summary = intersection_summary,
    gene_membership = membership,
    matrix_df = matrix_df,
    set_sizes = set_sizes
  )
}

ppi_upset_bar_plot <- function(intersection_summary, title, subtitle) {
  ggplot2::ggplot(intersection_summary, ggplot2::aes(x = Intersection, y = Intersection_Size)) +
    ggplot2::geom_col(width = 0.72, fill = "#2A6F97", color = "black", linewidth = 0.25) +
    ggplot2::geom_text(
      ggplot2::aes(label = Intersection_Size),
      vjust = -0.35,
      size = 4.2,
      family = "Arial",
      fontface = "bold",
      color = "black"
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.16))) +
    ggplot2::labs(x = NULL, y = "Intersection size", title = title, subtitle = subtitle) +
    ppi_upset_theme() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey90", linewidth = 0.3)
    )
}

ppi_upset_matrix_plot <- function(matrix_df) {
  active_df <- matrix_df[matrix_df$Present, , drop = FALSE]
  connector_df <- do.call(rbind, lapply(split(active_df, active_df$Intersection), function(df) {
    if (nrow(df) < 2) {
      return(NULL)
    }
    y_values <- as.numeric(df$Method_Label)
    label_levels <- levels(df$Method_Label)
    data.frame(
      Intersection = df$Intersection[1],
      y_min = label_levels[min(y_values)],
      y_max = label_levels[max(y_values)],
      stringsAsFactors = FALSE
    )
  }))

  p <- ggplot2::ggplot(matrix_df, ggplot2::aes(x = Intersection, y = Method_Label))
  if (!is.null(connector_df) && nrow(connector_df) > 0) {
    p <- p + ggplot2::geom_segment(
      data = connector_df,
      ggplot2::aes(x = Intersection, xend = Intersection, y = y_min, yend = y_max),
      inherit.aes = FALSE,
      linewidth = 1.0,
      color = "#334E68",
      lineend = "round"
    )
  }

  p +
    ggplot2::geom_point(
      data = matrix_df[!matrix_df$Present, , drop = FALSE],
      shape = 21,
      size = 4.2,
      fill = "#E7EEF5",
      color = "#CAD6E0",
      stroke = 0.25
    ) +
    ggplot2::geom_point(
      data = active_df,
      shape = 21,
      size = 5.0,
      fill = "#F28E2B",
      color = "black",
      stroke = 0.35
    ) +
    ggplot2::labs(x = "Intersection", y = NULL) +
    ppi_upset_theme() +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "grey92", linewidth = 0.3),
      axis.text.x = ggplot2::element_text(family = "Arial", size = 14, color = "black", face = "bold"),
      axis.text.y = ggplot2::element_text(family = "Arial", size = 14, color = "black")
    )
}

ppi_upset_theme <- function() {
  ggplot2::theme_classic(base_size = 14, base_family = "Arial") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial", color = "black"),
      plot.title = ggplot2::element_text(family = "Arial", size = 22, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = ggplot2::element_text(family = "Arial", size = 14, hjust = 0.5, color = "#334E68"),
      axis.title = ggplot2::element_text(family = "Arial", size = 16, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(family = "Arial", size = 14, color = "black"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.45),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.35),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

ppi_default_output_dir <- function(ppi_result, fallback) {
  if (is.list(ppi_result) && !is.null(ppi_result$output_dir) && !is.na(ppi_result$output_dir)) {
    return(ppi_result$output_dir)
  }
  if (is.list(ppi_result) && !is.null(ppi_result$ppi) && !is.null(ppi_result$ppi$output_dir) && !is.na(ppi_result$ppi$output_dir)) {
    return(ppi_result$ppi$output_dir)
  }
  fallback
}
