# cytoHubba-style PPI hub screening with 11 topology algorithms.

ppi_methods <- function() {
  c(
    "Degree",
    "MNC",
    "DMNC",
    "MCC",
    "BottleNeck",
    "EcCentricity",
    "Closeness",
    "Radiality",
    "Betweenness",
    "Stress",
    "EPC"
  )
}

auto_ppi_analysis <- function(edge_df,
                              from_col = 1,
                              to_col = 2,
                              output_dir = "PPI_screening_results",
                              top_n = 10,
                              consensus_min_methods = NULL,
                              epsilon = 1.7,
                              epc_threshold = 0.1,
                              epc_n_sim = 1000,
                              seed = 123,
                              include_self_radiality = TRUE,
                              write_outputs = TRUE,
                              plot_outputs = TRUE,
                              verbose = TRUE) {
  ppi_need_package("igraph")

  top_n <- as.integer(top_n)
  if (is.na(top_n) || top_n < 1) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }
  if (!is.null(consensus_min_methods)) {
    consensus_min_methods <- as.integer(consensus_min_methods)
    if (is.na(consensus_min_methods) || consensus_min_methods < 1) {
      stop("consensus_min_methods must be a positive integer.", call. = FALSE)
    }
  }

  graph <- ppi_make_graph(edge_df, from_col = from_col, to_col = to_col)
  if (igraph::vcount(graph) == 0) {
    stop("The PPI network has no valid nodes after preprocessing.", call. = FALSE)
  }

  if (verbose) {
    message("PPI network: ", igraph::vcount(graph), " nodes, ", igraph::ecount(graph), " edges, ",
            igraph::components(graph)$no, " connected component(s).")
  }

  local_scores <- ppi_calc_mnc_dmnc(graph, epsilon = epsilon)
  shortest_scores <- ppi_calc_shortest_path_scores(graph, include_self_radiality = include_self_radiality)
  degree_score <- igraph::degree(graph, mode = "all")
  names(degree_score) <- igraph::V(graph)$name

  score_df <- data.frame(
    Node = igraph::V(graph)$name,
    Degree = as.numeric(degree_score[igraph::V(graph)$name]),
    MNC = as.numeric(local_scores$MNC[igraph::V(graph)$name]),
    DMNC = as.numeric(local_scores$DMNC[igraph::V(graph)$name]),
    MCC = as.numeric(ppi_calc_mcc(graph)[igraph::V(graph)$name]),
    BottleNeck = as.numeric(ppi_calc_bottleneck(graph)[igraph::V(graph)$name]),
    EcCentricity = as.numeric(shortest_scores$EcCentricity[igraph::V(graph)$name]),
    Closeness = as.numeric(shortest_scores$Closeness[igraph::V(graph)$name]),
    Radiality = as.numeric(shortest_scores$Radiality[igraph::V(graph)$name]),
    Betweenness = as.numeric(shortest_scores$Betweenness[igraph::V(graph)$name]),
    Stress = as.numeric(ppi_calc_stress(graph)[igraph::V(graph)$name]),
    EPC = as.numeric(ppi_calc_epc(
      graph,
      threshold = epc_threshold,
      n_sim = epc_n_sim,
      seed = seed,
      average = TRUE
    )[igraph::V(graph)$name]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  top_tables <- ppi_make_top_tables(score_df, top_n = top_n)
  top_nodes <- lapply(top_tables, function(x) x$Node)
  top_combined <- ppi_bind_top_tables(top_tables)
  hub_frequency <- ppi_make_hub_frequency(top_nodes)
  integrated_rank <- ppi_make_integrated_rank(score_df)

  if (is.null(consensus_min_methods)) {
    consensus_min_methods <- min(length(ppi_methods()), max(3, ceiling(0.30 * length(ppi_methods()))))
  }
  consensus_hubs <- hub_frequency[hub_frequency$Frequency_in_top_lists >= consensus_min_methods, , drop = FALSE]
  intersections <- ppi_calculate_algorithm_intersections(top_nodes)

  if (write_outputs || plot_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }
  if (write_outputs) {
    ppi_write_outputs(
      score_df = score_df,
      top_combined = top_combined,
      hub_frequency = hub_frequency,
      integrated_rank = integrated_rank,
      consensus_hubs = consensus_hubs,
      intersections = intersections,
      output_dir = output_dir,
      top_n = top_n,
      consensus_min_methods = consensus_min_methods
    )
  }
  if (plot_outputs) {
    ppi_plot_outputs(
      graph = graph,
      score_df = score_df,
      hub_frequency = hub_frequency,
      consensus_hubs = consensus_hubs,
      output_dir = output_dir
    )
  }

  structure(
    list(
      graph = graph,
      scores = score_df,
      top_tables = top_tables,
      top_nodes = top_nodes,
      top_combined = top_combined,
      hub_frequency = hub_frequency,
      integrated_rank = integrated_rank,
      consensus_hubs = consensus_hubs,
      intersections = intersections,
      method_formula_check = ppi_formula_check_table(),
      output_dir = if (write_outputs || plot_outputs) normalizePath(output_dir, mustWork = FALSE) else NA_character_
    ),
    class = "auto_ppi_analysis"
  )
}

ppi_make_graph <- function(edge_df, from_col = 1, to_col = 2) {
  get_col <- function(df, col) {
    if (is.numeric(col)) {
      return(df[[col]])
    }
    if (!col %in% colnames(df)) {
      stop("Column not found: ", col, call. = FALSE)
    }
    df[[col]]
  }

  edges <- data.frame(
    from = as.character(get_col(edge_df, from_col)),
    to = as.character(get_col(edge_df, to_col)),
    stringsAsFactors = FALSE
  )
  edges <- edges[!is.na(edges$from) & !is.na(edges$to), , drop = FALSE]
  edges <- edges[edges$from != "" & edges$to != "", , drop = FALSE]
  edges <- edges[edges$from != edges$to, , drop = FALSE]
  if (nrow(edges) == 0) {
    stop("No valid PPI edges were found.", call. = FALSE)
  }

  graph <- igraph::graph_from_data_frame(edges, directed = FALSE)
  igraph::simplify(graph, remove.multiple = TRUE, remove.loops = TRUE)
}

ppi_calc_mnc_dmnc <- function(graph, epsilon = 1.7) {
  n <- igraph::vcount(graph)
  mnc <- numeric(n)
  dmnc <- numeric(n)

  for (i in seq_len(n)) {
    nb <- igraph::neighbors(graph, i, mode = "all")
    if (length(nb) == 0) {
      next
    }
    subg <- igraph::induced_subgraph(graph, nb)
    comp <- igraph::components(subg)
    mc_nodes <- which(comp$membership == which.max(comp$csize))
    mc_subg <- igraph::induced_subgraph(subg, mc_nodes)
    nv <- igraph::vcount(mc_subg)
    ne <- igraph::ecount(mc_subg)
    mnc[i] <- nv
    dmnc[i] <- ifelse(nv > 0, ne / (nv ^ epsilon), 0)
  }

  names(mnc) <- igraph::V(graph)$name
  names(dmnc) <- igraph::V(graph)$name
  list(MNC = mnc, DMNC = dmnc)
}

ppi_calc_mcc <- function(graph) {
  n <- igraph::vcount(graph)
  mcc <- numeric(n)
  cliques <- igraph::max_cliques(graph, min = 2)
  if (length(cliques) > 0) {
    for (clique in cliques) {
      mcc[as.integer(clique)] <- mcc[as.integer(clique)] + factorial(length(clique) - 1)
    }
  }
  names(mcc) <- igraph::V(graph)$name
  mcc
}

ppi_calc_shortest_path_scores <- function(graph, include_self_radiality = TRUE) {
  n <- igraph::vcount(graph)
  dmat <- as.matrix(igraph::distances(graph, v = igraph::V(graph), to = igraph::V(graph), mode = "all", weights = NA))

  closeness_score <- rowSums(ifelse(is.finite(dmat) & dmat > 0, 1 / dmat, 0), na.rm = TRUE)
  comp <- igraph::components(graph)
  membership <- comp$membership
  comp_ids <- unique(membership)
  eccentricity_score <- numeric(n)
  radiality_score <- numeric(n)

  component_diameter <- stats::setNames(numeric(length(comp_ids)), as.character(comp_ids))
  for (cid in comp_ids) {
    idx <- which(membership == cid)
    subd <- dmat[idx, idx, drop = FALSE]
    finite_d <- subd[is.finite(subd)]
    component_diameter[as.character(cid)] <- ifelse(length(finite_d) > 0, max(finite_d), 0)
  }

  for (i in seq_len(n)) {
    cid <- membership[i]
    idx <- which(membership == cid)
    comp_size <- length(idx)
    d <- dmat[i, idx]
    d <- d[is.finite(d)]
    max_dist_from_i <- max(d)
    eccentricity_score[i] <- ifelse(max_dist_from_i > 0, (comp_size / n) * (1 / max_dist_from_i), 0)

    delta_c <- component_diameter[as.character(cid)]
    if (delta_c > 0) {
      d_for_radiality <- if (include_self_radiality) d else d[d > 0]
      radiality_score[i] <- (comp_size / n) * sum((delta_c + 1 - d_for_radiality) / delta_c)
    }
  }

  betweenness_score <- igraph::betweenness(graph, v = igraph::V(graph), directed = FALSE, weights = NA, normalized = FALSE)
  names(closeness_score) <- igraph::V(graph)$name
  names(eccentricity_score) <- igraph::V(graph)$name
  names(radiality_score) <- igraph::V(graph)$name
  names(betweenness_score) <- igraph::V(graph)$name
  list(Closeness = closeness_score, EcCentricity = eccentricity_score, Radiality = radiality_score, Betweenness = betweenness_score)
}

ppi_calc_stress <- function(graph) {
  n <- igraph::vcount(graph)
  adj <- igraph::as_adj_list(graph, mode = "all")
  stress_score <- numeric(n)

  for (source in seq_len(n)) {
    stack <- integer(0)
    pred <- vector("list", n)
    for (i in seq_len(n)) pred[[i]] <- integer(0)
    sigma <- numeric(n)
    sigma[source] <- 1
    dist <- rep(-1L, n)
    dist[source] <- 0L
    queue <- integer(n)
    head <- 1L
    tail <- 1L
    queue[tail] <- source

    while (head <= tail) {
      v <- queue[head]
      head <- head + 1L
      stack <- c(stack, v)
      for (w in as.integer(adj[[v]])) {
        if (dist[w] < 0L) {
          tail <- tail + 1L
          queue[tail] <- w
          dist[w] <- dist[v] + 1L
        }
        if (dist[w] == dist[v] + 1L) {
          sigma[w] <- sigma[w] + sigma[v]
          pred[[w]] <- c(pred[[w]], v)
        }
      }
    }

    dependency <- numeric(n)
    for (w in rev(stack)) {
      for (v in pred[[w]]) {
        if (sigma[w] > 0) {
          dependency[v] <- dependency[v] + (sigma[v] / sigma[w]) * (sigma[w] + dependency[w])
        }
      }
      if (w != source) {
        stress_score[w] <- stress_score[w] + dependency[w]
      }
    }
  }

  stress_score <- stress_score / 2
  names(stress_score) <- igraph::V(graph)$name
  stress_score
}

ppi_calc_bottleneck <- function(graph) {
  n <- igraph::vcount(graph)
  adj <- igraph::as_adj_list(graph, mode = "all")
  bottleneck_score <- numeric(n)

  for (source in seq_len(n)) {
    parent <- rep(NA_integer_, n)
    dist <- rep(-1L, n)
    dist[source] <- 0L
    queue <- integer(n)
    head <- 1L
    tail <- 1L
    queue[tail] <- source

    while (head <= tail) {
      v <- queue[head]
      head <- head + 1L
      for (w in as.integer(adj[[v]])) {
        if (dist[w] < 0L) {
          tail <- tail + 1L
          queue[tail] <- w
          dist[w] <- dist[v] + 1L
          parent[w] <- v
        }
      }
    }

    reachable <- which(dist >= 0L)
    tree_size <- length(reachable)
    meet_count <- numeric(n)
    for (target in reachable) {
      if (target == source) next
      x <- parent[target]
      while (!is.na(x) && x != source) {
        meet_count[x] <- meet_count[x] + 1
        x <- parent[x]
      }
    }
    bottleneck_score <- bottleneck_score + as.numeric(meet_count > tree_size / 4)
  }

  names(bottleneck_score) <- igraph::V(graph)$name
  bottleneck_score
}

ppi_calc_epc <- function(graph, threshold = 0.1, n_sim = 1000, seed = 123, average = TRUE) {
  set.seed(seed)
  n <- igraph::vcount(graph)
  m <- igraph::ecount(graph)
  epc_score <- numeric(n)
  if (m == 0) {
    epc_score <- rep(1 / n, n)
    names(epc_score) <- igraph::V(graph)$name
    return(epc_score)
  }

  for (k in seq_len(n_sim)) {
    keep_edges <- stats::runif(m) >= threshold
    if (utils::packageVersion("igraph") >= "2.1.0") {
      gk <- igraph::subgraph_from_edges(graph = graph, eids = igraph::E(graph)[keep_edges], delete.vertices = FALSE)
    } else {
      gk <- igraph::subgraph.edges(graph = graph, eids = igraph::E(graph)[keep_edges], delete.vertices = FALSE)
    }
    comp <- igraph::components(gk)
    epc_score <- epc_score + comp$csize[comp$membership] / n
  }
  if (average) {
    epc_score <- epc_score / n_sim
  }
  names(epc_score) <- igraph::V(graph)$name
  epc_score
}

ppi_make_top_tables <- function(score_df, top_n) {
  out <- list()
  for (method in ppi_methods()) {
    ordered_index <- order(-score_df[[method]], score_df$Node)
    tmp <- score_df[ordered_index, c("Node", method), drop = FALSE]
    tmp <- head(tmp, top_n)
    tmp$Rank <- seq_len(nrow(tmp))
    out[[method]] <- tmp[, c("Rank", "Node", method), drop = FALSE]
  }
  out
}

ppi_bind_top_tables <- function(top_tables) {
  do.call(rbind, lapply(names(top_tables), function(method) {
    tmp <- top_tables[[method]]
    data.frame(Method = method, Rank = tmp$Rank, Node = tmp$Node, Score = tmp[[method]], stringsAsFactors = FALSE)
  }))
}

ppi_make_hub_frequency <- function(top_nodes) {
  freq <- sort(table(unlist(top_nodes)), decreasing = TRUE)
  data.frame(Node = names(freq), Frequency_in_top_lists = as.integer(freq), stringsAsFactors = FALSE)
}

ppi_make_integrated_rank <- function(score_df) {
  rank_df <- score_df
  for (method in ppi_methods()) {
    rank_df[[paste0(method, "_Rank")]] <- rank(-rank_df[[method]], ties.method = "min")
  }
  rank_cols <- paste0(ppi_methods(), "_Rank")
  rank_df$Mean_Rank <- rowMeans(rank_df[, rank_cols, drop = FALSE], na.rm = TRUE)
  rank_df$Median_Rank <- apply(rank_df[, rank_cols, drop = FALSE], 1, stats::median, na.rm = TRUE)
  rank_df[order(rank_df$Mean_Rank, rank_df$Median_Rank, rank_df$Node), , drop = FALSE]
}

ppi_calculate_algorithm_intersections <- function(top_nodes, keep_empty = TRUE) {
  method_names <- names(top_nodes)
  if (length(method_names) < 2) {
    stop("At least two algorithms are required for intersection analysis.", call. = FALSE)
  }
  intersection_list <- list()
  gene_long_list <- list()
  idx <- 1L

  for (k in 2:length(method_names)) {
    for (comb in utils::combn(method_names, k, simplify = FALSE)) {
      common_genes <- sort(unique(Reduce(intersect, top_nodes[comb])))
      if (!keep_empty && length(common_genes) == 0) {
        next
      }
      combination_id <- paste(comb, collapse = " + ")
      intersection_list[[idx]] <- data.frame(
        Algorithm_Number = k,
        Combination_ID = combination_id,
        Common_Gene_Number = length(common_genes),
        Common_Genes = ifelse(length(common_genes) > 0, paste(common_genes, collapse = ";"), NA),
        stringsAsFactors = FALSE
      )
      if (length(common_genes) > 0) {
        gene_long_list[[idx]] <- data.frame(Algorithm_Number = k, Combination_ID = combination_id, Gene = common_genes, stringsAsFactors = FALSE)
      }
      idx <- idx + 1L
    }
  }

  intersection_summary <- do.call(rbind, intersection_list)
  intersection_summary <- intersection_summary[order(-intersection_summary$Algorithm_Number, -intersection_summary$Common_Gene_Number, intersection_summary$Combination_ID), , drop = FALSE]
  rownames(intersection_summary) <- NULL
  intersection_gene_long <- if (length(gene_long_list) > 0) do.call(rbind, gene_long_list) else data.frame(Algorithm_Number = integer(), Combination_ID = character(), Gene = character())
  list(intersection_summary = intersection_summary, intersection_gene_long = intersection_gene_long)
}

ppi_write_outputs <- function(score_df,
                              top_combined,
                              hub_frequency,
                              integrated_rank,
                              consensus_hubs,
                              intersections,
                              output_dir,
                              top_n,
                              consensus_min_methods) {
  utils::write.csv(score_df, file.path(output_dir, "PPI_cytohubba_11_scores.csv"), row.names = FALSE)
  utils::write.csv(top_combined, file.path(output_dir, paste0("PPI_top", top_n, "_by_each_method.csv")), row.names = FALSE)
  utils::write.csv(hub_frequency, file.path(output_dir, paste0("PPI_top", top_n, "_hub_frequency.csv")), row.names = FALSE)
  utils::write.csv(integrated_rank, file.path(output_dir, "PPI_all_nodes_11_methods_with_integrated_rank.csv"), row.names = FALSE)
  utils::write.csv(consensus_hubs, file.path(output_dir, paste0("PPI_consensus_hubs_frequency_ge_", consensus_min_methods, ".csv")), row.names = FALSE)
  utils::write.csv(intersections$intersection_summary, file.path(output_dir, paste0("PPI_top", top_n, "_algorithm_intersection_summary_2_to_11.csv")), row.names = FALSE)
  utils::write.csv(intersections$intersection_gene_long, file.path(output_dir, paste0("PPI_top", top_n, "_algorithm_intersection_gene_long_2_to_11.csv")), row.names = FALSE)
}

ppi_plot_outputs <- function(graph, score_df, hub_frequency, consensus_hubs, output_dir) {
  ppi_need_package("ggplot2")
  figdir <- file.path(output_dir, "figures")
  dir.create(figdir, showWarnings = FALSE, recursive = TRUE)

  p_freq <- ggplot2::ggplot(hub_frequency, ggplot2::aes(x = Frequency_in_top_lists, y = factor(Node, levels = rev(Node)))) +
    ggplot2::geom_col(width = 0.75, fill = "#1F78B4", color = "black", linewidth = 0.25) +
    ggplot2::scale_x_continuous(breaks = seq(0, length(ppi_methods()), 1), limits = c(0, length(ppi_methods())), expand = ggplot2::expansion(mult = c(0, 0.03))) +
    ggplot2::labs(x = "Frequency across top-ranked lists", y = NULL, title = "Consensus Hub Frequency Across 11 PPI Algorithms") +
    ppi_plot_theme()
  ppi_save_plot(p_freq, file.path(figdir, "PPI_hub_frequency_barplot"), 8, max(5, 0.32 * nrow(hub_frequency) + 2))

  heatmap_genes <- as.character(hub_frequency$Node)
  heatmap_df <- score_df[match(heatmap_genes, score_df$Node), c("Node", ppi_methods()), drop = FALSE]
  heat_long <- do.call(rbind, lapply(ppi_methods(), function(method) {
    x <- heatmap_df[[method]]
    scaled <- if (max(x, na.rm = TRUE) == min(x, na.rm = TRUE)) rep(0, length(x)) else (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
    data.frame(Node = heatmap_df$Node, Method = method, Score = scaled, stringsAsFactors = FALSE)
  }))
  heat_long$Node <- factor(heat_long$Node, levels = rev(heatmap_genes))
  heat_long$Method <- factor(heat_long$Method, levels = ppi_methods())
  p_heat <- ggplot2::ggplot(heat_long, ggplot2::aes(x = Method, y = Node, fill = Score)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.3) +
    ggplot2::scale_fill_gradient(low = "#F2F6FB", high = "#B2182B", name = "Scaled score") +
    ggplot2::labs(x = "PPI topology algorithm", y = NULL, title = "Normalized Topological Scores Across 11 Algorithms") +
    ppi_plot_theme() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(family = "Arial", size = 18, angle = 45, hjust = 1))
  ppi_save_plot(p_heat, file.path(figdir, "PPI_normalized_11_algorithm_score_heatmap"), max(9, 0.72 * length(ppi_methods())), max(6, 0.32 * length(heatmap_genes) + 2))

  ppi_plot_network(graph, consensus_hubs$Node, file.path(figdir, "PPI_network_consensus_hubs"))
}

ppi_plot_theme <- function() {
  ggplot2::theme_classic(base_size = 14, base_family = "Arial") +
    ggplot2::theme(
      text = ggplot2::element_text(family = "Arial", color = "black"),
      plot.title = ggplot2::element_text(family = "Arial", size = 24, face = "bold", hjust = 0.5, color = "black"),
      plot.subtitle = ggplot2::element_text(family = "Arial", size = 24, color = "black"),
      axis.title = ggplot2::element_text(family = "Arial", size = 20, face = "bold", color = "black"),
      axis.text = ggplot2::element_text(family = "Arial", size = 18, color = "black"),
      legend.title = ggplot2::element_text(family = "Arial", size = 18, face = "bold", color = "black"),
      legend.text = ggplot2::element_text(family = "Arial", size = 18, color = "black")
    )
}

ppi_save_plot <- function(plot, file_base, width, height, dpi = 600) {
  ggplot2::ggsave(paste0(file_base, ".pdf"), plot, width = width, height = height, bg = "white", device = grDevices::cairo_pdf)
  ggplot2::ggsave(paste0(file_base, ".png"), plot, width = width, height = height, dpi = dpi, bg = "white")
}

ppi_plot_network <- function(graph, consensus_nodes, file_base) {
  node_degree <- igraph::degree(graph)
  node_size <- 6 + 10 * (node_degree - min(node_degree)) / max(1, max(node_degree) - min(node_degree))
  node_color <- ifelse(igraph::V(graph)$name %in% consensus_nodes, "#C44E52", "#4C72B0")
  set.seed(123)
  layout_graph <- igraph::layout_with_fr(graph)

  grDevices::cairo_pdf(paste0(file_base, ".pdf"), width = 7, height = 6, family = "Arial")
  ppi_base_plot_par()
  plot(graph, layout = layout_graph, vertex.label = igraph::V(graph)$name, vertex.label.cex = 0.9, vertex.label.color = "black", vertex.size = node_size, vertex.color = node_color, vertex.frame.color = "white", edge.color = "grey70", edge.width = 1.2, main = "PPI Network Highlighting Consensus Hubs")
  graphics::legend("topleft", legend = c("Consensus hub", "Other node"), col = c("#C44E52", "#4C72B0"), pch = 19, pt.cex = 1.5, bty = "n")
  grDevices::dev.off()

  grDevices::png(paste0(file_base, ".png"), width = 4200, height = 3600, res = 600)
  ppi_base_plot_par()
  plot(graph, layout = layout_graph, vertex.label = igraph::V(graph)$name, vertex.label.cex = 0.9, vertex.label.color = "black", vertex.size = node_size, vertex.color = node_color, vertex.frame.color = "white", edge.color = "grey70", edge.width = 1.2, main = "PPI Network Highlighting Consensus Hubs")
  graphics::legend("topleft", legend = c("Consensus hub", "Other node"), col = c("#C44E52", "#4C72B0"), pch = 19, pt.cex = 1.5, bty = "n")
  grDevices::dev.off()
}

ppi_base_plot_par <- function() {
  graphics::par(family = "Arial", cex.axis = 1.5, cex.lab = 1.67, cex.main = 2, font.lab = 2, font.main = 2)
}

ppi_formula_check_table <- function() {
  data.frame(
    Method = ppi_methods(),
    Formula_Source = c(
      "Deg(v)=|N(v)|",
      "MNC(v)=|V(MC(v))|",
      "DMNC(v)=|E(MC(v))|/|V(MC(v))|^epsilon, epsilon=1.7",
      "MCC(v)=sum_{C in S(v)} (|C|-1)!",
      "BN(v)=sum_s p_s(v), p_s(v)=1 if more than |V(Ts)|/4 shortest-path-tree paths meet at v",
      "EC(v)=|V(C(v))|/|V| * 1/max_dist(v,w in C(v))",
      "Clo(v)=sum_w 1/dist(v,w)",
      "Rad(v)=|V(C(v))|/|V| * sum_w (Delta_C(v)+1-dist(v,w))/Delta_C(v)",
      "BC(v)=sum_{s!=t!=v} sigma_st(v)/sigma_st",
      "Str(v)=sum_{s!=t!=v} sigma_st(v)",
      "EPC(v)=1/|V| * sum_t delta_vt over 1000 edge-percolated networks"
    ),
    Implementation_Status = rep("Implemented", length(ppi_methods())),
    stringsAsFactors = FALSE
  )
}

ppi_need_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for PPI analysis.", call. = FALSE)
  }
}
