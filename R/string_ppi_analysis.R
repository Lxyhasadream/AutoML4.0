# Offline STRING PPI construction followed by cytoHubba-style hub screening.

auto_string_ppi_analysis <- function(genes,
                                     resource_dir = "resource",
                                     output_dir = "STRING_PPI_analysis_results",
                                     species_id = "9606",
                                     string_version = "v12.0",
                                     aliases_file = NULL,
                                     info_file = NULL,
                                     links_file = NULL,
                                     download_resources = TRUE,
                                     overwrite_resources = FALSE,
                                     score_threshold = 400,
                                     top_n = 10,
                                     consensus_min_methods = NULL,
                                     alias_chunk_size = 300000,
                                     links_chunk_size = 500000,
                                     epsilon = 1.7,
                                     epc_threshold = 0.1,
                                     epc_n_sim = 1000,
                                     seed = 123,
                                     include_self_radiality = TRUE,
                                     write_outputs = TRUE,
                                     plot_outputs = TRUE,
                                     verbose = TRUE) {
  string_ppi_need_package("data.table")
  ppi_need_package("igraph")

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    stop("genes must contain at least one non-empty gene symbol.", call. = FALSE)
  }

  score_threshold <- as.numeric(score_threshold)
  if (is.na(score_threshold) || score_threshold < 0) {
    stop("score_threshold must be a non-negative number.", call. = FALSE)
  }

  if (write_outputs || plot_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  if (verbose) {
    message("Input genes: ", length(genes))
    message("STRING score threshold: ", score_threshold)
    message("Resource directory: ", normalizePath(resource_dir, mustWork = FALSE))
  }

  string_result <- string_build_gene_ppi_edges(
    genes = genes,
    resource_dir = resource_dir,
    output_dir = output_dir,
    species_id = species_id,
    string_version = string_version,
    aliases_file = aliases_file,
    info_file = info_file,
    links_file = links_file,
    download_resources = download_resources,
    overwrite_resources = overwrite_resources,
    score_threshold = score_threshold,
    alias_chunk_size = alias_chunk_size,
    links_chunk_size = links_chunk_size,
    write_outputs = write_outputs,
    verbose = verbose
  )

  hub_result <- auto_ppi_analysis(
    edge_df = string_result$edge_df,
    from_col = "source",
    to_col = "target",
    output_dir = output_dir,
    top_n = top_n,
    consensus_min_methods = consensus_min_methods,
    epsilon = epsilon,
    epc_threshold = epc_threshold,
    epc_n_sim = epc_n_sim,
    seed = seed,
    include_self_radiality = include_self_radiality,
    write_outputs = write_outputs,
    plot_outputs = plot_outputs,
    verbose = verbose
  )

  strict_intersection_genes <- sort(unique(Reduce(intersect, hub_result$top_nodes)))
  consensus_hub_genes <- sort(unique(hub_result$consensus_hubs$Node))

  if (write_outputs) {
    utils::write.csv(
      data.frame(gene = strict_intersection_genes, stringsAsFactors = FALSE),
      file.path(output_dir, paste0("STRING_PPI_top", top_n, "_intersection_all_11_methods.csv")),
      row.names = FALSE
    )
    utils::write.csv(
      data.frame(gene = consensus_hub_genes, stringsAsFactors = FALSE),
      file.path(output_dir, "STRING_PPI_consensus_hub_genes.csv"),
      row.names = FALSE
    )
  }

  structure(
    list(
      genes = genes,
      resource_dir = normalizePath(resource_dir, mustWork = FALSE),
      resource_files = string_result$resource_files,
      score_threshold = score_threshold,
      mapping = string_result$mapping,
      unmapped_genes = string_result$unmapped_genes,
      raw_string_interactions = string_result$raw_string_interactions,
      edge_df = string_result$edge_df,
      ppi = hub_result,
      intersection_genes_all_11_methods = strict_intersection_genes,
      consensus_hub_genes = consensus_hub_genes,
      downstream_genes = strict_intersection_genes,
      output_dir = if (write_outputs || plot_outputs) normalizePath(output_dir, mustWork = FALSE) else NA_character_
    ),
    class = "auto_string_ppi_analysis"
  )
}

string_build_gene_ppi_edges <- function(genes,
                                        resource_dir = "resource",
                                        output_dir = "STRING_PPI_analysis_results",
                                        species_id = "9606",
                                        string_version = "v12.0",
                                        aliases_file = NULL,
                                        info_file = NULL,
                                        links_file = NULL,
                                        download_resources = TRUE,
                                        overwrite_resources = FALSE,
                                        score_threshold = 400,
                                        alias_chunk_size = 300000,
                                        links_chunk_size = 500000,
                                        write_outputs = TRUE,
                                        verbose = TRUE) {
  string_ppi_need_package("data.table")

  genes <- unique(as.character(genes))
  genes <- genes[!is.na(genes) & nzchar(genes)]
  if (length(genes) == 0) {
    stop("genes must contain at least one non-empty gene symbol.", call. = FALSE)
  }

  files <- string_resource_files(
    resource_dir = resource_dir,
    species_id = species_id,
    string_version = string_version,
    aliases_file = aliases_file,
    info_file = info_file,
    links_file = links_file,
    download_resources = download_resources,
    overwrite = overwrite_resources,
    verbose = verbose
  )
  if (write_outputs) {
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  }

  if (verbose) {
    message("Reading STRING protein.info file...")
  }
  info_dt <- string_safe_fread(files$info, sep = "\t")
  info_dt <- string_standardize_id_col(info_dt)
  if (!"preferred_name" %in% colnames(info_dt)) {
    stop("Column 'preferred_name' was not found in protein.info file.", call. = FALSE)
  }
  info_dt <- info_dt[, c("STRING_id", "preferred_name"), with = FALSE]
  data.table::set(info_dt, j = "STRING_id", value = as.character(info_dt$STRING_id))
  data.table::set(info_dt, j = "preferred_name", value = as.character(info_dt$preferred_name))
  info_dt <- unique(info_dt)
  data.table::set(info_dt, j = "preferred_upper", value = toupper(info_dt$preferred_name))

  gene_map <- data.table::data.table(gene_input = genes, gene_upper = toupper(genes))
  map_preferred <- merge(
    gene_map,
    info_dt,
    by.x = "gene_upper",
    by.y = "preferred_upper",
    all = FALSE,
    allow.cartesian = TRUE
  )

  if (nrow(map_preferred) > 0) {
    data.table::set(map_preferred, j = "mapping_source", value = "preferred_name")
    map_preferred <- map_preferred[, c("gene_input", "STRING_id", "preferred_name", "mapping_source"), with = FALSE]
    data.table::setnames(map_preferred, "preferred_name", "matched_name")
  } else {
    map_preferred <- data.table::data.table(
      gene_input = character(),
      STRING_id = character(),
      matched_name = character(),
      mapping_source = character()
    )
  }

  genes_for_alias <- setdiff(genes, unique(map_preferred$gene_input))
  if (verbose) {
    message("Genes mapped by preferred_name: ", data.table::uniqueN(map_preferred$gene_input))
    message("Genes requiring alias search: ", length(genes_for_alias))
  }

  map_alias <- string_stream_alias_mapping(
    aliases_file = files$aliases,
    genes_for_alias = genes_for_alias,
    gene_map = gene_map,
    chunk_size = alias_chunk_size,
    verbose = verbose
  )

  mapped <- data.table::rbindlist(list(map_preferred, map_alias), use.names = TRUE, fill = TRUE)
  mapped <- unique(mapped)
  mapped <- merge(mapped, info_dt[, c("STRING_id", "preferred_name"), with = FALSE], by = "STRING_id", all.x = TRUE)
  data.table::set(mapped, j = "mapping_priority", value = ifelse(mapped$mapping_source == "preferred_name", 1, 2))
  data.table::setorder(mapped, gene_input, mapping_priority, STRING_id)

  mapped_genes <- unique(mapped$gene_input)
  unmapped_genes <- setdiff(genes, mapped_genes)

  if (write_outputs) {
    utils::write.csv(mapped, file.path(output_dir, "STRING_offline_mapped_genes.csv"), row.names = FALSE)
    utils::write.csv(data.frame(unmapped_gene = unmapped_genes), file.path(output_dir, "STRING_offline_unmapped_genes.csv"), row.names = FALSE)
  }

  if (verbose) {
    message("Genes mapped by alias: ", data.table::uniqueN(map_alias$gene_input))
    message("Mapped genes: ", length(mapped_genes), "/", length(genes))
  }
  if (nrow(mapped) == 0) {
    stop("No genes were mapped to STRING IDs. Please check gene symbols and species.", call. = FALSE)
  }

  query_string_ids <- unique(mapped$STRING_id)
  if (verbose) {
    message("Extracting interactions from local STRING links file...")
    message("Mapped STRING IDs: ", length(query_string_ids))
  }

  ppi_raw <- string_stream_links_among_query_ids(
    links_file = files$links,
    query_string_ids = query_string_ids,
    score_threshold = score_threshold,
    chunk_size = links_chunk_size,
    verbose = verbose
  )

  if (write_outputs) {
    utils::write.csv(ppi_raw, file.path(output_dir, "STRING_offline_raw_query_interactions.csv"), row.names = FALSE)
  }
  if (nrow(ppi_raw) == 0) {
    stop(
      "No STRING interactions were found among the mapped query proteins at the current threshold. ",
      "Try score_threshold = 150 or expand the input gene list.",
      call. = FALSE
    )
  }

  edge_df <- string_collapse_edges_to_genes(ppi_raw, mapped)
  if (write_outputs) {
    utils::write.csv(edge_df, file.path(output_dir, "STRING_offline_edge_df_for_cytohubba.csv"), row.names = FALSE)
  }
  if (nrow(edge_df) == 0) {
    stop(
      "edge_df is empty. No downstream PPI topology analysis can be performed. ",
      "Consider score_threshold = 150 or expanding the input gene list.",
      call. = FALSE
    )
  }

  if (verbose) {
    message("Final gene-level edges: ", nrow(edge_df))
  }

  list(
    resource_files = files,
    mapping = mapped,
    unmapped_genes = unmapped_genes,
    raw_string_interactions = ppi_raw,
    edge_df = edge_df
  )
}

download_string_resources <- function(resource_dir = "resource",
                                      species_id = "9606",
                                      string_version = "v12.0",
                                      overwrite = FALSE,
                                      timeout = 600,
                                      verbose = TRUE) {
  dir.create(resource_dir, showWarnings = FALSE, recursive = TRUE)
  files <- string_resource_paths(resource_dir, species_id, string_version)
  urls <- string_resource_urls(species_id, string_version)

  old_timeout <- getOption("timeout")
  options(timeout = max(timeout, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)

  for (nm in names(files)) {
    if (file.exists(files[[nm]]) && !overwrite) {
      if (verbose) {
        message("STRING resource exists: ", files[[nm]])
      }
      next
    }
    if (verbose) {
      message("Downloading STRING ", nm, " resource from: ", urls[[nm]])
    }
    status <- tryCatch(
      utils::download.file(
        url = urls[[nm]],
        destfile = files[[nm]],
        mode = "wb",
        quiet = !verbose
      ),
      error = function(e) {
        if (file.exists(files[[nm]])) {
          unlink(files[[nm]])
        }
        stop("Failed to download STRING resource: ", urls[[nm]], "\n", conditionMessage(e), call. = FALSE)
      }
    )
    if (!identical(status, 0L) || !file.exists(files[[nm]])) {
      if (file.exists(files[[nm]])) {
        unlink(files[[nm]])
      }
      stop("Failed to download STRING resource: ", urls[[nm]], call. = FALSE)
    }
  }

  invisible(files)
}

string_resource_files <- function(resource_dir,
                                  species_id = "9606",
                                  string_version = "v12.0",
                                  aliases_file = NULL,
                                  info_file = NULL,
                                  links_file = NULL,
                                  download_resources = TRUE,
                                  overwrite = FALSE,
                                  verbose = TRUE) {
  manual_files <- string_manual_resource_files(
    aliases_file = aliases_file,
    info_file = info_file,
    links_file = links_file
  )
  if (!is.null(manual_files)) {
    if (verbose) {
      message("Using manually specified local STRING resource files.")
    }
    return(manual_files)
  }

  files <- string_resource_paths(resource_dir, species_id, string_version)
  missing <- names(files)[!file.exists(unlist(files, use.names = FALSE))]

  if (length(missing) > 0 && download_resources) {
    if (verbose) {
      message("Missing STRING resource file(s); downloading required files.")
    }
    download_string_resources(
      resource_dir = resource_dir,
      species_id = species_id,
      string_version = string_version,
      overwrite = overwrite,
      verbose = verbose
    )
  }

  missing <- names(files)[!file.exists(unlist(files, use.names = FALSE))]
  if (length(missing) > 0) {
    stop(
      "Missing STRING resource file(s): ",
      paste(sprintf("%s=%s", missing, unlist(files[missing], use.names = FALSE)), collapse = "; "),
      "\nRun download_string_resources(resource_dir = '", resource_dir, "', species_id = '", species_id, "') or set download_resources = TRUE.",
      call. = FALSE
    )
  }
  files
}

string_manual_resource_files <- function(aliases_file = NULL, info_file = NULL, links_file = NULL) {
  supplied <- !vapply(list(aliases_file, info_file, links_file), is.null, logical(1))
  if (!any(supplied)) {
    return(NULL)
  }
  if (!all(supplied)) {
    stop(
      "When manually specifying STRING resource files, aliases_file, info_file, and links_file must all be provided.",
      call. = FALSE
    )
  }

  files <- list(
    aliases = aliases_file,
    info = info_file,
    links = links_file
  )
  missing <- names(files)[!file.exists(unlist(files, use.names = FALSE))]
  if (length(missing) > 0) {
    stop(
      "Manually specified STRING file(s) do not exist: ",
      paste(sprintf("%s=%s", missing, unlist(files[missing], use.names = FALSE)), collapse = "; "),
      call. = FALSE
    )
  }
  files
}

string_resource_paths <- function(resource_dir, species_id = "9606", string_version = "v12.0") {
  files <- list(
    aliases = file.path(resource_dir, paste0(species_id, ".protein.aliases.", string_version, ".txt.gz")),
    info = file.path(resource_dir, paste0(species_id, ".protein.info.", string_version, ".txt.gz")),
    links = file.path(resource_dir, paste0(species_id, ".protein.links.", string_version, ".txt.gz"))
  )
  files
}

string_resource_urls <- function(species_id = "9606", string_version = "v12.0") {
  file_stems <- c(
    aliases = "protein.aliases",
    info = "protein.info",
    links = "protein.links"
  )
  stats::setNames(
    vapply(file_stems, function(stem) {
      file_name <- paste0(species_id, ".", stem, ".", string_version, ".txt.gz")
      paste0("https://stringdb-downloads.org/download/", stem, ".", string_version, "/", file_name)
    }, character(1)),
    names(file_stems)
  )
}

string_standardize_id_col <- function(dt) {
  id_col <- grep("string_protein_id", colnames(dt), value = TRUE)
  id_col <- if (length(id_col) == 0) colnames(dt)[1] else id_col[1]
  data.table::setnames(dt, id_col, "STRING_id")
  dt
}

string_safe_fread <- function(file, sep = "\t") {
  data.table::fread(
    file,
    sep = sep,
    header = TRUE,
    quote = "",
    data.table = TRUE,
    showProgress = FALSE
  )
}

string_stream_alias_mapping <- function(aliases_file,
                                        genes_for_alias,
                                        gene_map,
                                        chunk_size = 300000,
                                        verbose = TRUE) {
  genes_for_alias <- unique(as.character(genes_for_alias))
  if (length(genes_for_alias) == 0) {
    return(data.table::data.table(
      gene_input = character(),
      STRING_id = character(),
      matched_name = character(),
      mapping_source = character()
    ))
  }

  target_upper <- unique(toupper(genes_for_alias))
  con <- gzfile(aliases_file, open = "rt")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1)
  header_cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  id_col <- grep("string_protein_id", header_cols, value = TRUE)
  id_col <- if (length(id_col) == 0) header_cols[1] else id_col[1]
  if (!"alias" %in% header_cols) {
    stop("Column 'alias' was not found in aliases file header.", call. = FALSE)
  }

  res_list <- list()
  chunk_id <- 1L
  total_hits <- 0L

  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) {
      break
    }

    dt <- data.table::fread(
      text = paste(lines, collapse = "\n"),
      sep = "\t",
      header = FALSE,
      col.names = header_cols,
      quote = "",
      fill = TRUE,
      data.table = TRUE,
      showProgress = FALSE
    )
    data.table::setnames(dt, id_col, "STRING_id")
    dt <- dt[, c("STRING_id", "alias"), with = FALSE]
    data.table::set(dt, j = "STRING_id", value = as.character(dt$STRING_id))
    data.table::set(dt, j = "alias", value = as.character(dt$alias))
    data.table::set(dt, j = "alias_upper", value = toupper(dt$alias))
    dt <- dt[alias_upper %in% target_upper]

    if (nrow(dt) > 0) {
      tmp <- merge(
        gene_map[gene_input %in% genes_for_alias],
        dt,
        by.x = "gene_upper",
        by.y = "alias_upper",
        all = FALSE,
        allow.cartesian = TRUE
      )
      if (nrow(tmp) > 0) {
        tmp <- tmp[, c("gene_input", "STRING_id", "alias"), with = FALSE]
        data.table::setnames(tmp, "alias", "matched_name")
        data.table::set(tmp, j = "mapping_source", value = "alias")
        res_list[[length(res_list) + 1L]] <- tmp
        total_hits <- total_hits + nrow(tmp)
      }
    }

    if (verbose && chunk_id %% 10L == 0L) {
      message("Alias chunks processed: ", chunk_id, " | hits: ", total_hits)
    }
    chunk_id <- chunk_id + 1L
  }

  if (length(res_list) == 0) {
    return(data.table::data.table(
      gene_input = character(),
      STRING_id = character(),
      matched_name = character(),
      mapping_source = character()
    ))
  }
  unique(data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE))
}

string_stream_links_among_query_ids <- function(links_file,
                                                query_string_ids,
                                                score_threshold = 400,
                                                chunk_size = 500000,
                                                verbose = TRUE) {
  query_string_ids <- unique(as.character(query_string_ids))
  con <- gzfile(links_file, open = "rt")
  on.exit(close(con), add = TRUE)

  header <- readLines(con, n = 1)
  header_cols <- strsplit(header, " ", fixed = TRUE)[[1]]
  header_cols <- header_cols[nzchar(header_cols)]
  required_cols <- c("protein1", "protein2", "combined_score")
  if (!all(required_cols %in% header_cols)) {
    stop(
      "The links file header does not contain expected columns: ",
      paste(required_cols, collapse = ", "),
      call. = FALSE
    )
  }

  res_list <- list()
  chunk_id <- 1L
  total_edges <- 0L

  repeat {
    lines <- readLines(con, n = chunk_size)
    if (length(lines) == 0) {
      break
    }

    dt <- data.table::fread(
      text = paste(lines, collapse = "\n"),
      sep = " ",
      header = FALSE,
      col.names = header_cols,
      quote = "",
      fill = TRUE,
      data.table = TRUE,
      showProgress = FALSE
    )
    dt <- dt[, c("protein1", "protein2", "combined_score"), with = FALSE]
    data.table::set(dt, j = "protein1", value = as.character(dt$protein1))
    data.table::set(dt, j = "protein2", value = as.character(dt$protein2))
    data.table::set(dt, j = "combined_score", value = as.numeric(dt$combined_score))
    dt <- dt[
      protein1 %in% query_string_ids &
        protein2 %in% query_string_ids &
        combined_score >= score_threshold
    ]

    if (nrow(dt) > 0) {
      res_list[[length(res_list) + 1L]] <- dt
      total_edges <- total_edges + nrow(dt)
    }
    if (verbose && chunk_id %% 10L == 0L) {
      message("Link chunks processed: ", chunk_id, " | matched edges: ", total_edges)
    }
    chunk_id <- chunk_id + 1L
  }

  if (length(res_list) == 0) {
    return(data.table::data.table(
      protein1 = character(),
      protein2 = character(),
      combined_score = numeric()
    ))
  }
  unique(data.table::rbindlist(res_list, use.names = TRUE, fill = TRUE))
}

string_collapse_edges_to_genes <- function(ppi_raw, mapped) {
  id_to_gene <- mapped[
    ,
    list(
      gene_symbol = sort(unique(gene_input))[1],
      all_mapped_input_genes = paste(sort(unique(gene_input)), collapse = ";"),
      preferred_name = sort(unique(preferred_name))[1]
    ),
    by = "STRING_id"
  ]

  edge_dt <- merge(
    ppi_raw,
    id_to_gene[, c("STRING_id", "gene_symbol", "preferred_name"), with = FALSE],
    by.x = "protein1",
    by.y = "STRING_id",
    all.x = TRUE
  )
  data.table::setnames(edge_dt, c("gene_symbol", "preferred_name"), c("source", "source_preferred_name"))

  edge_dt <- merge(
    edge_dt,
    id_to_gene[, c("STRING_id", "gene_symbol", "preferred_name"), with = FALSE],
    by.x = "protein2",
    by.y = "STRING_id",
    all.x = TRUE
  )
  data.table::setnames(edge_dt, c("gene_symbol", "preferred_name"), c("target", "target_preferred_name"))
  edge_dt <- edge_dt[!is.na(source) & !is.na(target) & source != target]
  data.table::set(edge_dt, j = "geneA", value = pmin(edge_dt$source, edge_dt$target))
  data.table::set(edge_dt, j = "geneB", value = pmax(edge_dt$source, edge_dt$target))

  edge_df <- edge_dt[
    ,
    list(
      combined_score = max(combined_score, na.rm = TRUE),
      protein_pair_count = .N,
      STRING_protein_pairs = paste(unique(paste(protein1, protein2, sep = "--")), collapse = ";")
    ),
    by = list(source = geneA, target = geneB)
  ]
  data.table::setorder(edge_df, -combined_score, source, target)
  as.data.frame(edge_df)
}

string_ppi_need_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required for offline STRING PPI analysis.", call. = FALSE)
  }
}
