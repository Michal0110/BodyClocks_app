server <- function(input, output, session) {

  # --------- Circadian loaders (dynamic) ----------
  load_network_circ <- function(ds, string_threshold = 700,
                                enrichment_background = "all_string") {
    meta <- CIRC_DATASETS[[ds]]
    if (is.null(meta)) return(NULL)

    if (!enrichment_background %in% ENRICHMENT_BACKGROUND_CODES) {
      stop(
        "Unknown enrichment background '", enrichment_background,
        "' for dataset '", ds, "'.", call. = FALSE
      )
    }

    folder <- meta$folder
    prefix <- meta$prefix

    edge_file <- file.path(folder, paste0(prefix, "_edges_string", string_threshold, ".rds"))
    if (!file.exists(edge_file) && string_threshold == 700)
      edge_file <- file.path(folder, paste0(prefix, "_edges.rds"))

    positions_file <- file.path(folder, paste0(prefix, "_positions_string", string_threshold, ".rds"))
    if (!file.exists(positions_file) && string_threshold == 700)
      positions_file <- file.path(folder, paste0(prefix, "_positions.rds"))

    legacy_nodes_file <- file.path(
      folder, paste0(prefix, "_nodes_with_tooltips.rds")
    )
    enrichment_bundle_file <- file.path(
      folder, paste0(prefix, "_enrichment_backgrounds.rds")
    )
    nodes_with_tooltips <- readRDS(legacy_nodes_file)
    using_legacy_enrichment <- FALSE

    if (file.exists(enrichment_bundle_file)) {
      bundle <- tryCatch(
        readRDS(enrichment_bundle_file),
        error = function(e) {
          stop(
            "Could not read enrichment backgrounds for dataset '", ds,
            "': ", conditionMessage(e), call. = FALSE
          )
        }
      )
      enrichment <- select_enrichment_background(
        bundle, ds, enrichment_background
      )
      nodes_with_tooltips <- apply_enrichment_node_titles(
        nodes_with_tooltips, enrichment$node_titles, ds,
        enrichment_background
      )
    } else if (identical(enrichment_background, "all_string")) {
      # Backward compatibility for deployments created before background-
      # specific enrichment bundles were introduced.
      enrichment <- normalise_legacy_enrichment(
        readRDS(file.path(
          folder, paste0(prefix, "_cluster_descriptions.rds")
        )),
        readRDS(file.path(
          folder, paste0(prefix, "_cluster_for_plot.rds")
        ))
      )
      using_legacy_enrichment <- TRUE
    } else {
      dataset_label <- if (!is.null(meta$dataset_label) &&
                           nzchar(meta$dataset_label)) {
        meta$dataset_label
      } else {
        ds
      }
      stop(
        "Dataset-tested enrichment is not available for '", dataset_label,
        "'. Expected ", enrichment_bundle_file,
        ". Rebuild and deploy this dataset's background-specific enrichment ",
        "bundle, or select 'All STRING proteins'.",
        call. = FALSE
      )
    }

    base_data <- list(
      nodes_with_tooltips  = nodes_with_tooltips,
      edges                = readRDS(edge_file),
      cluster_to_genes     = enrichment$cluster_to_genes,
      cluster_descriptions = enrichment$cluster_descriptions,
      enrichment_background = enrichment_background
    )

    # Load cached positions if available
    if (file.exists(positions_file)) {
      base_data$positions <- readRDS(positions_file)
    }

    if (using_legacy_enrichment &&
        !is.null(base_data$nodes_with_tooltips) &&
        "title" %in% names(base_data$nodes_with_tooltips)) {
      base_data$nodes_with_tooltips$title <- vapply(
        base_data$nodes_with_tooltips$title,
        limit_go_terms_in_tooltip,
        character(1)
      )
    }
    base_data
  }

  load_plot_circ <- function(ds) {
    meta <- CIRC_DATASETS[[ds]]
    if (is.null(meta)) return(NULL)

    if (meta$col_type == "single") {
      data_file <- file.path(meta$folder, paste0(meta$prefix, "_data.rds"))
      list(
        data            = attach_curve_cache(readRDS(data_file), data_file, meta$period, meta$t),
        period          = meta$period,
        t               = meta$t,
        table           = readRDS(file.path(meta$folder, paste0(meta$prefix, "_display_dt.rds"))),
        columnsToFormat = meta$columnsToFormat
      )
    } else {
      # wt_primary or ko_primary — shared folder with explicit file paths
      list(
        data_wt         = attach_curve_cache(readRDS(meta$wt_data_file), meta$wt_data_file, meta$period, meta$t),
        data_ko         = attach_curve_cache(readRDS(meta$ko_data_file), meta$ko_data_file, meta$period, meta$t),
        period          = meta$period,
        t               = meta$t,
        table           = readRDS(meta$display_file),
        columnsToFormat = meta$columnsToFormat
      )
    }
  }

  generate_plot_circ <- function(data_info, genes, input, session) {
    plots <- list()
    for (gene in genes) {
      ds   <- input$dataset
      meta <- CIRC_DATASETS[[ds]]
      if (is.null(meta)) next
      p <- switch(meta$plot_fn,
        osmo    = fitted_plot_osmo(data_info$data, gene, data_info$period, data_info$t,
                                   session, input$show_fit, input$show_line),
        atlas   = fitted_plot_atlas(data_info$data, NULL, gene, data_info$period,
                                    session, input$show_fit, input$show_line),
        cart_wt = fitted_plot_cart(
                                   if (!is.null(data_info$data_wt)) data_info$data_wt else data_info$data,
                                   data_info$data_ko, gene,
                                   data_info$period, data_info$t, session,
                                   input$show_fit, input$show_line, input$show_comparison,
                                   is_primary_wt = TRUE),
        cart_ko = fitted_plot_cart(data_info$data_ko, data_info$data_wt, gene,
                                   data_info$period, data_info$t, session,
                                   input$show_fit, input$show_line, input$show_comparison,
                                   is_primary_wt = FALSE),
        gobs         = fitted_plot_gobs(data_info$data, gene, data_info$period, data_info$t,
                                        session, input$show_fit, input$show_line),
        zt24         = fitted_plot_zt24(data_info$data, gene, data_info$period, data_info$t,
                                        session, input$show_fit, input$show_line),
        baboon_atlas = fitted_plot_baboon_atlas(data_info$data, NULL, gene, data_info$period,
                                               session, input$show_fit, input$show_line),
        NULL
      )
      if (!is.null(p)) plots[[gene]] <- p
    }
    plots
  }

  # --------- Start explorer (namespaced) ----------
  explorerServer(
    id                         = "circ",
    species_choices            = CIRC_SPECIES_CHOICES,
    tissue_choices_by_species  = CIRC_TISSUE_CHOICES_BY_SPECIES,
    dataset_choices_by_species = CIRC_DATASET_CHOICES_BY_SPECIES,
    load_network_fn            = load_network_circ,
    load_plotdata_fn           = load_plot_circ,
    generate_plot_fn           = generate_plot_circ,
    legend_html                = "www/legend.html"
  )

  comparisonServer(
    id                         = "comp",
    species_choices            = CIRC_SPECIES_CHOICES,
    tissue_choices_by_species  = CIRC_TISSUE_CHOICES_BY_SPECIES,
    dataset_choices_by_species = CIRC_DATASET_CHOICES_BY_SPECIES,
    load_network_fn            = load_network_circ,
    load_plotdata_fn           = load_plot_circ
  )
}
