explorerUI <- function(id, title,
                       species_choices,
                       tissue_choices_by_species,
                       dataset_choices_by_species,
                       legend_html = "www/legend.html",
                       show_plot_controls = TRUE) {
  ns <- NS(id)
  default_species <- species_choices[[1]]
  default_tissues <- tissue_choices_by_species[[default_species]]
  tagList(
    fluidRow(
      column(
        width = 3,
        div(class = "explorer-sidebar bc-sidebar bc-panel",
          # Species section
          div(class = "section-header", "Species:"),
          selectInput(
            ns("species"),
            NULL,
            choices  = species_choices,
            selected = default_species
          ),

          # Tissue section
          div(class = "section-header", "Tissue:",
              tags$span(icon("question-circle"), id = ns("tissue_help"))
          ),
          selectInput(
            ns("tissue"),
            NULL,
            choices = c("Select a tissue..." = "", default_tissues),
            selected = ""
          ),

          # Dataset section (choices updated by server based on tissue)
          div(class = "section-header", "Dataset:",
              tags$span(icon("question-circle"), id = ns("dataset_help"))
          ),
          selectInput(
            ns("dataset"),
            NULL,
            choices = c("Select a dataset..." = ""),
            selected = ""
          ),

          # STRING network stringency
          div(class = "section-header", "STRING Network:"),
          radioButtons(
            ns("string_threshold"),
            NULL,
            choices  = c("Inclusive (700)" = "700", "Stringent (900)" = "900"),
            selected = "700",
            inline   = TRUE
          ),

          # Gene
          div(class = "section-header", "Gene:",
              tags$span(icon("question-circle"), id = ns("gene_help"))
          ),
          textInput(ns("gene"), NULL, placeholder = "Enter gene name or list separated by commas"),
          actionButton(ns("plot"), "Plot", class = "btn-primary action-button"),
          br(), br(),

          # Enrichment controls (NO landing-section here)
          div(class = "section-header", "Functional Enrichment",
              tags$span(icon("question-circle"), id = ns("enrichment_help"))
          ),
          div(class = "enrichment-controls",
            radioButtons(
              ns("enrichment_background"),
              "Background:",
              choices = c(
                "Genes tested in this dataset" = "dataset_tested",
                "Entire genome" = "all_string"
              ),
              selected = "all_string",
              inline = FALSE
            ),
            selectInput(ns("enrichment_category"), "Category:",
                        choices = c("Select category" = "",
                                    "GO Biological Process" = "Process",
                                    "KEGG Pathway" = "KEGG",
                                    "Reactome" = "RCTM",
                                    "WikiPathways" = "WikiPathways")),
            div(class = "scrollable-checkbox",
              checkboxGroupInput(ns("selected_term"), NULL, choices = NULL, selected = NULL)
            ),
            actionButton(ns("plot_enrichment"), "Plot selected genes", class = "btn-primary action-button")
          ),
          br(), br(),

          # Data table (light blue container, NO landing-section)
          div(class = "data-table-container",
            div(class = "section-header data-table-title", "Data Table",
                tags$span(icon("question-circle"), id = ns("data_table_help"))
            ),
            DTOutput(ns("dataTable")),
            fluidRow(
              column(6, numericInput(ns("pval_cutoff"), "pVal \u2264",
                                    value = 1, min = 0, max = 1, step = 0.001)),
              column(6, numericInput(ns("bhq_cutoff"),  "BH.Q \u2264",
                                    value = 1, min = 0, max = 1, step = 0.001))
            )
          ),

          tags$script(HTML(sprintf(
            "$(document).on('keydown', '#%s', function(e) {
                if (e.key === 'Enter') { $('#%s').click(); }
            });",
            ns('gene'), ns('plot')
          ))),

          actionButton(ns("clearSelection"), "Clear Selection", class = "btn-secondary action-button")
        )
      ),
      column(
        width = 9,
        div(class = "explorer-main",
          # "Coming soon" placeholder shown when selected species has no data
          uiOutput(ns("coming_soon")),
          div(class = "network-container bc-network-panel bc-panel",
              visNetworkOutput(ns("network"), height = "600px")),
          # Dynamically switched legend
          div(class = "bc-panel bc-legend-panel", uiOutput(ns("legend"))),
           # Show/hide plot controls
          if (isTRUE(show_plot_controls)) {
           div(class = "plot-controls-container bc-toolbar bc-panel",
             div(class = "plot-controls",
               span(class = "toolbar-label", "Graph Controls"),
               checkboxInput(ns("show_fit"), "Show Sine Fitting", value = FALSE),
               checkboxInput(ns("show_line"), "Show Connecting Line", value = TRUE),
               checkboxInput(ns("show_comparison"), "Show Comparison Dataset", value = FALSE),
               div(class = "tooltip-icon", icon("question-circle"), id = ns("graph_help"))
             )
           )
         },
          div(class = "plots-container bc-plots-panel bc-panel", uiOutput(ns("plotsUI")))
        )
      )
    )
  )
}

explorerServer <- function(id,
                           species_choices,
                           tissue_choices_by_species,
                           dataset_choices_by_species,
                           load_network_fn,
                           load_plotdata_fn,
                           generate_plot_fn,
                           legend_html = "www/legend.html"
) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values
    selected_cluster_genes <- reactiveVal(NULL)
    gene_list <- reactiveVal(NULL)
    plot_batch_size <- 4L
    plot_page_size <- 32L
    plot_slot_count <- reactiveVal(0L)
    visible_plot_count <- reactiveVal(0L)

    # Whether the selected species has any datasets
    species_has_data <- reactive({
      req(input$species)
      length(tissue_choices_by_species[[input$species]]) > 0
    })

    # "Coming soon" panel for species without data yet
    output$coming_soon <- renderUI({
      if (species_has_data()) return(NULL)
      sp_label <- species_choices[[input$species]]
      div(class = "bc-panel bc-empty-state",
          h3(paste(sp_label, "datasets")),
          p("Data for this species is not yet available.",
            class = "bc-empty-copy"),
          p(tags$em("Coming soon."), class = "bc-empty-copy")
      )
    })

    # When species changes, reset tissue / dataset dropdowns
    observeEvent(input$species, {
      tc <- tissue_choices_by_species[[input$species]]
      if (length(tc) == 0) tc <- character(0)
      updateSelectInput(session, "tissue",
                        choices  = c("Select a tissue..." = "", tc),
                        selected = "")
      updateSelectInput(session, "dataset",
                        choices  = c("Select a dataset..." = ""),
                        selected = "")
      selected_cluster_genes(NULL)
      gene_list(NULL)
    }, ignoreInit = TRUE)

    # When tissue changes, update dataset choices
    observeEvent(input$tissue, {
      if (!nzchar(input$tissue)) {
        updateSelectInput(session, "dataset", choices = c("Select a dataset..." = ""), selected = "")
        return()
      }
      ds_choices <- dataset_choices_by_species[[input$species]][[input$tissue]]
      if (is.null(ds_choices) || length(ds_choices) == 0) {
        ds_choices <- c("No datasets available" = "")
      }
      # Auto-select when there is only one option (e.g. baboon tissues)
      auto_selected <- if (length(ds_choices) == 1) unname(ds_choices)[[1]] else ""
      updateSelectInput(session, "dataset",
                        choices  = c("Select a dataset..." = "", ds_choices),
                        selected = auto_selected)
    }, ignoreInit = TRUE)

    # Reset enrichment/plots when tissue changes
    observeEvent(input$tissue, {
      selected_cluster_genes(NULL)
      gene_list(NULL)
      updateSelectInput(session, "enrichment_category", selected = "")
      updateCheckboxGroupInput(session, "selected_term", choices = NULL, selected = NULL)
      session$sendCustomMessage(type = "clearCheckboxGroup", message = ns("selected_term"))
    }, priority = 100, ignoreInit = TRUE)

    # The STRING threshold controls edges and positions only. Enrichment terms
    # remain selected when switching between the 700 and 900 networks.

    # Reset enrichment controls/highlighting when its statistical background
    # changes, without clearing expression plots or the background selection
    # itself. The radio selection naturally persists across dataset changes.
    observeEvent(input$enrichment_background, {
      selected_cluster_genes(NULL)
      updateSelectInput(session, "enrichment_category", selected = "")
      updateCheckboxGroupInput(session, "selected_term", choices = NULL, selected = NULL)
      session$sendCustomMessage(type = "clearCheckboxGroup", message = ns("selected_term"))
    }, ignoreInit = TRUE)

    # Reset enrichment UI when dataset changes (do nothing if empty)
    observeEvent(input$dataset, {
      if (!nzchar(input$dataset)) {
        # Clear local state when user clears selection
        selected_cluster_genes(NULL)
        gene_list(NULL)
        updateSelectInput(session, "enrichment_category", selected = "")
        updateCheckboxGroupInput(session, "selected_term", choices = NULL, selected = NULL)
        session$sendCustomMessage(type = "clearCheckboxGroup", message = ns("selected_term"))
        return()
      }

      # Dataset changed to a real id: reset UI/state and restore node styling
      selected_cluster_genes(NULL)
      # Clear plots before loading the new network. If a legacy deployment is
      # missing the requested background bundle, selected_network_data() below
      # raises a user-facing validation error; old plots must not survive it.
      gene_list(NULL)
      updateSelectInput(session, "enrichment_category", selected = "")
      updateCheckboxGroupInput(session, "selected_term", choices = NULL, selected = NULL)
      session$sendCustomMessage(type = "clearCheckboxGroup", message = ns("selected_term"))

      isolate({
        data <- selected_network_data()
        if (!is.null(data) && !is.null(data$nodes_with_tooltips)) {
          visNetworkProxy(ns("network")) %>%
            visUpdateNodes(nodes = data.frame(
              id = data$nodes_with_tooltips$id,
              borderWidth = 1,
              color = data$nodes_with_tooltips$color
            ))
        }
      })
    }, priority = 100, ignoreInit = TRUE)

    output$legend <- renderUI({
      ds   <- input$dataset
      meta <- CIRC_DATASETS[[ds]]
      html_file <- if (!is.null(meta) && meta$plot_fn %in% c("atlas", "baboon_atlas"))
        "www/legend_atlas.html"
      else
        legend_html
      includeHTML(html_file)
    })

    # Load network/plot data only when a dataset is selected and species has data
    selected_network_data <- reactive({
      req(species_has_data(), nzchar(input$dataset))
      req(input$string_threshold, input$enrichment_background)
      network_data <- tryCatch(
        load_network_fn(
          input$dataset,
          string_threshold = as.integer(input$string_threshold),
          enrichment_background = input$enrichment_background
        ),
        error = function(e) {
          shiny::validate(shiny::need(FALSE, conditionMessage(e)))
        }
      )
      shiny::validate(shiny::need(
        !is.null(network_data),
        paste0("Network data are unavailable for dataset '", input$dataset, "'.")
      ))
      network_data
    })

    selected_plot_data <- reactive({
      req(species_has_data(), nzchar(input$dataset))
      load_plotdata_fn(input$dataset)
    })

    # Enrichment term choices
    enrichment_terms <- reactive({
      req(nzchar(input$dataset), input$enrichment_category)
      data <- selected_network_data()
      format_enrichment_term_choices(
        data$cluster_descriptions, input$enrichment_category
      )
    })

    # When category changes, reset and populate the terms
    observeEvent(input$enrichment_category, {
      if (!nzchar(input$dataset)) return()
      data <- selected_network_data()

      updateCheckboxGroupInput(session, "selected_term", choices = NULL, selected = character(0))

      visNetworkProxy(ns("network")) %>%
        visUpdateNodes(nodes = data.frame(
          id = data$nodes_with_tooltips$id,
          borderWidth = 1,
          color = data$nodes_with_tooltips$color
        ))

      selected_cluster_genes(NULL)

      if (!is.null(input$enrichment_category) && input$enrichment_category != "") {
        term_choices <- enrichment_terms()
        if (length(term_choices) > 0) {
          updateCheckboxGroupInput(session, "selected_term", choices = term_choices, selected = character(0))
        } else {
          updateCheckboxGroupInput(session, "selected_term", choices = c("No terms available" = "none"), selected = character(0))
        }
      }
    }, ignoreNULL = FALSE, ignoreInit = TRUE)

    # Render network
    output$network <- renderVisNetwork({
      req(nzchar(input$dataset))
      data <- selected_network_data()

      known_ids <- data$nodes_with_tooltips$id
      safe_edges <- data$edges[data$edges$from %in% known_ids & data$edges$to %in% known_ids, ]

      # STRING combined scores are retained in the edge RDS files. Pass them
      # to visNetwork as edge values so widths are scaled within the selected
      # threshold-specific network. Older edge files remain usable.
      if ("combined_score" %in% names(safe_edges)) {
        safe_edges$value <- safe_edges$combined_score
      }

      nodes_positioned <- data$nodes_with_tooltips

      # Use pre-computed positions if available
      if (!is.null(data$positions)) {
        pos <- data$positions
        nodes_positioned$x <- pos$x[match(nodes_positioned$id, pos$id)]
        nodes_positioned$y <- pos$y[match(nodes_positioned$id, pos$id)]
      }

      visNetwork(nodes_positioned, safe_edges) %>%
        visEdges(scaling = list(min = 1, max = 6)) %>%
        visPhysics(enabled = FALSE) %>%
        visOptions(highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE)) %>%
        visInteraction(navigationButtons = TRUE) %>%
        visEvents(selectNode = sprintf(
          "function(nodes) { Shiny.setInputValue('%s', nodes.nodes[0], {priority: 'event'}); }",
          ns('selected_node')
        )) %>%
        visInteraction(tooltipStyle = "position: fixed; visibility: hidden; padding: 8px; font-family: Inter, Arial, sans-serif; font-size: 13px; color: #e0e3e5; background-color: rgba(15, 23, 42, 0.96); border: 1px solid #334155; border-radius: 6px; box-shadow: 0 12px 28px rgba(0,0,0,0.32); z-index: 1000; line-height: 1.4;") %>%
        visNodes(font = list(color = "#ecf0f1"))
    })

    # Zoom buttons
    observeEvent(input$zoom_in, {
      req(nzchar(input$dataset))
      visNetworkProxy(ns("network")) %>%
        visZoom(scale = 1.2, animation = list(duration = 500, easingFunction = "easeInOutQuad"))
    })
    observeEvent(input$zoom_out, {
      req(nzchar(input$dataset))
      visNetworkProxy(ns("network")) %>%
        visZoom(scale = 0.8, animation = list(duration = 500, easingFunction = "easeInOutQuad"))
    })
    # Highlight nodes when a term is selected
    observe({
      req(nzchar(input$dataset))
      selected_term <- input$selected_term
      data <- selected_network_data()

      visNetworkProxy(ns("network")) %>%
        visUpdateNodes(nodes = data.frame(
          id = data$nodes_with_tooltips$id, 
          borderWidth = 1, 
          color = data$nodes_with_tooltips$color
        ))

      if (length(selected_term) == 0) {
        selected_cluster_genes(NULL)
        return()
      }
      selected_term <- selected_term[1]
      if (selected_term == "none") {
        selected_cluster_genes(NULL)
        return()
      }

      req(data$cluster_to_genes, input$enrichment_category)
      selected_genes_df <- data$cluster_to_genes[
        data$cluster_to_genes$category == input$enrichment_category &
          data$cluster_to_genes$cluster == selected_term,
        , drop = FALSE
      ]
      if (nrow(selected_genes_df) == 0) {
        selected_cluster_genes(NULL)
        return()
      }

      selected_genes <- unique(unlist(selected_genes_df$genes))
      term_nodes <- data$nodes_with_tooltips$id[data$nodes_with_tooltips$id %in% selected_genes]
      all_nodes <- data$nodes_with_tooltips$id
      non_term_nodes <- setdiff(all_nodes, term_nodes)

      if (length(term_nodes) > 0) {
        visNetworkProxy(ns("network")) %>%
          visUpdateNodes(nodes = data.frame(id = term_nodes, borderWidth = 5)) %>%
          visUpdateNodes(nodes = data.frame(id = non_term_nodes, color = "rgba(200,200,200,0.5)"))
      }
      selected_cluster_genes(selected_genes)
    })

    # Data table
    proxy <- dataTableProxy("dataTable", session = session)
    output$dataTable <- renderDT({
      req(nzchar(input$dataset))
      data_info <- selected_plot_data()
      table_data <- data_info$table
      columnsToFormat <- data_info$columnsToFormat

      # Apply pVal / BH.Q cutoff filters
      pval_cols <- grep("^pVal", names(table_data), value = TRUE)
      bhq_cols  <- grep("^BH\\.Q", names(table_data), value = TRUE)
      if (length(pval_cols) > 0 && isTruthy(input$pval_cutoff) && input$pval_cutoff < 1) {
        keep <- apply(as.data.frame(table_data[pval_cols]), 1,
                      function(x) min(x, na.rm = TRUE)) <= input$pval_cutoff
        table_data <- table_data[keep, , drop = FALSE]
      }
      if (length(bhq_cols) > 0 && isTruthy(input$bhq_cutoff) && input$bhq_cutoff < 1) {
        keep <- apply(as.data.frame(table_data[bhq_cols]), 1,
                      function(x) min(x, na.rm = TRUE)) <= input$bhq_cutoff
        table_data <- table_data[keep, , drop = FALSE]
      }

      datatable(
        table_data,
        rownames = FALSE,
        selection = "multiple",     # make sure row selection is enabled
        options = list(
            pageLength = 5,
            autoWidth = TRUE,
            scrollX = TRUE,
            rowCallback = JS(
            "function(row, data, index) {",
            "  var columnsToFormat = ", toJSON(columnsToFormat), ";",
            "  columnsToFormat.forEach(function(i) {",
            "    if (!isNaN(data[i]) && data[i] !== null) {",
            "      var num = parseFloat(data[i]);",
            "      if (Math.abs(num) >= 1000 || Math.abs(num) < 0.001) {",
            "        $('td:eq(' + i + ')', row).html(num.toExponential(2));",
            "      } else {",
            "        $('td:eq(' + i + ')', row).html(num.toFixed(3));",
            "      }",
            "    }",
            "  });",
            "}"
            )
      ))
    })

    observeEvent(input$dataTable_rows_selected, {
      if (!nzchar(input$dataset)) return()
      selected_row <- input$dataTable_rows_selected
      if (length(selected_row) > 0) {
        data_info <- selected_plot_data()
        selected_gene <- data_info$table[selected_row, "symbol"]
        gene_list(c(selected_gene))
      }
    })
    observeEvent(input$clearSelection, {
      if (!nzchar(input$dataset)) return()
      # 1) Clear DT row selection
      selectRows(proxy, NULL)
      # 2) Clear plotted genes and the gene input box
      gene_list(NULL)
      updateTextInput(session, "gene", value = "")
      # 3) Clear enrichment selection (keep choices)
      updateCheckboxGroupInput(session, "selected_term", selected = character(0))
      selected_cluster_genes(NULL)
      # 4) Restore network node styling
      isolate({
        data <- selected_network_data()
        if (!is.null(data) && !is.null(data$nodes_with_tooltips)) {
          visNetworkProxy(ns("network")) %>%
            visUpdateNodes(nodes = data.frame(
              id = data$nodes_with_tooltips$id,
              borderWidth = 1,
              color = data$nodes_with_tooltips$color
            ))
          }
        })
      })

    # Network click -> plot
    observeEvent(input$selected_node, {
      if (!nzchar(input$dataset)) return()
      req(input$selected_node)

      node_id <- as.character(input$selected_node)
      data <- selected_network_data()
      nodes <- data$nodes_with_tooltips

      # Map to a gene symbol: prefer 'symbol', then 'label', else use id
      symbol <- NULL
      if (!is.null(nodes)) {
        if ("symbol" %in% names(nodes)) {
          hit <- nodes$symbol[nodes$id == node_id]
          if (length(hit) && !is.na(hit[1])) symbol <- as.character(hit[1])
        }
        if (is.null(symbol) && "label" %in% names(nodes)) {
          hit <- nodes$label[nodes$id == node_id]
          if (length(hit) && !is.na(hit[1])) symbol <- as.character(hit[1])
        }
      }
      if (is.null(symbol) || !nzchar(symbol)) symbol <- node_id

      # Set the plotted genes and reflect in the text box
      gene_list(c(symbol))
      updateTextInput(session, "gene", value = symbol)
    })

    # Plot button
    formatGeneInput <- function(str, species) {
      genes <- strsplit(str, ",")[[1]]
      genes <- trimws(genes)

      if (identical(species, "baboon")) {
        return(toupper(genes))
      }

      if (identical(species, "mouse")) {
        return(vapply(genes, function(gene) {
          paste(toupper(substring(gene, 1, 1)), tolower(substring(gene, 2)), sep = "")
        }, character(1), USE.NAMES = FALSE))
      }

      genes
    }
    observeEvent(input$plot, {
      if (!nzchar(input$dataset)) return()
      if (!is.null(input$gene) && nchar(input$gene) > 0) {
        genes <- formatGeneInput(input$gene, input$species)
        gene_list(genes)
        updateTextInput(session, "gene", value = paste(genes, collapse = ", "))
      }
    })

    # Plot selected genes from enrichment
    observeEvent(input$plot_enrichment, {
      if (!nzchar(input$dataset)) return()
      req(selected_cluster_genes())
      gene_list(selected_cluster_genes())
    })

    # Reserve stable plot slots, then let the renderers fill them in batches.
    observe({
      genes <- gene_list()
      has_genes <- !is.null(genes) && length(genes) > 0 && isTRUE(nzchar(input$dataset))
      loaded <- if (has_genes) min(plot_page_size, length(genes)) else 0L
      plot_slot_count(loaded)
      visible_plot_count(if (loaded > 0) min(plot_batch_size, loaded) else 0L)
    })

    observeEvent(
      list(input$show_fit, input$show_line, input$show_comparison),
      {
        genes <- gene_list()
        if (is.null(genes) || length(genes) == 0 || !isTRUE(nzchar(input$dataset))) return()
        loaded <- min(length(genes), plot_slot_count())
        visible_plot_count(if (loaded > 0) min(plot_batch_size, loaded) else 0L)
      },
      ignoreInit = TRUE
    )

    observeEvent(input$load_more_plots, {
      genes <- gene_list()
      if (is.null(genes) || length(genes) == 0 || !isTRUE(nzchar(input$dataset))) return()

      current_loaded <- min(length(genes), plot_slot_count())
      next_loaded <- min(length(genes), current_loaded + plot_page_size)
      if (next_loaded <= current_loaded) return()

      plot_slot_count(next_loaded)
      visible_plot_count(min(next_loaded, max(visible_plot_count(), current_loaded) + plot_batch_size))
    }, ignoreInit = TRUE)

    observe({
      genes <- gene_list()
      if (is.null(genes) || length(genes) == 0 || !isTRUE(nzchar(input$dataset))) return()

      shown <- visible_plot_count()
      total <- min(length(genes), plot_slot_count())
      if (shown <= 0 || shown >= total) return()

      genes_snapshot <- as.character(genes)
      session$onFlushed(function() {
        isolate({
          current_genes <- gene_list()
          if (!is.null(current_genes) && identical(as.character(current_genes), genes_snapshot)) {
            current_loaded <- min(length(current_genes), plot_slot_count())
            visible_plot_count(min(current_loaded, shown + plot_batch_size))
          }
        })
      }, once = TRUE)
    })

    # Plot UI and renderers
    output$plotsUI <- renderUI({
      req(nzchar(input$dataset), gene_list())
      genes <- gene_list()
      num_genes <- length(genes)
      slot_count <- min(num_genes, plot_slot_count())
      visible_count <- min(num_genes, visible_plot_count())
      if (slot_count <= 0) return(NULL)

      num_cols <- 4
      num_rows <- ceiling(slot_count / num_cols)

      plot_output_list <- lapply(seq_len(slot_count), function(i) {
        plotname <- paste0("plot", i)
        plotOutput(ns(plotname), height = "260px", width = "100%")
      })

      rows <- lapply(seq_len(num_rows), function(row) {
        fluidRow(class = "plot-row",
                 lapply(1:num_cols, function(col) {
                   idx <- (row - 1) * num_cols + col
                   if (idx <= slot_count) column(3, class = "plot-column", plot_output_list[[idx]])
                 }))
      })

      loading_status <- if (visible_count < slot_count) {
        div(
          class = "plot-loading-status",
          sprintf("Rendering plots %d of %d loaded...", visible_count, slot_count)
        )
      } else {
        NULL
      }

      footer <- div(
        class = "plot-pagination-footer",
        div(
          class = "plot-pagination-status",
          sprintf("Showing %d out of %d genes", slot_count, num_genes)
        ),
        if (slot_count < num_genes) {
          actionButton(ns("load_more_plots"), "Load more", class = "btn btn-default btn-sm")
        } else {
          NULL
        }
      )

      tagList(loading_status, do.call(tagList, rows), footer)
    })

    observe({
      req(nzchar(input$dataset), gene_list())
      genes <- gene_list()
      slot_count <- min(length(genes), plot_slot_count())
      if (slot_count <= 0) return()
      data_info <- selected_plot_data()

      lapply(seq_len(slot_count), function(i) {
        local({
          plot_index <- i
          gene <- genes[[i]]

          output[[paste0("plot", plot_index)]] <- renderCachedPlot({
            req(plot_index <= visible_plot_count(), cancelOutput = TRUE)
            plots <- generate_plot_fn(data_info, gene, input, session)
            plots[[gene]]
          }, cacheKeyExpr = {
            req(plot_index <= visible_plot_count(), cancelOutput = TRUE)
            list(
              scope           = "explorer-expression",
              dataset         = input$dataset,
              gene            = gene,
              show_fit        = input$show_fit,
              show_line       = input$show_line,
              show_comparison = input$show_comparison,
              data_signature  = plot_data_cache_signature(data_info)
            )
          }, cache = "session")
        })
      })
    })
  })
}
