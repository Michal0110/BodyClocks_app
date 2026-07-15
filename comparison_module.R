# comparison_module.R
# Shiny module for the Dataset Comparison tab.
# Provides comparisonUI() and comparisonServer().

# ============================================================
# UI
# ============================================================
comparisonUI <- function(id,
                         species_choices,
                         tissue_choices_by_species,
                         dataset_choices_by_species) {
  ns <- NS(id)

  # Default initial choices derived from first species's first tissue
  default_species  <- species_choices[[1]]
  default_tissues  <- tissue_choices_by_species[[default_species]]
  init_choices1    <- if (length(default_tissues) > 0)
                        dataset_choices_by_species[[default_species]][[default_tissues[[1]]]]
                      else character(0)
  init_choices2    <- init_choices1

  fluidPage(
    fluidRow(
      # ---- Sidebar ----
      column(
        width = 3,
        div(
          class = "explorer-sidebar bc-sidebar bc-panel",

          # --- Species selector ---
          div(class = "section-header", "Species:"),
          selectInput(
            inputId  = ns("species"),
            label    = NULL,
            choices  = species_choices,
            selected = default_species
          ),

          # --- Dataset 1 ---
          div(class = "section-header", "Dataset 1"),
          selectInput(
            inputId  = ns("tissue1"),
            label    = NULL,
            choices  = default_tissues,
            selected = if (length(default_tissues) > 0) default_tissues[[1]] else ""
          ),
          selectInput(
            inputId  = ns("dataset1"),
            label    = NULL,
            choices  = init_choices1,
            selected = init_choices1[[1]]
          ),

          # --- Dataset 2 ---
          div(class = "section-header", "Dataset 2"),
          selectInput(
            inputId  = ns("tissue2"),
            label    = NULL,
            choices  = default_tissues,
            selected = if (length(default_tissues) > 0) default_tissues[[1]] else ""
          ),
          selectInput(
            inputId  = ns("dataset2"),
            label    = NULL,
            choices  = init_choices2,
            selected = init_choices2[[1]]
          ),

          # --- BH.Q cutoff ---
          div(class = "section-header", "BH.Q cutoff"),
          sliderInput(
            inputId = ns("bhq_cutoff"),
            label   = NULL,
            min     = 0.001, max = 0.2,
            value   = 0.05,
            step    = 0.001,
            ticks   = FALSE
          ),

          # --- Compare button ---
          actionButton(
            inputId = ns("compare"),
            label   = "Compare",
            class   = "btn btn-primary btn-block action-button"
          ),

          hr(),

          # --- Overlapping genes table (shown after Compare) ---
          conditionalPanel(
            condition = paste0("output['", ns("has_result"), "'] === true"),

            # --- Enrichment overlap section ---
            div(class = "section-header", "Enrichment Overlap"),
            selectInput(
              inputId  = ns("enrichment_category"),
              label    = NULL,
              choices  = c("Select category" = "",
                           "GO Biological Process" = "Process",
                           "KEGG Pathway" = "KEGG",
                           "Reactome" = "RCTM",
                           "WikiPathways" = "WikiPathways"),
              selected = ""
            ),
            div(
              class = "data-table-container",
              DT::DTOutput(ns("enrichment_table"))
            ),
            div(
              class = "bc-action-row",
              actionButton(ns("enrich_plot_all"), "Plot All Genes",
                           class = "btn btn-primary btn-sm bc-action-button"),
              actionButton(ns("enrich_plot_common"), "Plot Common Rhythmic",
                           class = "btn btn-default btn-sm bc-action-button")
            ),

            hr(),

            div(
              class = "data-table-container",
              div(class = "section-header data-table-title", "Overlapping Genes"),
              DT::DTOutput(ns("overlap_table")),
              fluidRow(
                column(6, numericInput(ns("overlap_pval_cutoff"), "pVal \u2264",
                                       value = 1, min = 0, max = 1, step = 0.001)),
                column(6, numericInput(ns("overlap_bhq_cutoff"), "BH.Q \u2264",
                                       value = 1, min = 0, max = 1, step = 0.001))
              )
            ),
            br(),
            actionButton(
              inputId = ns("clear"),
              label   = "Clear",
              class   = "btn btn-default btn-sm btn-block"
            )
          )
        )
      ),

      # ---- Main panel ----
      column(
        width = 9,
        div(
          class = "explorer-main",
          uiOutput(ns("main_content"))
        )
      )
    )
  )
}

# ============================================================
# Server
# ============================================================
comparisonServer <- function(id,
                              species_choices,
                              tissue_choices_by_species,
                              dataset_choices_by_species,
                              load_network_fn  = NULL,
                              load_plotdata_fn,
                              generate_plot_fn = NULL) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # ---- Reactive state ----
    comp_result <- reactiveVal(NULL)   # full comparison output
    gene_list   <- reactiveVal(NULL)   # genes selected for expression plots
    plot_batch_size <- 4L
    plot_page_size <- 32L
    plot_slot_count <- reactiveVal(0L)
    visible_plot_count <- reactiveVal(0L)

    # Flag output so conditionalPanel can query it
    output$has_result <- reactive({ !is.null(comp_result()) })
    outputOptions(output, "has_result", suspendWhenHidden = FALSE)

    observe({
      gl <- gene_list()
      cr <- comp_result()
      has_genes <- !is.null(gl) && length(gl) > 0 && !is.null(cr)
      loaded <- if (has_genes) min(plot_page_size, length(gl)) else 0L
      plot_slot_count(loaded)
      visible_plot_count(if (loaded > 0) min(plot_batch_size, loaded) else 0L)
    })

    observeEvent(input$load_more_plots, {
      gl <- gene_list()
      cr <- comp_result()
      if (is.null(gl) || length(gl) == 0 || is.null(cr)) return()

      current_loaded <- min(length(gl), plot_slot_count())
      next_loaded <- min(length(gl), current_loaded + plot_page_size)
      if (next_loaded <= current_loaded) return()

      plot_slot_count(next_loaded)
      visible_plot_count(min(next_loaded, max(visible_plot_count(), current_loaded) + plot_batch_size))
    })

    observe({
      gl <- gene_list()
      cr <- comp_result()
      if (is.null(gl) || length(gl) == 0 || is.null(cr)) return()

      shown <- visible_plot_count()
      total <- min(length(gl), plot_slot_count())
      if (shown <= 0 || shown >= total) return()

      genes_snapshot <- as.character(gl)
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

    # ---- Species change: reset both tissue and dataset dropdowns ----
    observeEvent(input$species, {
      tc <- tissue_choices_by_species[[input$species]]
      if (is.null(tc) || length(tc) == 0) tc <- character(0)
      updateSelectInput(session, "tissue1", choices = tc,
                        selected = if (length(tc) > 0) tc[[1]] else "")
      updateSelectInput(session, "tissue2", choices = tc,
                        selected = if (length(tc) > 0) tc[[1]] else "")
      comp_result(NULL)
      gene_list(NULL)
    }, ignoreInit = TRUE)

    # ---- Tissue → dataset cascades ----
    observeEvent(input$tissue1, {
      choices <- dataset_choices_by_species[[input$species]][[input$tissue1]]
      if (is.null(choices)) choices <- character(0)
      updateSelectInput(session, "dataset1",
                        choices  = choices,
                        selected = if (length(choices) > 0) choices[[1]] else "")
    }, ignoreInit = TRUE)

    observeEvent(input$tissue2, {
      choices <- dataset_choices_by_species[[input$species]][[input$tissue2]]
      if (is.null(choices)) choices <- character(0)
      updateSelectInput(session, "dataset2",
                        choices  = choices,
                        selected = if (length(choices) > 0) choices[[1]] else "")
    }, ignoreInit = TRUE)

    # ---- Compare button ----
    observeEvent(input$compare, {
      req(input$dataset1, input$dataset2)

      withProgress(message = "Running comparison…", value = 0, {

        setProgress(0.1, message = "Loading dataset 1…")
        plotdata1 <- tryCatch(
          load_plotdata_fn(input$dataset1),
          error = function(e) {
            showNotification(paste("Error loading", input$dataset1, ":", e$message),
                             type = "error")
            NULL
          }
        )
        req(plotdata1)

        setProgress(0.3, message = "Loading dataset 2…")
        plotdata2 <- tryCatch(
          load_plotdata_fn(input$dataset2),
          error = function(e) {
            showNotification(paste("Error loading", input$dataset2, ":", e$message),
                             type = "error")
            NULL
          }
        )
        req(plotdata2)

        setProgress(0.5, message = "Computing gene overlap…")

        gd1_all <- extract_gene_data(plotdata1$table, input$dataset1)
        gd2_all <- extract_gene_data(plotdata2$table, input$dataset2)

        if (is.null(gd1_all) || is.null(gd2_all)) {
          showNotification("Could not parse gene data for one or both datasets.",
                           type = "error")
          return()
        }

        # Filter by BH.Q cutoff for the comparison itself
        gd1 <- dplyr::filter(gd1_all, !is.na(bhq) & bhq <= input$bhq_cutoff)
        gd2 <- dplyr::filter(gd2_all, !is.na(bhq) & bhq <= input$bhq_cutoff)

        # Hypergeometric test.
        # Universe = genes measured in BOTH datasets; genes absent from one
        # dataset cannot be rhythmic in it and must not inflate the background.
        shared_genes <- intersect(gd1_all$symbol, gd2_all$symbol)
        genome_size  <- length(shared_genes)
        ht <- hypergeometric_test(
          intersect(gd1$symbol, shared_genes),
          intersect(gd2$symbol, shared_genes),
          genome_size
        )

        # Overlap table (BH.Q + phase from both datasets)
        overlap_syms <- intersect(gd1$symbol, gd2$symbol)
        overlap_tbl <- dplyr::inner_join(
          dplyr::filter(gd1, symbol %in% overlap_syms) %>%
            dplyr::rename(pval_ds1 = pval, bhq_ds1 = bhq, phase_ds1 = phase),
          dplyr::filter(gd2, symbol %in% overlap_syms) %>%
            dplyr::rename(pval_ds2 = pval, bhq_ds2 = bhq, phase_ds2 = phase),
          by = "symbol"
        ) %>% dplyr::arrange(bhq_ds1)

        setProgress(0.7, message = "Phase coherence analysis…")

        # Locate clock anchor gene phases: Bmal1 (mouse) or ARNTL (baboon)
        find_clock_row <- function(gd) {
          r <- gd[gd$symbol == "Bmal1" & !is.na(gd$phase), ]
          if (nrow(r) == 0) r <- gd[gd$symbol == "ARNTL" & !is.na(gd$phase), ]
          r
        }
        clock_gene_name <- function(gd) {
          if (any(gd$symbol == "Bmal1")) "Bmal1" else "ARNTL"
        }

        bmal1_r1 <- find_clock_row(gd1_all)
        bmal1_r2 <- find_clock_row(gd2_all)

        phase_info <- if (nrow(bmal1_r1) == 0) {
          list(available = FALSE,
               message = paste0("Clock-relative phase analysis unavailable: ",
                                clock_gene_name(gd1_all), " not detected in dataset 1 (",
                                input$dataset1, ")."))
        } else if (nrow(bmal1_r2) == 0) {
          list(available = FALSE,
               message = paste0("Clock-relative phase analysis unavailable: ",
                                clock_gene_name(gd2_all), " not detected in dataset 2 (",
                                input$dataset2, ")."))
        } else if (nrow(overlap_tbl) < 3) {
          list(available = FALSE,
               message = "Clock-relative phase analysis unavailable: fewer than 3 overlapping rhythmic genes.")
        } else {
          crp_df <- compute_crp_for_overlap(
            overlap_tbl %>% dplyr::rename(phase1 = phase_ds1, phase2 = phase_ds2),
            bmal1_phase1 = bmal1_r1$phase[1],
            bmal1_phase2 = bmal1_r2$phase[1]
          )
          if (is.null(crp_df) || nrow(crp_df) < 3) {
            list(available = FALSE,
                 message = "Clock-relative phase analysis unavailable: insufficient phase data in overlapping genes.")
          } else {
            perm <- run_crp_permutation_test(crp_df$crp1, crp_df$crp2)
            list(
              available    = TRUE,
              crp_df       = crp_df,
              perm         = perm,
              bmal1_phase1 = bmal1_r1$phase[1],
              bmal1_phase2 = bmal1_r2$phase[1]
            )
          }
        }

        gene_list(NULL)

        # Load enrichment data if network loader is available
        enrich1 <- NULL
        enrich2 <- NULL
        if (!is.null(load_network_fn)) {
          # Dataset Comparison retains the historical all-STRING enrichment
          # definition; the background toggle is intentionally Explorer-only.
          net1 <- tryCatch(
            load_network_fn(
              input$dataset1, enrichment_background = "all_string"
            ),
            error = function(e) NULL
          )
          net2 <- tryCatch(
            load_network_fn(
              input$dataset2, enrichment_background = "all_string"
            ),
            error = function(e) NULL
          )
          if (!is.null(net1)) {
            enrich1 <- list(desc = net1$cluster_descriptions, genes = net1$cluster_to_genes)
          }
          if (!is.null(net2)) {
            enrich2 <- list(desc = net2$cluster_descriptions, genes = net2$cluster_to_genes)
          }
        }

        comp_result(list(
          ht          = ht,
          overlap_tbl = overlap_tbl,
          plotdata1   = plotdata1,
          plotdata2   = plotdata2,
          ds1         = input$dataset1,
          ds2         = input$dataset2,
          phase       = phase_info,
          enrich1     = enrich1,
          enrich2     = enrich2
        ))

        setProgress(1, message = "Done")
      })
    })

    # ---- Overlap DT ----
    filtered_overlap_table <- reactive({
      cr <- comp_result()
      req(cr)

      df <- cr$overlap_tbl

      min_or_inf <- function(x) {
        x <- suppressWarnings(as.numeric(x))
        if (all(is.na(x))) return(Inf)
        min(x, na.rm = TRUE)
      }

      pval_cols <- intersect(c("pval_ds1", "pval_ds2"), names(df))
      if (length(pval_cols) > 0 &&
          isTruthy(input$overlap_pval_cutoff) &&
          input$overlap_pval_cutoff < 1) {
        keep <- apply(as.data.frame(df[pval_cols]), 1, min_or_inf) <= input$overlap_pval_cutoff
        df <- df[keep, , drop = FALSE]
      }

      bhq_cols <- intersect(c("bhq_ds1", "bhq_ds2"), names(df))
      if (length(bhq_cols) > 0 &&
          isTruthy(input$overlap_bhq_cutoff) &&
          input$overlap_bhq_cutoff < 1) {
        keep <- apply(as.data.frame(df[bhq_cols]), 1, min_or_inf) <= input$overlap_bhq_cutoff
        df <- df[keep, , drop = FALSE]
      }

      df
    })

    output$overlap_table <- DT::renderDT({
      cr <- comp_result()
      req(cr)
      df     <- filtered_overlap_table()
      label1 <- get_dataset_label(cr$ds1)
      label2 <- get_dataset_label(cr$ds2)

      tbl <- data.frame(
        Symbol  = df$symbol,
        pVal_1  = df$pval_ds1,
        BH.Q_1  = df$bhq_ds1,
        Phase_1 = round(df$phase_ds1, 1),
        pVal_2  = df$pval_ds2,
        BH.Q_2  = df$bhq_ds2,
        Phase_2 = round(df$phase_ds2, 1),
        stringsAsFactors = FALSE
      )

      DT::datatable(
        tbl,
        colnames  = c("Symbol",
                      paste0("pVal (", label1, ")"),
                      paste0("BH.Q (", label1, ")"),
                      paste0("Phase (", label1, ")"),
                      paste0("pVal (", label2, ")"),
                      paste0("BH.Q (", label2, ")"),
                      paste0("Phase (", label2, ")")),
        selection = "multiple",
        rownames  = FALSE,
        options   = list(
          pageLength = 5,
          autoWidth  = TRUE,
          scrollX    = TRUE,
          rowCallback = DT::JS(
            "function(row, data, index) {",
            "  [1, 2, 3, 4, 5, 6].forEach(function(i) {",
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
        )
      )
    })

    # Table row click → gene list
    observeEvent(input$overlap_table_rows_selected, {
      cr <- comp_result()
      req(cr)
      df    <- filtered_overlap_table()
      idx   <- input$overlap_table_rows_selected
      if (length(idx) == 0) return()
      genes <- df$symbol[idx]
      gene_list(genes)
      DT::dataTableProxy("enrichment_table", session = session) %>% DT::selectRows(NULL)
    })

    # Clear selection and plots
    observeEvent(input$clear, {
      gene_list(NULL)
      DT::dataTableProxy("overlap_table", session = session) %>% DT::selectRows(NULL)
      DT::dataTableProxy("enrichment_table", session = session) %>% DT::selectRows(NULL)
    })

    # ---- Plotly scatter click → replace gene list ----
    observe({
      ev <- plotly::event_data("plotly_click", source = ns("phase_scatter"))
      if (is.null(ev)) return()
      cd <- ev$customdata
      if (!is.null(cd) && length(cd) > 0 && nzchar(cd[1])) {
        new_genes <- strsplit(cd[1], ",")[[1]]
        gene_list(new_genes)
        DT::dataTableProxy("overlap_table", session = session) %>% DT::selectRows(NULL)
        DT::dataTableProxy("enrichment_table", session = session) %>% DT::selectRows(NULL)
      }
    })

    # ---- Enrichment overlap reactive ----
    enrichment_overlap <- reactive({
      cr  <- comp_result()
      cat <- input$enrichment_category
      if (is.null(cr) || is.null(cat) || cat == "") return(NULL)
      if (is.null(cr$enrich1) || is.null(cr$enrich2)) return(NULL)
      overlap_syms <- if (!is.null(cr$overlap_tbl)) cr$overlap_tbl$symbol else character(0)
      compute_enrichment_overlap(
        cr$enrich1$desc, cr$enrich1$genes,
        cr$enrich2$desc, cr$enrich2$genes,
        cat,
        overlap_syms
      )
    })

    # ---- Enrichment DT ----
    output$enrichment_table <- DT::renderDT({
      eo <- enrichment_overlap()
      if (is.null(eo) || nrow(eo) == 0) return(NULL)

      cr     <- comp_result()
      label1 <- get_dataset_label(cr$ds1)
      label2 <- get_dataset_label(cr$ds2)

      tbl <- data.frame(
        Term       = eo$term,
        Genes_DS1  = eo$n_genes_ds1,
        Genes_DS2  = eo$n_genes_ds2,
        Common     = eo$n_common,
        stringsAsFactors = FALSE
      )

      DT::datatable(
        tbl,
        colnames  = c("Term",
                       paste0("# genes (", label1, ")"),
                       paste0("# genes (", label2, ")"),
                       "# common rhythmic"),
        selection = "single",
        rownames  = FALSE,
        options   = list(
          pageLength = 5,
          autoWidth  = TRUE,
          scrollX    = TRUE,
          dom        = "ltp"
        )
      ) %>%
        DT::formatStyle(columns = 1:4, fontSize = "13px")
    })

    # ---- Enrichment "Plot All Genes" button ----
    observeEvent(input$enrich_plot_all, {
      eo  <- enrichment_overlap()
      idx <- input$enrichment_table_rows_selected
      if (is.null(eo) || is.null(idx) || length(idx) == 0) {
        showNotification("Select an enrichment term first.", type = "warning")
        return()
      }
      row <- eo[idx, ]
      g1 <- if (nzchar(row$genes_ds1)) strsplit(row$genes_ds1, ",")[[1]] else character(0)
      g2 <- if (nzchar(row$genes_ds2)) strsplit(row$genes_ds2, ",")[[1]] else character(0)
      all_genes <- unique(c(g1, g2))
      if (length(all_genes) == 0) return()
      gene_list(all_genes)
      DT::dataTableProxy("overlap_table", session = session) %>% DT::selectRows(NULL)
    })

    # ---- Enrichment "Plot Common Rhythmic" button ----
    observeEvent(input$enrich_plot_common, {
      eo  <- enrichment_overlap()
      idx <- input$enrichment_table_rows_selected
      if (is.null(eo) || is.null(idx) || length(idx) == 0) {
        showNotification("Select an enrichment term first.", type = "warning")
        return()
      }
      row <- eo[idx, ]
      common <- if (nzchar(row$common_genes)) strsplit(row$common_genes, ",")[[1]] else character(0)
      if (length(common) == 0) {
        showNotification("No common rhythmic genes for this term.", type = "warning")
        return()
      }
      gene_list(common)
      DT::dataTableProxy("overlap_table", session = session) %>% DT::selectRows(NULL)
    })

    # ---- Rose plot ----
    output$rose_plot <- renderPlot({
      cr <- comp_result()
      req(cr, nrow(cr$overlap_tbl) > 0)

      label1 <- get_dataset_label(cr$ds1)
      label2 <- get_dataset_label(cr$ds2)
      df     <- cr$overlap_tbl

      # Add +1 offset so 2h bins centre on even CT hours and CT0 bar
      # doesn't straddle the polar origin (x = 0)
      phase1_plot <- (df$phase_ds1[!is.na(df$phase_ds1)] + 1) %% 24
      phase2_plot <- (df$phase_ds2[!is.na(df$phase_ds2)] + 1) %% 24

      breaks <- seq(0, 24, by = 2)
      bins1  <- cut(phase1_plot, breaks = breaks, right = FALSE, include.lowest = TRUE)
      bins2  <- cut(phase2_plot, breaks = breaks, right = FALSE, include.lowest = TRUE)

      counts1 <- table(factor(bins1, levels = levels(bins1)))
      counts2 <- table(factor(bins2, levels = levels(bins2)))

      # Midpoints in offset space: 1, 3, 5, ..., 23
      bin_mids_plot <- seq(1, 23, by = 2)

      rose_df <- data.frame(
        x_plot  = rep(bin_mids_plot, 2),
        count   = c(as.integer(counts1), as.integer(counts2)),
        dataset = factor(rep(c(label1, label2), each = length(bin_mids_plot)),
                         levels = c(label1, label2)),
        stringsAsFactors = FALSE
      )

      # CT axis labels: offset positions for CT0, CT4, CT8, CT12, CT16, CT20
      ct_breaks <- seq(1, 21, by = 4)  # 1=CT0, 5=CT4, 9=CT8, 13=CT12, 17=CT16, 21=CT20
      ct_labels <- paste0("CT", seq(0, 20, by = 4))

      ggplot2::ggplot(rose_df,
             ggplot2::aes(x = x_plot, y = count, fill = dataset)) +
        ggplot2::geom_bar(stat = "identity", alpha = 0.6, colour = "black",
                          width = 1, position = "identity") +
        # start = -pi/12: compensates for the +1 offset so CT0 sits at 12 o'clock
        ggplot2::coord_polar(theta = "x", start = -pi / 12) +
        ggplot2::scale_x_continuous(breaks = ct_breaks, labels = ct_labels,
                                    limits = c(0, 24), expand = c(0, 0)) +
        ggplot2::scale_fill_manual(values = stats::setNames(
          c("#1565C0", "#BF360C"), c(label1, label2))) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          plot.title       = ggplot2::element_text(size = 11, face = "bold", hjust = 0.5),
          axis.text.x      = ggplot2::element_text(size = 8, face = "bold"),
          axis.text.y      = ggplot2::element_blank(),
          axis.title        = ggplot2::element_blank(),
          panel.grid.major  = ggplot2::element_line(colour = "grey80", linewidth = 0.4),
          panel.grid.minor  = ggplot2::element_blank(),
          legend.position   = "bottom",
          legend.title      = ggplot2::element_blank()
        ) +
        ggplot2::labs(title = "Phase distribution\n(overlapping rhythmic genes)")
    })

    # ---- CRP bubble scatter (plotly, interactive) ----
    output$crp_scatter <- plotly::renderPlotly({
      cr <- comp_result()
      req(cr)
      pi <- cr$phase
      req(pi$available)

      crp_df <- pi$crp_df
      label1 <- get_dataset_label(cr$ds1)
      label2 <- get_dataset_label(cr$ds2)

      # Aggregate genes sharing the same (crp1, crp2) phase intersection
      bubble_df <- crp_df %>%
        dplyr::group_by(crp1, crp2) %>%
        dplyr::summarise(
          count  = dplyr::n(),
          genes  = paste(symbol, collapse = ","),
          n_cc   = sum(classification == "Clock-controlled"),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          # Scale bubble area: sqrt so area is proportional to count
          bsize  = pmax(sqrt(count) * 8, 6),
          bcolor = ifelse(n_cc > count / 2, "steelblue", "#9E9E9E"),
          hover  = paste0(
            "<b>", count, " gene", ifelse(count > 1, "s", ""), "</b><br>",
            "CRP (", label1, "): ", round(crp1, 1), " h<br>",
            "CRP (", label2, "): ", round(crp2, 1), " h<br>",
            "<i>", ifelse(nchar(genes) > 80,
                         paste0(substr(genes, 1, 80), "\u2026"),
                         genes), "</i>"
          )
        )

      plotly::plot_ly(
        data   = bubble_df,
        source = ns("phase_scatter")
      ) %>%
        plotly::add_trace(
          x            = ~crp1,
          y            = ~crp2,
          type         = "scatter",
          mode         = "markers",
          customdata   = ~genes,
          text         = ~hover,
          hoverinfo    = "text",
          marker       = list(
            size    = bubble_df$bsize,
            color   = bubble_df$bcolor,
            opacity = 0.7,
            line    = list(color = "white", width = 0.5)
          ),
          showlegend = FALSE
        ) %>%
        plotly::layout(
          autosize = TRUE,
          xaxis = list(
            title         = paste0("CRP in ", label1, " (h from Bmal1)"),
            range         = c(-12.5, 16),
            tickvals      = seq(-12, 12, by = 4),
            ticktext      = as.character(seq(-12, 12, by = 4)),
            zeroline      = TRUE,
            zerolinewidth = 1,
            zerolinecolor = "#aaa"
          ),
          yaxis = list(
            title         = paste0("CRP in ", label2, " (h from Bmal1)"),
            range         = c(-12.5, 16),
            tickvals      = seq(-12, 12, by = 4),
            ticktext      = as.character(seq(-12, 12, by = 4)),
            zeroline      = TRUE,
            zerolinewidth = 1,
            zerolinecolor = "#aaa"
          ),
          shapes = list(
            list(
              type = "line",
              x0 = -12, x1 = 12, y0 = -12, y1 = 12,
              line = list(color = "red", dash = "dash", width = 1)
            )
          ),
          margin = list(l = 60, r = 20, t = 40, b = 70)
        ) %>%
        plotly::config(displayModeBar = FALSE) %>%
        plotly::event_register("plotly_click")
    })

    # ---- Main content (rendered on Compare) ----
    output$main_content <- renderUI({
      cr <- comp_result()

      if (is.null(cr)) {
        return(
          div(
            class = "bc-empty-state bc-empty-state-large",
            tags$i(class = "fa fa-arrow-left bc-empty-icon"),
            h4("Select two datasets in the sidebar and click Compare"),
            p("Gene overlap statistics, clock-relative phase coherence, and",
              "interactive expression plots will appear here.")
          )
        )
      }

      ht     <- cr$ht
      label1 <- get_dataset_label(cr$ds1)
      label2 <- get_dataset_label(cr$ds2)

      # Hypergeometric p-value formatting
      p_fmt <- if (is.nan(ht$p_value) || ht$p_value < 0.001)
        "< 0.001"
      else
        formatC(ht$p_value, format = "e", digits = 2)

      # ---- Overlap stats box ----
      overlap_box <- div(
        class = "comp-stats-panel",
        h5("Gene overlap", class = "comp-panel-title"),
        p(class = "comp-panel-note",
          "Rhythmic genes shared between both datasets.",
          "p-value from hypergeometric test."),
        div(
          class = "comp-stats-row",
          div(class = "comp-stat",
            div(class = "comp-stat-value", ht$overlap),
            div(class = "comp-stat-label", "overlapping genes")
          ),
          div(class = "comp-stat",
            div(class = "comp-stat-value", ht$size1),
            div(class = "comp-stat-label", paste0("rhythmic in ", label1))
          ),
          div(class = "comp-stat",
            div(class = "comp-stat-value", ht$size2),
            div(class = "comp-stat-label", paste0("rhythmic in ", label2))
          ),
          div(class = "comp-stat",
            div(class = "comp-stat-value", round(ht$expected, 1)),
            div(class = "comp-stat-label", "expected by chance")
          ),
          div(class = "comp-stat",
            div(class = "comp-stat-value",
                if (!is.na(ht$fold_enrichment)) round(ht$fold_enrichment, 2) else "\u2014"),
            div(class = "comp-stat-label", "fold enrichment")
          ),
          div(class = "comp-stat",
            div(class = "comp-stat-value", p_fmt),
            div(class = "comp-stat-label", "p (hypergeometric)")
          )
        )
      )

      # ---- Phase coherence section ----
      phase_section <- if (!cr$phase$available) {
        div(
          class = "alert alert-info bc-alert",
          tags$i(class = "fa fa-info-circle"), " ", cr$phase$message
        )
      } else {
        perm   <- cr$phase$perm
        p_perm <- if (!is.na(perm$p_value) && perm$p_value == 0)
          paste0("< ", formatC(1 / 1000, format = "e", digits = 1))
        else if (is.na(perm$p_value))
          "NA"
        else
          formatC(perm$p_value, format = "e", digits = 2)

        tagList(
          div(
            class = "comp-stats-panel",
            h5("Clock-relative phase coherence",
               class = "comp-panel-title"),
            p(class = "comp-panel-note",
              "Genes classified as clock-controlled if their clock-relative phases",
              "agree within 4 h. p-value from permutation test (n = 1000, seed = 42)."),
            div(
              class = "comp-stats-row",
              div(class = "comp-stat",
                div(class = "comp-stat-value", perm$n_coherent),
                div(class = "comp-stat-label", "n clock-controlled")
              ),
              div(class = "comp-stat",
                div(class = "comp-stat-value", paste0(round(perm$pct_coherent, 1), "%")),
                div(class = "comp-stat-label", "% clock-controlled")
              ),
              div(class = "comp-stat",
                div(class = "comp-stat-value", paste0(round(perm$pct_null, 1), "%")),
                div(class = "comp-stat-label", "null (permutation) %")
              ),
              div(class = "comp-stat",
                div(class = "comp-stat-value", p_perm),
                div(class = "comp-stat-label", "p-value (permutation)")
              )
            ),
            div(
              class = "comp-phase-layout",
              div(
                class = "comp-rose-pane",
                plotOutput(ns("rose_plot"), height = "480px", width = "100%")
              ),
              div(
                class = "comp-scatter-pane",
                div(
                  class = "comp-scatter-wrap",
                  plotly::plotlyOutput(ns("crp_scatter"), height = "480px", width = "100%")
                ),
                div(
                  class = "comp-interaction-note",
                  tags$i(class = "fa fa-mouse-pointer"),
                  " Click a point to add the gene to expression plots below.",
                  " Genes can also be selected from the table on the left."
                )
              )
            )
          )
        )
      }

      # ---- Expression plots (shown when genes selected) ----
      genes_section <- if (!is.null(gene_list()) && length(gene_list()) > 0) {
        div(
          class = "plots-container bc-plots-panel bc-panel comp-expression-section",
          div(
            class = "section-header",
            paste0("Expression profiles: ", label1, " (left) vs ", label2, " (right)")
          ),
          uiOutput(ns("plots_ui"))
        )
      } else {
        NULL
      }

      tagList(overlap_box, phase_section, genes_section)
    })

    # ---- Expression plots UI (one row per gene) ----
    output$plots_ui <- renderUI({
      gl <- gene_list()
      cr <- comp_result()
      if (is.null(gl) || length(gl) == 0 || is.null(cr)) return(NULL)

      label1 <- get_dataset_label(cr$ds1)
      label2 <- get_dataset_label(cr$ds2)
      slot_count <- min(length(gl), plot_slot_count())
      visible_count <- min(length(gl), visible_plot_count())
      if (slot_count <= 0) return(NULL)
      visible_genes <- gl[seq_len(slot_count)]

      rows <- lapply(visible_genes, function(gene) {
        g_safe <- make.names(gene)
        tagList(
          div(
            class = "comp-expression-row",
            div(
              class = "comp-expression-plot",
              p(class = "comp-expression-label", label1),
              plotOutput(ns(paste0("plt1_", g_safe)), height = "260px", width = "260px")
            ),
            div(
              class = "comp-expression-plot",
              p(class = "comp-expression-label", label2),
              plotOutput(ns(paste0("plt2_", g_safe)), height = "260px", width = "260px")
            )
          )
        )
      })

      loading_status <- if (visible_count < slot_count) {
        div(
          class = "plot-loading-status",
          sprintf("Rendering expression rows %d of %d loaded...", visible_count, slot_count)
        )
      } else {
        NULL
      }

      footer <- div(
        class = "plot-pagination-footer",
        div(
          class = "plot-pagination-status",
          sprintf("Showing %d out of %d genes", slot_count, length(gl))
        ),
        if (slot_count < length(gl)) {
          actionButton(ns("load_more_plots"), "Load more", class = "btn btn-default btn-sm")
        } else {
          NULL
        }
      )

      tagList(loading_status, do.call(tagList, rows), footer)
    })

    # ---- Render individual gene plots ----
    observe({
      gl <- gene_list()
      cr <- comp_result()
      if (is.null(gl) || length(gl) == 0 || is.null(cr)) return()
      slot_count <- min(length(gl), plot_slot_count())
      if (slot_count <= 0) return()

      for (i in seq_len(slot_count)) {
        local({
          plot_index <- i
          g      <- gl[[i]]
          g_safe <- make.names(g)

          output[[paste0("plt1_", g_safe)]] <- renderCachedPlot({
            req(plot_index <= visible_plot_count())
            generate_comparison_plot(
              data_info = cr$plotdata1,
              ds_id     = cr$ds1,
              gene      = g,
              session   = session
            )
          }, cacheKeyExpr = {
            req(plot_index <= visible_plot_count())
            list(
              scope          = "comparison-expression",
              side           = "dataset1",
              dataset        = cr$ds1,
              gene           = g,
              show_fit       = FALSE,
              show_line      = TRUE,
              show_comparison = FALSE,
              data_signature = plot_data_cache_signature(cr$plotdata1)
            )
          }, cache = "session", width = 260, height = 260)

          output[[paste0("plt2_", g_safe)]] <- renderCachedPlot({
            req(plot_index <= visible_plot_count())
            generate_comparison_plot(
              data_info = cr$plotdata2,
              ds_id     = cr$ds2,
              gene      = g,
              session   = session
            )
          }, cacheKeyExpr = {
            req(plot_index <= visible_plot_count())
            list(
              scope          = "comparison-expression",
              side           = "dataset2",
              dataset        = cr$ds2,
              gene           = g,
              show_fit       = FALSE,
              show_line      = TRUE,
              show_comparison = FALSE,
              data_signature = plot_data_cache_signature(cr$plotdata2)
            )
          }, cache = "session", width = 260, height = 260)
        })
      }
    })

  }) # end moduleServer
}
