# precompute_positions.R
#
# Pre-computes visNetwork node positions using the original layout_nicely +
# vis.js physics simulation. Runs a temporary Shiny app that renders each
# network, lets physics stabilise, captures final positions, and saves them.
#
# Usage:
#   source("precompute_positions.R")
#   # Browser opens. Networks compute one by one. Positions saved. Window closes.
#
# To force recompute existing positions:
#   FORCE_RECOMPUTE <- TRUE
#   source("precompute_positions.R")

library(shiny)
library(visNetwork)
library(jsonlite)

force_recompute <- exists("FORCE_RECOMPUTE") && isTRUE(FORCE_RECOMPUTE)

# STRING score thresholds to precompute positions for
STRING_THRESHOLDS <- c(700, 900)

# ---- Dataset registry --------------------------------------------------------

# Derive circadian datasets from the central registry in datasets.R.
# This ensures precompute_positions.R always covers the same set as the app.
source("datasets.R")
circ_datasets <- lapply(CIRC_DATASETS, function(m) list(folder = m$folder, prefix = m$prefix))

# Treatment datasets (unprefixed files: nodes_with_tooltips.rds, etc.)
discover_treatment_datasets <- function() {
  root <- "data/comparisons"
  if (!dir.exists(root)) return(list())
  nf <- list.files(root, pattern = "^nodes_with_tooltips\\.rds$",
                   recursive = TRUE, full.names = TRUE)
  lapply(nf, function(f) list(folder = dirname(f), prefix = NULL))
}

treat_datasets <- discover_treatment_datasets()

# Helper: file paths for a dataset entry.
# For circadian entries ds$threshold selects threshold-specific files;
# for treatment entries (prefix = NULL) or entries without a threshold, the
# unsuffixed legacy files are used.
ds_files <- function(ds) {
  threshold <- ds$threshold
  if (is.null(ds$prefix)) {
    # Treatment dataset — no threshold variant
    list(
      nodes = file.path(ds$folder, "nodes_with_tooltips.rds"),
      edges = file.path(ds$folder, "edges.rds"),
      pos   = file.path(ds$folder, "positions.rds"),
      label = basename(ds$folder)
    )
  } else if (!is.null(threshold)) {
    # Circadian dataset with a specific STRING score threshold
    edge_file <- file.path(ds$folder, paste0(ds$prefix, "_edges_string", threshold, ".rds"))
    # Fall back to the unsuffixed edges file when the threshold-700 file is absent
    if (!file.exists(edge_file) && threshold == 700) {
      edge_file <- file.path(ds$folder, paste0(ds$prefix, "_edges.rds"))
    }
    list(
      nodes = file.path(ds$folder, paste0(ds$prefix, "_nodes_with_tooltips.rds")),
      edges = edge_file,
      pos   = file.path(ds$folder, paste0(ds$prefix, "_positions_string", threshold, ".rds")),
      label = paste0(ds$prefix, "_string", threshold)
    )
  } else {
    # Circadian dataset without threshold (legacy unsuffixed files)
    list(
      nodes = file.path(ds$folder, paste0(ds$prefix, "_nodes_with_tooltips.rds")),
      edges = file.path(ds$folder, paste0(ds$prefix, "_edges.rds")),
      pos   = file.path(ds$folder, paste0(ds$prefix, "_positions.rds")),
      label = ds$prefix
    )
  }
}

# Expand circadian datasets across all STRING thresholds; treatment datasets are
# threshold-agnostic (prefix = NULL, threshold = NULL).
circ_entries  <- unlist(lapply(circ_datasets, function(ds) {
  lapply(STRING_THRESHOLDS, function(thr) {
    list(folder = ds$folder, prefix = ds$prefix, threshold = thr)
  })
}), recursive = FALSE)

treat_entries <- lapply(treat_datasets, function(ds) {
  list(folder = ds$folder, prefix = NULL, threshold = NULL)
})

all_entries <- c(circ_entries, treat_entries)

# Build queue — skip entries whose edge file is missing or positions already exist
queue <- Filter(function(ds) {
  f <- ds_files(ds)
  if (!file.exists(f$nodes)) { message("Skipping (missing nodes): ", f$nodes); return(FALSE) }
  if (!file.exists(f$edges)) { message("Skipping (missing edges): ", f$edges); return(FALSE) }
  if (!force_recompute && file.exists(f$pos)) { message("Skipping (exists): ", f$pos); return(FALSE) }
  TRUE
}, all_entries)

# ---- Run ---------------------------------------------------------------------

if (length(queue) == 0) {
  message("All positions already computed. Set FORCE_RECOMPUTE <- TRUE to recompute.")
} else {
  message(sprintf("Computing positions for %d dataset(s)...\n", length(queue)))

  ui <- fluidPage(
    h3("Precomputing network positions"),
    textOutput("status"),
    visNetworkOutput("network", height = "700px", width = "100%")
  )

  server <- function(input, output, session) {
    idx <- reactiveVal(1)

    output$status <- renderText({
      i <- idx()
      if (i > length(queue)) return("All done! Closing...")
      f <- ds_files(queue[[i]])
      paste0("Processing ", i, "/", length(queue), ": ", f$label)
    })

    output$network <- renderVisNetwork({
      i <- idx()
      req(i <= length(queue))
      ds <- queue[[i]]
      f  <- ds_files(ds)

      nodes <- readRDS(f$nodes)
      edges <- readRDS(f$edges)
      edges <- edges[edges$from %in% nodes$id & edges$to %in% nodes$id, ]

      message(sprintf("[%d/%d] %s: %d nodes, %d edges",
                      i, length(queue), f$label, nrow(nodes), nrow(edges)))

      visNetwork(nodes, edges) %>%
        visIgraphLayout(layout = "layout_nicely", physics = TRUE,
                        smooth = TRUE, randomSeed = 69420) %>%
        visPhysics(stabilization = list(iterations = 200)) %>%
        visEvents(stabilizationIterationsDone = paste0(
          "function() {",
          "  this.setOptions({physics: false});",
          "  var positions = this.getPositions();",
          "  Shiny.setInputValue('node_positions',",
          "    JSON.stringify(positions), {priority: 'event'});",
          "}"
        ))
    })

    observeEvent(input$node_positions, {
      i  <- idx()
      ds <- queue[[i]]
      f  <- ds_files(ds)

      pos_raw <- jsonlite::fromJSON(input$node_positions)
      pos_df  <- data.frame(
        id = names(pos_raw),
        x  = vapply(pos_raw, function(p) p$x, numeric(1)),
        y  = vapply(pos_raw, function(p) p$y, numeric(1)),
        stringsAsFactors = FALSE
      )

      saveRDS(pos_df, f$pos)
      message("  Saved: ", f$pos, " (", nrow(pos_df), " nodes)")

      next_i <- i + 1
      if (next_i > length(queue)) {
        message("\nAll positions computed!")
        Sys.sleep(1)
        stopApp()
      } else {
        idx(next_i)
      }
    })
  }

  runApp(shinyApp(ui, server), launch.browser = TRUE)
}
