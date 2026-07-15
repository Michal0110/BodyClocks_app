# precompute_curves.R
#
# Pre-computes smooth fitted expression curves for all raw genes in the
# CIRC_DATASETS registry. Cache files are saved beside each raw expression
# source, e.g. cartilage_wt_data.rds -> cartilage_wt_curve_cache.rds.
#
# Usage:
#   source("precompute_curves.R")
#
# To force recompute existing caches:
#   FORCE_RECOMPUTE <- TRUE
#   source("precompute_curves.R")
#
# To process selected datasets or run a smoke test:
#   DATASET_IDS <- c("tendon")
#   MAX_GENES <- 3
#   FORCE_RECOMPUTE <- TRUE
#   source("precompute_curves.R")
#
# To control parallel workers:
#   N_WORKERS <- 4
#   source("precompute_curves.R")

source("datasets.R")
source("functions.R")

force_recompute <- exists("FORCE_RECOMPUTE") && isTRUE(FORCE_RECOMPUTE)
selected_dataset_ids <- if (exists("DATASET_IDS")) as.character(DATASET_IDS) else names(CIRC_DATASETS)
max_genes <- if (exists("MAX_GENES")) as.integer(MAX_GENES) else NA_integer_
n_workers <- if (exists("N_WORKERS")) as.integer(N_WORKERS) else 4L

if (!is.na(max_genes) && max_genes <= 0) {
  stop("MAX_GENES must be a positive integer when supplied.")
}

if (is.na(n_workers) || n_workers <= 0) {
  stop("N_WORKERS must be a positive integer when supplied.")
}

n_workers <- max(1L, n_workers)

unknown_dataset_ids <- setdiff(selected_dataset_ids, names(CIRC_DATASETS))
if (length(unknown_dataset_ids) > 0) {
  stop("Unknown DATASET_IDS: ", paste(unknown_dataset_ids, collapse = ", "))
}

dataset_source_entries <- function(ds_id, meta) {
  if (identical(meta$col_type, "single")) {
    list(list(
      dataset_ids = ds_id,
      label       = ds_id,
      source_file = file.path(meta$folder, paste0(meta$prefix, "_data.rds")),
      period      = meta$period,
      t           = meta$t
    ))
  } else {
    list(
      list(
        dataset_ids = ds_id,
        label       = paste0(ds_id, ":wt"),
        source_file = meta$wt_data_file,
        period      = meta$period,
        t           = meta$t
      ),
      list(
        dataset_ids = ds_id,
        label       = paste0(ds_id, ":ko"),
        source_file = meta$ko_data_file,
        period      = meta$period,
        t           = meta$t
      )
    )
  }
}

dedupe_curve_sources <- function(entries) {
  result <- list()

  for (entry in entries) {
    key <- normalize_curve_cache_path(entry$source_file)
    existing <- result[[key]]

    if (is.null(existing)) {
      result[[key]] <- entry
      next
    }

    same_period <- identical(as.numeric(existing$period), as.numeric(entry$period))
    same_t <- same_curve_t(existing$t, entry$t)
    if (!same_period || !same_t) {
      stop(
        "Conflicting period/t metadata for source file ",
        entry$source_file,
        ". Cache filenames are based on source file names, so this must be resolved in datasets.R."
      )
    }

    existing$dataset_ids <- unique(c(existing$dataset_ids, entry$dataset_ids))
    existing$label <- paste(unique(c(existing$label, entry$label)), collapse = ", ")
    result[[key]] <- existing
  }

  unname(result)
}

all_entries <- unlist(
  Map(dataset_source_entries, selected_dataset_ids, CIRC_DATASETS[selected_dataset_ids]),
  recursive = FALSE
)
queue <- dedupe_curve_sources(all_entries)

message(sprintf("Discovered %d unique expression source(s).", length(queue)))

precompute_one_source <- function(entry, source_index, n_sources) {
  source_file <- entry$source_file
  cache_file <- curve_cache_path(source_file)

  if (!file.exists(source_file)) {
    message(sprintf(
      "[%d/%d] Skipping missing source: %s",
      source_index, n_sources, source_file
    ))
    return(list(source_file = source_file, cache_file = cache_file, status = "missing"))
  }

  message(sprintf(
    "[%d/%d] Loading %s",
    source_index, n_sources, source_file
  ))

  source_data <- readRDS(source_file)
  gene_values <- as.character(source_data$Gene)
  valid_gene <- !is.na(gene_values) & nzchar(gene_values)
  rows_by_gene <- split(which(valid_gene), gene_values[valid_gene])
  genes <- names(rows_by_gene)

  if (!is.na(max_genes)) {
    genes <- head(genes, max_genes)
    rows_by_gene <- rows_by_gene[genes]
  }

  expected_n_genes <- length(genes)
  current_cache <- load_curve_cache(source_file, entry$period, entry$t, expected_n_genes)
  if (!force_recompute && !is.null(current_cache)) {
    message(sprintf(
      "  Skipping current cache: %s (%d gene(s))",
      cache_file, expected_n_genes
    ))
    return(list(source_file = source_file, cache_file = cache_file, status = "skipped"))
  }

  message(sprintf(
    "  Computing %d gene curve(s) for %s",
    expected_n_genes, paste(entry$dataset_ids, collapse = ", ")
  ))

  curves <- vector("list", expected_n_genes)
  names(curves) <- genes
  failures <- character(0)

  for (i in seq_along(genes)) {
    gene <- genes[[i]]
    filtered_data <- source_data[rows_by_gene[[gene]], , drop = FALSE]
    fit_error <- NULL

    result <- tryCatch(
      suppressWarnings(fit_filtered_gene_curve(filtered_data, gene, entry$period, entry$t, quiet = TRUE)),
      error = function(e) {
        fit_error <<- e$message
        NULL
      }
    )

    if (is.null(result) || is.null(result$curve_data)) {
      failures <- c(failures, if (is.null(fit_error)) gene else paste0(gene, ": ", fit_error))
    } else {
      curves[[gene]] <- list(curve_data = result$curve_data)
    }

    if (i %% 500 == 0 || i == expected_n_genes) {
      message(sprintf(
        "    %d/%d gene(s) complete; %d failed",
        i, expected_n_genes, length(failures)
      ))
    }
  }

  curves <- Filter(Negate(is.null), curves)
  sig <- curve_source_signature(source_file)

  cache <- list(
    version      = CURVE_CACHE_VERSION,
    source_file  = sig$source_file,
    source_size  = sig$source_size,
    source_mtime = sig$source_mtime,
    period       = as.numeric(entry$period),
    t            = if (is.null(entry$t)) NULL else as.numeric(entry$t),
    generated_at = Sys.time(),
    n_genes      = expected_n_genes,
    failures     = failures,
    curves       = curves
  )

  saveRDS(cache, cache_file)
  message(sprintf(
    "  Saved: %s (%d cached, %d failed)",
    cache_file, length(curves), length(failures)
  ))

  list(source_file = source_file, cache_file = cache_file, status = "computed")
}

run_precompute_queue <- function(queue, n_workers) {
  if (length(queue) == 0) {
    return(list())
  }

  n_workers <- min(n_workers, length(queue))
  if (n_workers <= 1L) {
    message("Running with 1 worker.")
    return(lapply(seq_along(queue), function(i) {
      precompute_one_source(queue[[i]], i, length(queue))
    }))
  }

  if (.Platform$OS.type != "windows") {
    message(sprintf("Running with %d forked workers.", n_workers))
    return(parallel::mclapply(seq_along(queue), function(i) {
      precompute_one_source(queue[[i]], i, length(queue))
    }, mc.cores = n_workers))
  }

  message(sprintf("Running with %d parallel workers.", n_workers))
  cl <- tryCatch(
    parallel::makeCluster(n_workers, outfile = "", setup_timeout = 15),
    error = function(e) {
      message("Parallel worker setup failed; falling back to 1 worker. Reason: ", e$message)
      NULL
    }
  )

  if (is.null(cl)) {
    return(lapply(seq_along(queue), function(i) {
      precompute_one_source(queue[[i]], i, length(queue))
    }))
  }

  on.exit(parallel::stopCluster(cl), add = TRUE)

  parallel::clusterExport(
    cl,
    varlist = c(
      "queue",
      "force_recompute",
      "max_genes",
      "precompute_one_source",
      "CURVE_CACHE_VERSION",
      "curve_cache_path",
      "curve_source_signature",
      "fit_filtered_gene_curve",
      "is_curve_cache_current",
      "load_curve_cache",
      "normalize_curve_cache_path",
      "same_curve_t"
    ),
    envir = environment()
  )

  parallel::parLapplyLB(cl, seq_along(queue), function(i) {
    precompute_one_source(queue[[i]], i, length(queue))
  })
}

results <- run_precompute_queue(queue, n_workers)

status_counts <- table(vapply(results, `[[`, character(1), "status"))
message("Curve precompute complete.")
message(paste(capture.output(print(status_counts)), collapse = "\n"))

invisible(results)
