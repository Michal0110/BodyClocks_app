# comparison_functions.R
# Analysis helper functions for the Dataset Comparison module.
# No Shiny dependencies — pure R.

# ============================================================
# Hypergeometric test for gene overlap significance
# Ported from comparison_analysis.rmd
# ============================================================
hypergeometric_test <- function(list1, list2, genome_size) {
  overlap <- length(intersect(list1, list2))
  size1   <- length(list1)
  size2   <- length(list2)

  p_value <- phyper(overlap - 1, size1, genome_size - size1, size2,
                    lower.tail = FALSE)

  expected <- (size1 * size2) / genome_size
  fold_enrichment <- if (expected > 0) overlap / expected else NA_real_
  jaccard <- if (length(union(list1, list2)) > 0)
    overlap / length(union(list1, list2))
  else 0

  list(
    overlap         = overlap,
    size1           = size1,
    size2           = size2,
    expected        = expected,
    fold_enrichment = fold_enrichment,
    jaccard         = jaccard,
    p_value         = p_value
  )
}

# ============================================================
# Extract a normalised gene table (symbol, pval, bhq, phase) from display_dt.
# Handles single-condition and WT/KO display_dt formats.
# display_dt formats.
# ============================================================
extract_gene_data <- function(table_data, ds_id) {
  if (is.null(table_data)) return(NULL)

  # Drop rows with missing symbol (padding rows in some display_dt files)
  table_data <- table_data[!is.na(table_data$symbol) & nzchar(table_data$symbol), ]
  if (nrow(table_data) == 0) return(NULL)

  meta <- CIRC_DATASETS[[ds_id]]
  if (is.null(meta)) return(NULL)

  switch(meta$col_type,
    single = data.frame(
      symbol = table_data$symbol,
      pval   = if ("pVal" %in% names(table_data)) table_data$pVal
               else rep(NA_real_, nrow(table_data)),
      bhq    = table_data$BH.Q,
      phase  = if ("phase" %in% names(table_data)) table_data$phase
               else rep(NA_real_, nrow(table_data)),
      stringsAsFactors = FALSE
    ),
    wt_primary = data.frame(
      symbol = table_data$symbol,
      pval   = if ("pVal_wt" %in% names(table_data)) table_data$pVal_wt
               else rep(NA_real_, nrow(table_data)),
      bhq    = table_data$BH.Q_wt,
      phase  = if ("phase_wt" %in% names(table_data)) table_data$phase_wt
               else rep(NA_real_, nrow(table_data)),
      stringsAsFactors = FALSE
    ),
    ko_primary = data.frame(
      symbol = table_data$symbol,
      pval   = if ("pVal_ko" %in% names(table_data)) table_data$pVal_ko
               else rep(NA_real_, nrow(table_data)),
      bhq    = table_data$BH.Q_ko,
      phase  = if ("phase_ko" %in% names(table_data)) table_data$phase_ko
               else rep(NA_real_, nrow(table_data)),
      stringsAsFactors = FALSE
    ),
    NULL
  )
}

# ============================================================
# Generate a single expression plot for a specific dataset.
# Wraps fitted_plot_* without requiring input$dataset from a module.
# ============================================================
generate_comparison_plot <- function(data_info, ds_id, gene,
                                     show_fit        = FALSE,
                                     show_line       = TRUE,
                                     show_comparison = FALSE,
                                     session) {
  meta <- CIRC_DATASETS[[ds_id]]
  if (is.null(meta)) return(NULL)

  tryCatch({
    switch(meta$plot_fn,
      osmo    = fitted_plot_osmo(data_info$data, gene, data_info$period, data_info$t,
                                 session, show_fit, show_line),
      atlas   = fitted_plot_atlas(data_info$data, NULL, gene, data_info$period,
                                  session, show_fit, show_line),
      cart_wt = fitted_plot_cart(
                                 if (!is.null(data_info$data_wt)) data_info$data_wt else data_info$data,
                                 data_info$data_ko, gene,
                                 data_info$period, data_info$t, session,
                                 show_fit, show_line, show_comparison,
                                 is_primary_wt = TRUE),
      cart_ko = fitted_plot_cart(
                                 if (!is.null(data_info$data_ko)) data_info$data_ko else data_info$data,
                                 data_info$data_wt, gene,
                                 data_info$period, data_info$t, session,
                                 show_fit, show_line, show_comparison,
                                 is_primary_wt = FALSE),
      gobs    = fitted_plot_gobs(data_info$data, gene, data_info$period, data_info$t,
                                 session, show_fit, show_line),
      NULL
    )
  }, error = function(e) NULL)
}

# ============================================================
# Circular distance on a 24 h scale → result in [0, 12]
# ============================================================
circ_dist_24 <- function(a, b) {
  d <- abs(a - b) %% 24
  ifelse(d > 12, 24 - d, d)
}

# ============================================================
# Compute clock-relative phases (CRP) for overlapping genes.
# overlap_df must have columns phase1, phase2, symbol.
# Returns data.frame: symbol, phase1, phase2, crp1, crp2, dist, classification
# ============================================================
compute_crp_for_overlap <- function(overlap_df, bmal1_phase1, bmal1_phase2,
                                    threshold = 4) {
  df <- overlap_df[!is.na(overlap_df$phase1) & !is.na(overlap_df$phase2), ]
  if (nrow(df) == 0) return(NULL)

  # Subtract Bmal1 phase, map to [-12, 12)
  crp1 <- ((df$phase1 - bmal1_phase1) %% 24)
  crp1 <- ifelse(crp1 > 12, crp1 - 24, crp1)

  crp2 <- ((df$phase2 - bmal1_phase2) %% 24)
  crp2 <- ifelse(crp2 > 12, crp2 - 24, crp2)

  dist <- circ_dist_24(crp1, crp2)

  data.frame(
    symbol         = df$symbol,
    phase1         = df$phase1,
    phase2         = df$phase2,
    crp1           = crp1,
    crp2           = crp2,
    dist           = dist,
    classification = ifelse(dist <= threshold, "Clock-controlled", "Variable"),
    stringsAsFactors = FALSE
  )
}

# ============================================================
# Pairwise permutation test for CRP coherence.
# Null: shuffle crp2, recompute proportion within threshold.
# p-value = fraction of permutations >= observed proportion.
# ============================================================
run_crp_permutation_test <- function(crp1, crp2, threshold = 4,
                                     n_perm = 1000, seed = 42) {
  n_genes <- length(crp1)
  if (n_genes < 3) return(NULL)

  obs_n_cc <- sum(circ_dist_24(crp1, crp2) <= threshold, na.rm = TRUE)
  obs_pct  <- 100 * obs_n_cc / n_genes

  set.seed(seed)
  perm_pcts <- replicate(n_perm, {
    crp2_perm <- sample(crp2)
    100 * sum(circ_dist_24(crp1, crp2_perm) <= threshold, na.rm = TRUE) / n_genes
  })

  list(
    n_genes      = n_genes,
    n_coherent   = obs_n_cc,
    pct_coherent = obs_pct,
    pct_null     = mean(perm_pcts, na.rm = TRUE),
    p_value      = mean(perm_pcts >= obs_pct, na.rm = TRUE)
  )
}

# ============================================================
# Look up the human-readable label for a dataset ID via CIRC_DATASETS.
# ============================================================
get_dataset_label <- function(ds_id) {
  meta <- CIRC_DATASETS[[ds_id]]
  if (!is.null(meta) && nzchar(meta$dataset_label)) return(meta$dataset_label)
  ds_id  # fallback: return the id itself
}

# ============================================================
# Compute enrichment term overlap between two datasets.
# desc1/desc2:  cluster_descriptions data.frames (cluster, cluster_description, category)
# genes1/genes2: cluster_to_genes data.frames (cluster, genes [list column])
# category: category code (e.g. "Process", "KEGG", "RCTM", "WikiPathways")
# Returns data.frame: term, n_genes_ds1, n_genes_ds2, genes_ds1, genes_ds2
# ============================================================
compute_enrichment_overlap <- function(desc1, genes1, desc2, genes2, category,
                                       overlap_syms = character(0)) {
  if (is.null(desc1) || is.null(desc2)) return(data.frame())

  # Filter by category
  d1 <- desc1[desc1$category == category, , drop = FALSE]
  d2 <- desc2[desc2$category == category, , drop = FALSE]

  if (nrow(d1) == 0 || nrow(d2) == 0) return(data.frame())

  # Term identity is the category/STRING term ID key. Descriptions are display
  # labels only and are not guaranteed to be unique or identical over time.
  d1 <- d1[!duplicated(d1$cluster), , drop = FALSE]
  d2 <- d2[!duplicated(d2$cluster), , drop = FALSE]
  overlap_ids <- intersect(d1$cluster, d2$cluster)
  if (length(overlap_ids) == 0) return(data.frame())

  rows <- lapply(overlap_ids, function(cluster_id) {
    term_desc <- d1$cluster_description[match(cluster_id, d1$cluster)]

    # Look up gene lists — [row, "genes"] on a list-column returns a list, not a data.frame
    g1 <- genes1[
      genes1$category == category & genes1$cluster == cluster_id,
      "genes"
    ]
    g1_vec <- if (is.list(g1) && length(g1) > 0) g1[[1]] else character(0)

    g2 <- genes2[
      genes2$category == category & genes2$cluster == cluster_id,
      "genes"
    ]
    g2_vec <- if (is.list(g2) && length(g2) > 0) g2[[1]] else character(0)

    # Common rhythmic genes: term genes that are in the rhythmic overlap
    all_term_genes <- unique(c(g1_vec, g2_vec))
    common <- intersect(all_term_genes, overlap_syms)

    data.frame(
      term         = term_desc,
      n_genes_ds1  = length(g1_vec),
      n_genes_ds2  = length(g2_vec),
      n_common     = length(common),
      genes_ds1    = paste(g1_vec, collapse = ","),
      genes_ds2    = paste(g2_vec, collapse = ","),
      common_genes = paste(common, collapse = ","),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, rows)
}
