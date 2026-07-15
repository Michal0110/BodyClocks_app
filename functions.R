# ---- Curve cache helpers ----------------------------------------------------

CURVE_CACHE_VERSION <- 1L
LIVE_CURVE_CACHE <- new.env(parent = emptyenv())

normalize_curve_cache_path <- function(path) {
  sub("^\\./", "", gsub("\\\\", "/", path))
}

curve_cache_path <- function(source_file) {
  dir <- dirname(source_file)
  stem <- sub("\\.[Rr][Dd][Ss]$", "", basename(source_file))

  if (grepl("_data$", stem)) {
    stem <- sub("_data$", "_curve_cache", stem)
  } else {
    stem <- paste0(stem, "_curve_cache")
  }

  file.path(dir, paste0(stem, ".rds"))
}

curve_source_signature <- function(source_file) {
  info <- file.info(source_file)
  if (is.na(info$size)) return(NULL)

  list(
    source_file  = normalize_curve_cache_path(source_file),
    source_size  = unname(info$size),
    source_mtime = unname(as.numeric(info$mtime))
  )
}

same_curve_t <- function(x, y) {
  if (is.null(x) && is.null(y)) return(TRUE)
  if (is.null(x) || is.null(y)) return(FALSE)
  identical(as.numeric(x), as.numeric(y))
}

is_curve_cache_current <- function(cache, source_file, period, t,
                                   expected_n_genes = NULL) {
  sig <- curve_source_signature(source_file)
  if (is.null(sig) || is.null(cache) || !is.list(cache)) return(FALSE)

  isTRUE(cache$version == CURVE_CACHE_VERSION) &&
    identical(cache$source_file, sig$source_file) &&
    identical(as.numeric(cache$source_size), as.numeric(sig$source_size)) &&
    identical(as.numeric(cache$source_mtime), as.numeric(sig$source_mtime)) &&
    identical(as.numeric(cache$period), as.numeric(period)) &&
    same_curve_t(cache$t, t) &&
    (is.null(expected_n_genes) || identical(as.integer(cache$n_genes), as.integer(expected_n_genes))) &&
    is.list(cache$curves)
}

load_curve_cache <- function(source_file, period, t, expected_n_genes = NULL) {
  cache_file <- curve_cache_path(source_file)
  if (!file.exists(cache_file)) return(NULL)

  cache <- tryCatch(readRDS(cache_file), error = function(e) NULL)
  if (!is_curve_cache_current(cache, source_file, period, t, expected_n_genes)) {
    return(NULL)
  }

  cache
}

attach_curve_cache <- function(data, source_file, period, t, verbose = FALSE) {
  attr(data, "curve_source_signature") <- curve_source_signature(source_file)
  cache <- load_curve_cache(source_file, period, t)

  if (!is.null(cache)) {
    attr(data, "curve_cache") <- cache
    attr(data, "curve_cache_file") <- curve_cache_path(source_file)
    if (isTRUE(verbose)) {
      message("Loaded curve cache: ", curve_cache_path(source_file))
    }
  } else if (isTRUE(verbose)) {
    message("No current curve cache: ", curve_cache_path(source_file))
  }

  data
}

curve_data_is_valid <- function(curve_data) {
  is.data.frame(curve_data) && all(c("t", "y") %in% names(curve_data))
}

get_cached_curve_data <- function(data, gene) {
  cache <- attr(data, "curve_cache", exact = TRUE)
  if (is.null(cache) || !is.list(cache$curves)) return(NULL)

  entry <- cache$curves[[gene]]
  if (is.null(entry) || is.null(entry$curve_data)) return(NULL)

  curve_data <- entry$curve_data
  if (!curve_data_is_valid(curve_data)) return(NULL)

  curve_data
}

live_curve_cache_key <- function(data, gene, period, t) {
  sig <- attr(data, "curve_source_signature", exact = TRUE)
  if (is.null(sig)) return(NULL)

  t_key <- if (is.null(t)) "NULL" else paste(as.numeric(t), collapse = ",")
  paste(
    sig$source_file,
    sig$source_size,
    sig$source_mtime,
    as.numeric(period),
    t_key,
    gene,
    sep = "|"
  )
}

get_live_curve_data <- function(data, gene, period, t) {
  key <- live_curve_cache_key(data, gene, period, t)
  if (is.null(key) || !exists(key, envir = LIVE_CURVE_CACHE, inherits = FALSE)) {
    return(NULL)
  }

  curve_data <- get(key, envir = LIVE_CURVE_CACHE, inherits = FALSE)
  if (!curve_data_is_valid(curve_data)) return(NULL)

  curve_data
}

set_live_curve_data <- function(data, gene, period, t, curve_data) {
  if (!curve_data_is_valid(curve_data)) return(invisible(FALSE))

  key <- live_curve_cache_key(data, gene, period, t)
  if (is.null(key)) return(invisible(FALSE))

  assign(key, curve_data, envir = LIVE_CURVE_CACHE)
  invisible(TRUE)
}

plot_data_cache_signature <- function(data_info) {
  data_signature <- function(data) {
    sig <- attr(data, "curve_source_signature", exact = TRUE)
    if (is.null(sig)) return(NULL)

    list(
      source_file  = sig$source_file,
      source_size  = sig$source_size,
      source_mtime = sig$source_mtime
    )
  }

  list(
    period  = data_info$period,
    t       = data_info$t,
    data    = data_signature(data_info$data),
    data_wt = data_signature(data_info$data_wt),
    data_ko = data_signature(data_info$data_ko)
  )
}

# Function for fitting and generating curve from an already-filtered gene table.
fit_filtered_gene_curve <- function(filtered_data, gene, period, t, quiet = FALSE) {
  log_fit_message <- function(...) {
    if (!isTRUE(quiet)) cat(...)
  }

  y <- filtered_data$mean
  x <- as.numeric(filtered_data$time_point)  # Ensure x is numeric
  
  # Define the sine and cosine functions
  sine_func <- function(x, a, b, c, d, period) {
    a * sin(b * (x - c) * (2 * pi / period)) + d
  }
  
  cosine_func <- function(x, a, b, c, d, period) {
    a * cos(b * (x - c) * (2 * pi / period)) + d
  }
  
  # Initialize fits
  fit_sine <- NULL
  fit_cosine <- NULL
  fit_sine2 <- NULL
  fit_cosine2 <- NULL
  
  # Try fitting the first sine function
  fit_sine <- tryCatch({
    nls(y ~ sine_func(x, a, b, c, d, period),
        start = list(a = 1000, b = 2 * pi / period, c = 0, d = mean(y)),
        lower = c(a = -Inf, b = 0, c = -Inf, d = -Inf),
        upper = c(a = Inf, b = 2 * pi / (period / 2), c = Inf, d = Inf),
        control = nls.control(maxiter = 1000, minFactor = 1e-7))
  }, error = function(e) {
    log_fit_message("Fitting sine failed:", e$message, "\n")
    NULL
  })
  
  # Try fitting the first cosine function
  fit_cosine <- tryCatch({
    nls(y ~ cosine_func(x, a, b, c, d, period),
        start = list(a = 1000, b = 2 * pi / period, c = 0, d = mean(y)),
        lower = c(a = -Inf, b = 0, c = -Inf, d = -Inf),
        upper = c(a = Inf, b = 2 * pi / (period / 2), c = Inf, d = Inf),
        control = nls.control(maxiter = 1000, minFactor = 1e-7))
  }, error = function(e) {
    log_fit_message("Fitting cosine failed:", e$message, "\n")
    NULL
  })
  
  # Try fitting the second sine function
  fit_sine2 <- tryCatch({
    nls(y ~ sine_func(x, a, b, c, d, period),
        start = list(a = 1000, b = 1, c = 0, d = mean(y)),
        lower = c(a = -Inf, b = 0, c = -Inf, d = -Inf),
        upper = c(a = Inf, b = Inf, c = Inf, d = Inf),
        control = nls.control(maxiter = 1000, minFactor = 1e-7))
  }, error = function(e) {
    log_fit_message("Fitting sine 2 failed:", e$message, "\n")
    NULL
  })
  
  # Try fitting the second cosine function
  fit_cosine2 <- tryCatch({
    nls(y ~ cosine_func(x, a, b, c, d, period),
        start = list(a = 1000000, b = 1, c = 0, d = mean(y)),
        lower = c(a = -Inf, b = 0, c = -Inf, d = -Inf),
        upper = c(a = Inf, b = Inf, c = Inf, d = Inf),
        control = nls.control(maxiter = 1000, minFactor = 1e-7))
  }, error = function(e) {
    log_fit_message("Fitting cosine 2 failed:", e$message, "\n")
    NULL
  })
  
  # Choose the best fit
  fits <- list(fit_sine, fit_cosine, fit_sine2, fit_cosine2)
  fits <- Filter(Negate(is.null), fits)  # Remove NULL elements
  
  if (length(fits) == 0) {
    log_fit_message("All fitting attempts failed for", gene, "\n")
    # Fallback to a simple cosine model
    simple_cosine_func <- function(x, a, b, c, d) {
      a * cos(b * (x - c)) + d
    }
    fit <- tryCatch({
      nls(y ~ simple_cosine_func(x, a, b, c, d),
          start = list(a = 1000, b = 2 * pi / period, c = 0, d = mean(y)),
          control = nls.control(maxiter = 1000, minFactor = 1e-7))
    }, error = function(e) {
      log_fit_message("Fallback fitting cosine failed:", e$message, "\n")
      NULL
    })
    if (is.null(fit)) {
      # Fallback to a polynomial regression if cosine fitting fails
      fit <- tryCatch({
        lm(y ~ poly(x, 2))  # Quadratic polynomial
      }, error = function(e) {
        log_fit_message("Fallback fitting polynomial failed:", e$message, "\n")
        return(NULL)
      })
      if (is.null(fit)) {
        return(NULL)
      }
    }
    t_smooth <- seq(min(x), max(x), length.out = 100)
    y_smooth <- predict(fit, newdata = data.frame(x = t_smooth))
    curve_data <- data.frame(t = t_smooth, y = y_smooth)
    return(list(filtered_data = filtered_data, curve_data = curve_data, fit = fit))
  }
  
  deviances <- sapply(fits, function(fit) deviance(fit))
  
  # Calculate frequencies for each fit
  frequencies <- sapply(fits, function(fit) {
    if (inherits(fit, "nls")) {
      coef(fit)["b"] / (2 * pi / period)
    } else {
      NA
    }
  })
  
  # Calculate the difference between frequencies and target frequency (one cycle in the time series)
  frequency_diffs <- abs(frequencies - 1 / period)
  
  # Calculate a score based on deviance and frequency difference
  scores <- deviances * frequency_diffs
  
  # Select the fit with the lowest score
  selected_index <- which.min(scores)
  fit <- fits[[selected_index]]
  
  # Generate values for the smooth curve
  t_smooth <- seq(min(x), max(x), length.out = 100)
  y_smooth <- predict(fit, newdata = data.frame(x = t_smooth))
  
  curve_data <- data.frame(t = t_smooth, y = y_smooth)
  
  return(list(filtered_data = filtered_data, curve_data = curve_data, fit = fit))
}

# Function for filtering the raw expression points used by expression plots.
filter_gene_expression_data <- function(data, gene, quiet = FALSE) {
  if (is.null(data) || !"Gene" %in% names(data)) {
    if (!isTRUE(quiet)) cat("No expression data available for gene:", gene, "\n")
    return(NULL)
  }

  gene_match <- !is.na(data$Gene) & data$Gene == gene
  filtered_data <- data[gene_match, , drop = FALSE]

  # Check if the subset is empty
  if (nrow(filtered_data) == 0) {
    if (!isTRUE(quiet)) cat("No data found for gene:", gene, "\n")
    return(NULL)
  }

  filtered_data
}

# Function for preparing expression plot data. Curves are only fetched/fitted
# when include_fit is TRUE.
get_gene_plot_data <- function(data, gene, period, t, include_fit = FALSE, quiet = FALSE) {
  filtered_data <- filter_gene_expression_data(data, gene, quiet = quiet)
  if (is.null(filtered_data)) return(NULL)

  if (!isTRUE(include_fit)) {
    return(list(filtered_data = filtered_data, curve_data = NULL, fit = NULL))
  }

  curve_data <- get_cached_curve_data(data, gene)
  if (!is.null(curve_data)) {
    return(list(filtered_data = filtered_data, curve_data = curve_data, fit = NULL))
  }

  curve_data <- get_live_curve_data(data, gene, period, t)
  if (!is.null(curve_data)) {
    return(list(filtered_data = filtered_data, curve_data = curve_data, fit = NULL))
  }

  result <- fit_filtered_gene_curve(filtered_data, gene, period, t, quiet = quiet)
  if (!is.null(result) && !is.null(result$curve_data)) {
    set_live_curve_data(data, gene, period, t, result$curve_data)
  }

  result
}

# Function for fitting and generating curve
fit_and_generate_curve <- function(data, gene, period, t, quiet = FALSE) {
  get_gene_plot_data(data, gene, period, t, include_fit = TRUE, quiet = quiet)
}



# Updated fitted_plot_cart function with dynamic dataset labels
# Updated fitted_plot_cart function with dynamic dataset labels and proper sizing
fitted_plot_cart <- function(data1, data2, gene, period, t, session, show_fit = FALSE, show_line = TRUE, show_comparison = TRUE, is_primary_wt = TRUE) {
  result1 <- get_gene_plot_data(data1, gene, period, t, include_fit = show_fit, quiet = TRUE)
  result2 <- NULL
  
  if (show_comparison && !is.null(data2)) {
    result2 <- get_gene_plot_data(data2, gene, period, t, include_fit = show_fit, quiet = TRUE)
  }
  
  # Skip genes that are not available in the raw expression data.
  if (is.null(result1)) {
    return(NULL)
  }
  
  # Extract data for plotting
  filtered_data1 <- result1$filtered_data
  curve_data1 <- result1$curve_data
  
  # Extract BH.Q values
  bhq1 <- format(filtered_data1$BH.Q[1], scientific = FALSE, digits = 2)
  
  # Add the sizing code here
  # Adjust theme based on plot dimensions
  plot_width <- session$clientData$output_genePlot_width
  plot_height <- session$clientData$output_genePlot_height
  
  # Ensure plot dimensions are valid
  if (is.null(plot_width) || is.null(plot_height) || plot_width <= 0 || plot_height <= 0) {
    plot_width <- 500
    plot_height <- 375  # Maintain 3:4 aspect ratio
  } else {
    plot_height <- plot_width * 3 / 4  # Maintain 3:4 aspect ratio
  }
  
  # Calculate scaling factor
  scaling_factor <- min(plot_width, plot_height) / 200
  
  # Adjust sizes based on scaling factor
  text_size <- scaling_factor * 6
  line_size <- scaling_factor * 0.75
  point_size <- scaling_factor * 1.5
  annotate_text_size <- scaling_factor * 2
  
  # Set colors and labels based on which dataset is primary
  primary_color <- ifelse(is_primary_wt, "blue", "orange")
  comparison_color <- ifelse(is_primary_wt, "orange", "blue")
  primary_label <- ifelse(is_primary_wt, "WT", "KO")
  comparison_label <- ifelse(is_primary_wt, "KO", "WT")
  
  # Plot with primary dataset
  p <- ggplot() +
    theme_classic() +
    theme(
      plot.title = element_text(size = text_size * 1.5, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(size = text_size, face = "bold"),
      axis.text.y = element_text(size = text_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, title = gene) +
    geom_rect(
      aes(xmin = 3.5, xmax = 6.5,
          ymin = -Inf, ymax = Inf),
      fill = "grey90"
    ) +
    geom_rect(
      aes(xmin = 9.5, xmax = 12.5,
          ymin = -Inf, ymax = Inf),
      fill = "grey90"
    ) +
    geom_point(data = filtered_data1, aes(x = as.numeric(time_point), y = mean), 
               color = primary_color, size = point_size) +
    scale_x_continuous(breaks = as.numeric(filtered_data1$time_point),
                       labels = {
                         n   <- length(filtered_data1$time_point)
                         lbl <- rep("", n)
                         lbl[seq(1,        n, by = period)] <- "ZT2"
                         lbl[seq(1 + period/2, n, by = period)] <- "ZT14"
                         lbl
                       })
  
  # Add comparison dataset if enabled and available
  if (show_comparison && !is.null(result2)) {
    filtered_data2 <- result2$filtered_data
    curve_data2 <- result2$curve_data
    bhq2 <- format(filtered_data2$BH.Q[1], scientific = FALSE, digits = 2)
    
    p <- p + 
      geom_point(data = filtered_data2, aes(x = as.numeric(time_point), y = mean), 
                color = comparison_color, size = point_size) +
      annotate("text", x = Inf, y = Inf, 
              label = paste(primary_label, "bhq:", bhq1), 
              color = primary_color, hjust = 1.1, vjust = 2, 
              size = annotate_text_size, fontface = "bold") +
      annotate("text", x = Inf, y = Inf, 
              label = paste(comparison_label, "bhq:", bhq2), 
              color = comparison_color, hjust = 1.1, vjust = 3.5, 
              size = annotate_text_size, fontface = "bold")
  } else {
    # Just show the primary dataset label
    p <- p + annotate("text", x = Inf, y = Inf, 
                     label = paste(primary_label, "bhq:", bhq1), 
                     color = primary_color, hjust = 1.1, vjust = 2, 
                     size = annotate_text_size, fontface = "bold")
  }
  
  # Conditionally add the sine fitting line
  if (show_fit) {
    p <- p + geom_line(data = curve_data1, aes(x = t, y = y), 
                      color = primary_color, size = line_size)
    if (show_comparison && !is.null(result2)) {
      p <- p + geom_line(data = curve_data2, aes(x = t, y = y), 
                        color = comparison_color, size = line_size, 
                        linetype = "dashed")
    }
  }
  
  # Conditionally add the simple line connecting points
  if (show_line) {
    p <- p + geom_line(data = filtered_data1, aes(x = as.numeric(time_point), y = mean), 
                      color = primary_color, size = line_size)
    if (show_comparison && !is.null(result2)) {
      p <- p + geom_line(data = filtered_data2, aes(x = as.numeric(time_point), y = mean), 
                        color = comparison_color, size = line_size)
    }
  }
  
  return(p)
}



# Define the plotting function for gobs
fitted_plot_gobs <- function(data, gene, period, t, session, show_fit = FALSE, show_line = TRUE) {
  
  result <- get_gene_plot_data(data, gene, period, t, include_fit = show_fit, quiet = TRUE)
  # Skip genes that are not available in the raw expression data.
  if (is.null(result)) {
    return(NULL)
  }
  
  # Extract data for plotting
  filtered_data <- result$filtered_data
  curve_data <- result$curve_data
  
  # Extract time points for x-axis
  x <- filtered_data$time_point
  
  # Extract BH.Q values
  bhq <- format(filtered_data$BH.Q[1], scientific = FALSE, digits = 2)
  
  # Adjust theme based on plot dimensions
  plot_width <- session$clientData$output_genePlot_width
  plot_height <- session$clientData$output_genePlot_height
  
  # Ensure plot dimensions are valid
  if (is.null(plot_width) || is.null(plot_height) || plot_width <= 0 || plot_height <= 0) {
    plot_width <- 500
    plot_height <- 375  # Maintain 3:4 aspect ratio
  } else {
    plot_height <- plot_width * 3 / 4  # Maintain 3:4 aspect ratio
  }
  
  # Calculate scaling factor
  scaling_factor <- min(plot_width, plot_height) / 200
  
  # Adjust sizes based on scaling factor
  text_size <- scaling_factor * 6
  line_size <- scaling_factor * 0.75
  point_size <- scaling_factor * 1.5
  annotate_text_size <- scaling_factor * 2.5
  
  p <- ggplot() +
    theme_classic() +
    theme(
      plot.title = element_text(size = text_size * 1.5, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(size = text_size, face = "bold"),
      axis.text.y = element_text(size = text_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, title = gene) +
    geom_rect(
      aes(xmin = 3.5, xmax = 6.5,
          ymin = -Inf, ymax = Inf),
      fill = "grey90"
    ) +
    geom_point(data = filtered_data, aes(x = as.numeric(time_point), y = mean), color = "black", size = point_size) +
    geom_errorbar(data = filtered_data, aes(x = as.numeric(time_point), ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.9)) +
    scale_x_continuous(breaks = as.numeric(x), labels = c("ZT2", "6", "10", "14", "18", "22")) +
    annotate("text", x = Inf, y = Inf, label = paste("bhq:", bhq), color = "black", hjust = 1.1, vjust = 2, size = annotate_text_size, fontface = "bold")
  
  # Conditionally add the sine fitting line
  if (show_fit) {
    p <- p + geom_line(data = curve_data, aes(x = t, y = y), color = "black", size = line_size)
  }
  
  # Conditionally add the simple line connecting points
  if (show_line) {
    p <- p + geom_line(data = filtered_data, aes(x = as.numeric(time_point), y = mean), color = "blue", size = line_size)
  }
  return(p)
}



# Define the plotting function for osmo
fitted_plot_osmo <- function(data, gene, period, t, session, show_fit = FALSE, show_line = TRUE) {
  
  result <- get_gene_plot_data(data, gene, period, t, include_fit = show_fit, quiet = TRUE)
  # Skip genes that are not available in the raw expression data.
  if (is.null(result)) {
    return(NULL)
  }
  
  # Extract data for plotting
  filtered_data <- result$filtered_data
  curve_data <- result$curve_data
  
  # Extract BH.Q values
  bhq <- format(filtered_data$BH.Q[1], scientific = FALSE, digits = 2)
  
  # Adjust theme based on plot dimensions
  plot_width <- session$clientData$output_genePlot_width
  plot_height <- session$clientData$output_genePlot_height
  
  # Ensure plot dimensions are valid
  if (is.null(plot_width) || is.null(plot_height) || plot_width <= 0 || plot_height <= 0) {
    plot_width <- 500
    plot_height <- 375  # Maintain 3:4 aspect ratio
  } else {
    plot_height <- plot_width * 3 / 4  # Maintain 3:4 aspect ratio
  }
  
  # Calculate scaling factor
  scaling_factor <- min(plot_width, plot_height) / 200
  
  # Adjust sizes based on scaling factor
  text_size <- scaling_factor * 6
  line_size <- scaling_factor * 0.75
  point_size <- scaling_factor * 1.5
  annotate_text_size <- scaling_factor * 2
  
  # Plot both datasets and fits
  p <- ggplot() +
    theme_classic() +
    theme(
      plot.title = element_text(size = text_size * 1.5, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(size = text_size, face = "bold"),
      axis.text.y = element_text(size = text_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, title = gene) +
    geom_vline(xintercept = 6, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 12, linetype = "dashed", color = "grey") +
    geom_point(data = filtered_data, aes(x = as.numeric(time_point), y = mean), color = "black", size = point_size) +
    geom_errorbar(data = filtered_data, aes(x = as.numeric(time_point), ymin = mean - sd, ymax = mean + sd), width = 0.2, position = position_dodge(width = 0.9)) +
    scale_x_continuous(breaks = as.numeric(filtered_data$time_point), labels = c("CT4", "", "", "CT16", "", "", "CT28", "", "", "CT40", "", "")) +
    annotate("text", x = Inf, y = Inf, label = paste("bhq:", bhq), color = "black", hjust = 1.1, vjust = 2, size = annotate_text_size, fontface = "bold")
  
  # Conditionally add the sine fitting line
  if (show_fit) {
    p <- p + geom_line(data = curve_data, aes(x = t, y = y), color = "black", size = line_size)
  }
  
  # Conditionally add the simple line connecting points
  if (show_line) {
    p <- p + geom_line(data = filtered_data, aes(x = as.numeric(time_point), y = mean), color = "blue", size = line_size)
  }
  
  return(p)
}

fitted_plot_atlas <- function(data1, data2 = NULL, gene, period = 24, session,
                              show_fit = FALSE, show_line = TRUE, show_comparison = TRUE) {
  result1 <- get_gene_plot_data(data1, gene, period, t = NULL, include_fit = show_fit, quiet = TRUE)
  result2 <- NULL

  if (show_comparison && !is.null(data2)) {
    result2 <- get_gene_plot_data(data2, gene, period, t = NULL, include_fit = show_fit, quiet = TRUE)
  }

  # Skip genes that are not available in the raw expression data.
  if (is.null(result1)) {
    return(NULL)
  }

  filtered_data1 <- result1$filtered_data
  curve_data1    <- result1$curve_data
  bhq1 <- format(filtered_data1$BH.Q[1], scientific = FALSE, digits = 2)

  # Responsive sizing
  plot_width  <- session$clientData$output_genePlot_width
  plot_height <- session$clientData$output_genePlot_height

  if (is.null(plot_width) || is.null(plot_height) || plot_width <= 0 || plot_height <= 0) {
    plot_width  <- 500
    plot_height <- 375
  } else {
    plot_height <- plot_width * 3 / 4
  }

  scaling_factor    <- min(plot_width, plot_height) / 200
  text_size         <- scaling_factor * 6
  line_size         <- scaling_factor * 0.75
  point_size        <- scaling_factor * 1.5
  annotate_text_size <- scaling_factor * 2

  # X-axis labels at key ZT transitions (CT18=ZT18, CT24=ZT0, CT36=ZT12, CT48=ZT0, CT60=ZT12, CT64=ZT16)
  x_breaks <- seq(18, 64, by = 2)
  x_labels <- c("CT18", "", "", "CT24", "", "", "CT30", "", "", "CT36", "", "",
                "CT42", "", "", "CT48", "", "", "CT54", "", "", "CT60", "", "CT64")

  p <- ggplot() +
    theme_classic() +
    theme(
      plot.title = element_text(size = text_size * 1.5, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(size = text_size, face = "bold"),
      axis.text.y = element_text(size = text_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, title = gene) +
    # Dark periods: CT18-CT24 (start of data), CT36-CT48, CT60-CT64 (end of data)
    geom_rect(aes(xmin = 17, xmax = 24, ymin = -Inf, ymax = Inf), fill = "grey90") +
    geom_rect(aes(xmin = 36, xmax = 48, ymin = -Inf, ymax = Inf), fill = "grey90") +
    geom_rect(aes(xmin = 60, xmax = 65, ymin = -Inf, ymax = Inf), fill = "grey90") +
    geom_point(data = filtered_data1, aes(x = as.numeric(time_point), y = mean),
               color = "black", size = point_size) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels)

  if (show_comparison && !is.null(result2)) {
    filtered_data2 <- result2$filtered_data
    curve_data2    <- result2$curve_data
    bhq2 <- format(filtered_data2$BH.Q[1], scientific = FALSE, digits = 2)

    p <- p +
      geom_point(data = filtered_data2, aes(x = as.numeric(time_point), y = mean),
                 color = "navy", size = point_size) +
      annotate("text", x = Inf, y = Inf,
               label = paste("Tissue 1 BH.Q:", bhq1),
               color = "black", hjust = 1.1, vjust = 2,
               size = annotate_text_size, fontface = "bold") +
      annotate("text", x = Inf, y = Inf,
               label = paste("Tissue 2 BH.Q:", bhq2),
               color = "navy", hjust = 1.1, vjust = 3.5,
               size = annotate_text_size, fontface = "bold")
  } else {
    p <- p +
      annotate("text", x = Inf, y = Inf,
               label = paste("BH.Q:", bhq1),
               color = "black", hjust = 1.1, vjust = 2,
               size = annotate_text_size, fontface = "bold")
  }

  if (show_fit) {
    p <- p + geom_line(data = curve_data1, aes(x = t, y = y),
                       color = "black", size = line_size)
    if (show_comparison && !is.null(result2)) {
      p <- p + geom_line(data = curve_data2, aes(x = t, y = y),
                         color = "navy", size = line_size, linetype = "dashed")
    }
  }

  if (show_line) {
    p <- p + geom_line(data = filtered_data1, aes(x = as.numeric(time_point), y = mean),
                       color = "black", size = line_size)
    if (show_comparison && !is.null(result2)) {
      p <- p + geom_line(data = filtered_data2, aes(x = as.numeric(time_point), y = mean),
                         color = "navy", size = line_size)
    }
  }

  return(p)
}


# Plotting function for 24h ZT datasets (ZT00-ZT22, single sample per timepoint)
# Light on: ZT00-ZT12, light off: ZT12-ZT22
# Primary tissue: black; optional comparison tissue: navy
fitted_plot_baboon_atlas <- function(data1, data2 = NULL, gene, period = 24, session,
                           show_fit = FALSE, show_line = TRUE, show_comparison = TRUE) {
  result1 <- get_gene_plot_data(data1, gene, period, t = NULL, include_fit = show_fit, quiet = TRUE)
  result2 <- NULL

  if (show_comparison && !is.null(data2)) {
    result2 <- get_gene_plot_data(data2, gene, period, t = NULL, include_fit = show_fit, quiet = TRUE)
  }

  # Skip genes that are not available in the raw expression data.
  if (is.null(result1)) {
    return(NULL)
  }

  filtered_data1 <- result1$filtered_data
  curve_data1    <- result1$curve_data
  bhq1 <- format(filtered_data1$BH.Q[1], scientific = FALSE, digits = 2)

  # Responsive sizing
  plot_width  <- session$clientData$output_genePlot_width
  plot_height <- session$clientData$output_genePlot_height

  if (is.null(plot_width) || is.null(plot_height) || plot_width <= 0 || plot_height <= 0) {
    plot_width  <- 500
    plot_height <- 375
  } else {
    plot_height <- plot_width * 3 / 4
  }

  scaling_factor     <- min(plot_width, plot_height) / 200
  text_size          <- scaling_factor * 6
  line_size          <- scaling_factor * 0.75
  point_size         <- scaling_factor * 1.5
  annotate_text_size <- scaling_factor * 2

  x_breaks <- seq(0, 22, by = 2)
  x_labels <- sprintf("ZT%02d", x_breaks)

  p <- ggplot() +
    theme_classic() +
    theme(
      plot.title = element_text(size = text_size * 1.5, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(size = text_size, face = "bold"),
      axis.text.y = element_text(size = text_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, title = gene) +
    # Dark period: ZT12-ZT22
    geom_rect(aes(xmin = 12, xmax = 23, ymin = -Inf, ymax = Inf), fill = "grey90") +
    geom_point(data = filtered_data1, aes(x = as.numeric(time_point), y = mean),
               color = "black", size = point_size) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels)

  if (show_comparison && !is.null(result2)) {
    filtered_data2 <- result2$filtered_data
    curve_data2    <- result2$curve_data
    bhq2 <- format(filtered_data2$BH.Q[1], scientific = FALSE, digits = 2)

    p <- p +
      geom_point(data = filtered_data2, aes(x = as.numeric(time_point), y = mean),
                 color = "navy", size = point_size) +
      annotate("text", x = Inf, y = Inf,
               label = paste("Tissue 1 BH.Q:", bhq1),
               color = "black", hjust = 1.1, vjust = 2,
               size = annotate_text_size, fontface = "bold") +
      annotate("text", x = Inf, y = Inf,
               label = paste("Tissue 2 BH.Q:", bhq2),
               color = "navy", hjust = 1.1, vjust = 3.5,
               size = annotate_text_size, fontface = "bold")
  } else {
    p <- p +
      annotate("text", x = Inf, y = Inf,
               label = paste("BH.Q:", bhq1),
               color = "black", hjust = 1.1, vjust = 2,
               size = annotate_text_size, fontface = "bold")
  }

  if (show_fit) {
    p <- p + geom_line(data = curve_data1, aes(x = t, y = y),
                       color = "black", size = line_size)
    if (show_comparison && !is.null(result2)) {
      p <- p + geom_line(data = curve_data2, aes(x = t, y = y),
                         color = "navy", size = line_size, linetype = "dashed")
    }
  }

  if (show_line) {
    p <- p + geom_line(data = filtered_data1, aes(x = as.numeric(time_point), y = mean),
                       color = "black", size = line_size)
    if (show_comparison && !is.null(result2)) {
      p <- p + geom_line(data = filtered_data2, aes(x = as.numeric(time_point), y = mean),
                         color = "navy", size = line_size)
    }
  }

  return(p)
}


# Plotting function for 24h ZT datasets sampled every 4h in triplicate
# (ZT00, ZT04, ZT08, ZT12, ZT16, ZT20 → time_point indices 1–6)
# Light on: ZT00–ZT12 (x 1–3.5); dark: ZT12–ZT24 (x 3.5–6.5)
fitted_plot_zt24 <- function(data, gene, period = 24, t, session, show_fit = FALSE, show_line = TRUE) {

  result <- get_gene_plot_data(data, gene, period, t, include_fit = show_fit, quiet = TRUE)
  if (is.null(result)) return(NULL)

  filtered_data <- result$filtered_data
  curve_data    <- result$curve_data
  bhq <- format(filtered_data$BH.Q[1], scientific = FALSE, digits = 2)

  # Responsive sizing
  plot_width  <- session$clientData$output_genePlot_width
  plot_height <- session$clientData$output_genePlot_height

  if (is.null(plot_width) || is.null(plot_height) || plot_width <= 0 || plot_height <= 0) {
    plot_width  <- 500
    plot_height <- 375
  } else {
    plot_height <- plot_width * 3 / 4
  }

  scaling_factor     <- min(plot_width, plot_height) / 200
  text_size          <- scaling_factor * 6
  line_size          <- scaling_factor * 0.75
  point_size         <- scaling_factor * 1.5
  annotate_text_size <- scaling_factor * 2

  x_breaks <- 1:6
  x_labels <- c("ZT00", "ZT04", "ZT08", "ZT12", "ZT16", "ZT20")

  p <- ggplot() +
    theme_classic() +
    theme(
      plot.title  = element_text(size = text_size * 1.5, face = "bold", hjust = 0.5, vjust = 0.5),
      axis.text.x = element_text(size = text_size, face = "bold"),
      axis.text.y = element_text(size = text_size, face = "bold")
    ) +
    labs(x = NULL, y = NULL, title = gene) +
    # Dark period: ZT12–ZT24
    geom_rect(aes(xmin = 3.5, xmax = 6.5, ymin = -Inf, ymax = Inf), fill = "grey90") +
    geom_point(data = filtered_data, aes(x = as.numeric(time_point), y = mean),
               color = "black", size = point_size) +
    geom_errorbar(data = filtered_data,
                  aes(x = as.numeric(time_point), ymin = mean - sd, ymax = mean + sd),
                  width = 0.2) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels) +
    annotate("text", x = Inf, y = Inf, label = paste("bhq:", bhq),
             color = "black", hjust = 1.1, vjust = 2,
             size = annotate_text_size, fontface = "bold")

  if (show_fit) {
    p <- p + geom_line(data = curve_data, aes(x = t, y = y),
                       color = "black", size = line_size)
  }

  if (show_line) {
    p <- p + geom_line(data = filtered_data, aes(x = as.numeric(time_point), y = mean),
                       color = "blue", size = line_size)
  }

  return(p)
}


# ---- Background-specific enrichment helpers -------------------------------

ENRICHMENT_BACKGROUND_CODES <- c("dataset_tested", "all_string")
ENRICHMENT_BUNDLE_SCHEMA_VERSION <- 1L

validate_enrichment_variant <- function(variant, dataset_id, background) {
  context <- paste0(
    "Enrichment background '", background, "' for dataset '", dataset_id, "'"
  )
  required_components <- c(
    "node_titles", "cluster_descriptions", "cluster_to_genes"
  )

  if (!is.list(variant)) {
    stop(context, " must be a list.", call. = FALSE)
  }
  missing_components <- setdiff(required_components, names(variant))
  if (length(missing_components) > 0) {
    stop(
      context, " is missing: ", paste(missing_components, collapse = ", "),
      ".", call. = FALSE
    )
  }

  node_titles <- variant$node_titles
  if (!is.data.frame(node_titles) ||
      !all(c("id", "title") %in% names(node_titles))) {
    stop(context, "$node_titles must contain character columns id and title.",
         call. = FALSE)
  }
  if (!is.character(node_titles$id) || !is.character(node_titles$title) ||
      anyNA(node_titles$id) || any(!nzchar(node_titles$id)) ||
      anyNA(node_titles$title) ||
      anyDuplicated(node_titles$id)) {
    stop(context, "$node_titles contains invalid or duplicate node IDs.",
         call. = FALSE)
  }

  descriptions <- variant$cluster_descriptions
  description_columns <- c("category", "cluster", "description", "fdr")
  if (!is.data.frame(descriptions) ||
      !all(description_columns %in% names(descriptions))) {
    stop(
      context, "$cluster_descriptions must contain: ",
      paste(description_columns, collapse = ", "), ".", call. = FALSE
    )
  }
  if (!all(vapply(descriptions[c("category", "cluster", "description")],
                  is.character, logical(1))) ||
      !is.numeric(descriptions$fdr)) {
    stop(context, "$cluster_descriptions has invalid column types.",
         call. = FALSE)
  }
  valid_categories <- c("Process", "KEGG", "RCTM", "WikiPathways")
  if (anyNA(descriptions$category) ||
      any(!descriptions$category %in% valid_categories) ||
      anyNA(descriptions$cluster) || any(!nzchar(descriptions$cluster)) ||
      anyNA(descriptions$description) ||
      any(!is.finite(descriptions$fdr)) ||
      any(descriptions$fdr < 0 | descriptions$fdr > 1) ||
      anyDuplicated(descriptions[c("category", "cluster")])) {
    stop(context, "$cluster_descriptions contains invalid or duplicate terms.",
         call. = FALSE)
  }

  cluster_to_genes <- variant$cluster_to_genes
  gene_columns <- c("category", "cluster", "genes")
  if (!is.data.frame(cluster_to_genes) ||
      !all(gene_columns %in% names(cluster_to_genes))) {
    stop(
      context, "$cluster_to_genes must contain: ",
      paste(gene_columns, collapse = ", "), ".", call. = FALSE
    )
  }
  if (!is.character(cluster_to_genes$category) ||
      !is.character(cluster_to_genes$cluster) ||
      !is.list(cluster_to_genes$genes) ||
      !all(vapply(cluster_to_genes$genes, is.character, logical(1))) ||
      anyNA(cluster_to_genes$category) ||
      any(!cluster_to_genes$category %in% valid_categories) ||
      anyNA(cluster_to_genes$cluster) ||
      any(!nzchar(cluster_to_genes$cluster)) ||
      any(vapply(cluster_to_genes$genes, anyNA, logical(1))) ||
      any(vapply(cluster_to_genes$genes, function(genes) {
        !length(genes) || any(!nzchar(genes)) || anyDuplicated(genes) > 0L
      }, logical(1))) ||
      anyDuplicated(cluster_to_genes[c("category", "cluster")])) {
    stop(context, "$cluster_to_genes has invalid columns or duplicate terms.",
         call. = FALSE)
  }

  description_keys <- paste(descriptions$category, descriptions$cluster, sep = "\r")
  gene_keys <- paste(cluster_to_genes$category, cluster_to_genes$cluster, sep = "\r")
  if (!setequal(description_keys, gene_keys)) {
    stop(
      context,
      " has inconsistent term keys in cluster_descriptions and cluster_to_genes.",
      call. = FALSE
    )
  }

  # Preserve the legacy names used by the Explorer and Comparison modules.
  descriptions$cluster_description <- descriptions$description
  list(
    node_titles = node_titles,
    cluster_descriptions = descriptions,
    cluster_to_genes = cluster_to_genes
  )
}

select_enrichment_background <- function(bundle, dataset_id, background) {
  if (!background %in% ENRICHMENT_BACKGROUND_CODES) {
    stop(
      "Unknown enrichment background '", background, "'. Expected one of: ",
      paste(ENRICHMENT_BACKGROUND_CODES, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!is.list(bundle) ||
      !identical(as.integer(bundle$schema_version),
                 ENRICHMENT_BUNDLE_SCHEMA_VERSION) ||
      !is.list(bundle$metadata) || !is.list(bundle$backgrounds)) {
    stop(
      "Invalid enrichment-background bundle for dataset '", dataset_id,
      "': expected schema version ", ENRICHMENT_BUNDLE_SCHEMA_VERSION, ".",
      call. = FALSE
    )
  }
  missing_backgrounds <- setdiff(
    ENRICHMENT_BACKGROUND_CODES, names(bundle$backgrounds)
  )
  if (length(missing_backgrounds) > 0) {
    stop(
      "Enrichment-background bundle for dataset '", dataset_id,
      "' is missing: ", paste(missing_backgrounds, collapse = ", "), ".",
      call. = FALSE
    )
  }

  validate_enrichment_variant(
    bundle$backgrounds[[background]], dataset_id, background
  )
}

apply_enrichment_node_titles <- function(nodes, node_titles, dataset_id,
                                         background) {
  if (!is.data.frame(nodes) || !"id" %in% names(nodes)) {
    stop("Network nodes for dataset '", dataset_id,
         "' do not contain an id column.", call. = FALSE)
  }
  title_index <- match(as.character(nodes$id), node_titles$id)
  if (anyNA(title_index)) {
    missing_ids <- as.character(nodes$id[is.na(title_index)])
    stop(
      "Enrichment background '", background, "' for dataset '", dataset_id,
      "' is missing tooltip titles for ", length(missing_ids), " node(s): ",
      paste(utils::head(missing_ids, 5), collapse = ", "),
      if (length(missing_ids) > 5) ", ..." else "", ".",
      call. = FALSE
    )
  }
  nodes$title <- node_titles$title[title_index]
  nodes
}

normalise_legacy_enrichment <- function(cluster_descriptions,
                                        cluster_to_genes) {
  required_description_columns <- c(
    "cluster", "cluster_description", "category"
  )
  if (!is.data.frame(cluster_descriptions) ||
      !all(required_description_columns %in% names(cluster_descriptions))) {
    stop("Legacy cluster descriptions have an unsupported schema.",
         call. = FALSE)
  }
  if (!is.data.frame(cluster_to_genes) ||
      !all(c("cluster", "genes") %in% names(cluster_to_genes))) {
    stop("Legacy cluster-to-gene data have an unsupported schema.",
         call. = FALSE)
  }

  cluster_descriptions$description <-
    cluster_descriptions$cluster_description
  cluster_descriptions$fdr <- NA_real_

  if (!"category" %in% names(cluster_to_genes)) {
    category_lookup <- unique(
      cluster_descriptions[c("cluster", "category")]
    )
    cluster_to_genes <- dplyr::left_join(
      cluster_to_genes, category_lookup, by = "cluster"
    )
  }

  list(
    cluster_descriptions = cluster_descriptions,
    cluster_to_genes = cluster_to_genes
  )
}

format_enrichment_term_choices <- function(cluster_descriptions, category) {
  if (!is.data.frame(cluster_descriptions) ||
      !all(c("category", "cluster", "cluster_description") %in%
           names(cluster_descriptions))) {
    return(character(0))
  }

  terms <- cluster_descriptions[
    cluster_descriptions$category == category,
    , drop = FALSE
  ]
  if (nrow(terms) == 0) return(character(0))

  terms <- terms[!duplicated(terms[c("category", "cluster")]), , drop = FALSE]
  fdr <- if ("fdr" %in% names(terms)) {
    suppressWarnings(as.numeric(terms$fdr))
  } else {
    rep(NA_real_, nrow(terms))
  }
  fdr[is.na(fdr)] <- Inf
  order_index <- order(
    fdr, tolower(terms$cluster_description), terms$cluster,
    na.last = TRUE
  )
  terms <- terms[order_index, , drop = FALSE]

  labels <- terms$cluster_description
  duplicate_labels <- duplicated(labels) | duplicated(labels, fromLast = TRUE)
  labels[duplicate_labels] <- paste0(
    labels[duplicate_labels], " [", terms$cluster[duplicate_labels], "]"
  )
  stats::setNames(as.character(terms$cluster), labels)
}

# Function to limit GO Biological Process terms in tooltip
limit_go_terms_in_tooltip <- function(tooltip_text, max_terms = 5) {
  if (is.na(tooltip_text) || tooltip_text == "") return(tooltip_text)
  
  # Check if the tooltip contains GO Biological Process section
  if (!grepl("<b>GO Biological Process:</b>", tooltip_text)) {
    return(tooltip_text)  # No BP terms, return unchanged
  }
  
  # Extract the GO Biological Process section using regex
  bp_section_pattern <- "<b>GO Biological Process:</b>([^<]*)"
  bp_section_match <- regexpr(bp_section_pattern, tooltip_text, perl = TRUE)
  
  if (bp_section_match == -1) {
    return(tooltip_text)  # Pattern not found, return unchanged
  }
  
  bp_section <- regmatches(tooltip_text, bp_section_match)
  
  # Extract all terms from the section
  terms_part <- sub("<b>GO Biological Process:</b>\\s*", "", bp_section)
  terms <- trimws(strsplit(terms_part, ",\\s*")[[1]])
  
  # Check if we need to limit the terms
  if (length(terms) <= max_terms) {
    return(tooltip_text)  # No need to limit
  }
  
  # Keep only the first max_terms
  limited_terms <- terms[1:max_terms]
  num_remaining <- length(terms) - max_terms
  
  # Create the new section with note about additional terms
  new_section <- paste0(
    "<b>GO Biological Process:</b> ",
    paste(limited_terms, collapse = ", "),
    ", ... and ", num_remaining, " more terms"
  )
  
  # Replace the original section with the new one
  result <- sub(bp_section_pattern, new_section, tooltip_text, perl = TRUE)
  
  return(result)
}
