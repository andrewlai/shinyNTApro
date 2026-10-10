# --- shinyNTApro ---
options(shiny.maxRequestSize = 100 * 1024^2) # 100 MB total request limit

# Libraries
library(shiny); library(shinythemes); library(readr); library(dplyr)
library(stringr); library(ggplot2); library(DT)
library(viridis); library(broom); library(ggpubr); library(shinyjs)
library(car); library(rmarkdown); library(lubridate); library(purrr)

# Helpers
`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
str_starts_compat <- function(x, pattern) grepl(paste0("^", pattern), x)

# Parse QC timestamps deterministically in the laboratory's local timezone.
# Using an explicit timezone prevents early-morning Brisbane records from being
# shifted to the previous UTC calendar date during daily grouping.
DEFAULT_QC_TIMEZONE <- "Australia/Brisbane"

parse_qc_datetime <- function(x, tz = DEFAULT_QC_TIMEZONE) {
  x <- stringr::str_squish(as.character(x))
  formats <- c(
    "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M",
    "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M",
    "%H:%M:%S %d/%m/%Y", "%H:%M %d/%m/%Y",
    "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M",
    "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M",
    "%H:%M:%S %m/%d/%Y", "%H:%M %m/%d/%Y",
    "%Y-%m-%d", "%d/%m/%Y", "%Y/%m/%d", "%m/%d/%Y"
  )

  parse_one <- function(value) {
    if (is.na(value) || !nzchar(value)) return(NA_real_)
    for (fmt in formats) {
      parsed <- suppressWarnings(as.POSIXct(strptime(value, format = fmt, tz = tz)))
      if (!is.na(parsed)) return(as.numeric(parsed))
    }
    NA_real_
  }

  values <- vapply(x, parse_one, numeric(1))
  as.POSIXct(values, origin = "1970-01-01", tz = tz)
}

clean_name_vector <- function(x) {
  x <- as.character(x)
  x <- gsub("mL", "ml", x, fixed = TRUE)
  x <- gsub("uL", "ul", x, fixed = TRUE)
  x <- gsub("µL", "ul", x, fixed = TRUE)
  x <- gsub("([A-Z]+)([A-Z][a-z])", "\\1_\\2", x, perl = TRUE)
  x <- gsub("([a-z0-9])([A-Z])", "\\1_\\2", x, perl = TRUE)
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x[!nzchar(x)] <- "column"
  make.unique(x, sep = "_")
}

clean_names_df <- function(df) {
  names(df) <- clean_name_vector(names(df))
  df
}

safe_file_stub <- function(x) {
  x <- clean_name_vector(x %||% "file")
  x[1]
}

cohens_d_unpooled <- function(data, outcome, group_col = "group") {
  groups <- levels(factor(data[[group_col]]))
  if (length(groups) != 2) return(NULL)
  values1 <- data[[outcome]][data[[group_col]] == groups[1]]
  values2 <- data[[outcome]][data[[group_col]] == groups[2]]
  values1 <- values1[!is.na(values1)]; values2 <- values2[!is.na(values2)]
  if (length(values1) < 2 || length(values2) < 2) return(NULL)
  denominator <- sqrt((stats::var(values1) + stats::var(values2)) / 2)
  estimate <- if (is.finite(denominator) && denominator > 0) (mean(values1) - mean(values2)) / denominator else NA_real_
  tibble::tibble(effect = "Cohen's d (unpooled SD)", estimate = estimate, group1 = groups[1], group2 = groups[2])
}

rank_biserial_from_wilcox <- function(test_obj, data, group_col = "group") {
  groups <- levels(factor(data[[group_col]]))
  if (length(groups) != 2) return(NULL)
  n1 <- sum(data[[group_col]] == groups[1], na.rm = TRUE)
  n2 <- sum(data[[group_col]] == groups[2], na.rm = TRUE)
  if (n1 == 0 || n2 == 0) return(NULL)
  u <- unname(test_obj$statistic)
  tibble::tibble(effect = "Rank-biserial correlation", estimate = 2 * u / (n1 * n2) - 1, group1 = groups[1], group2 = groups[2])
}

anova_effect_sizes <- function(aov_obj) {
  tab <- broom::tidy(aov_obj)
  group_row <- tab[tab$term == "group", , drop = FALSE]
  residual_row <- tab[tab$term == "Residuals", , drop = FALSE]
  if (nrow(group_row) != 1 || nrow(residual_row) != 1) stop("One-way ANOVA terms unavailable for effect sizes.")
  total_ss <- sum(tab$sumsq)
  mse <- residual_row$sumsq / residual_row$df
  eta <- if (is.finite(total_ss) && total_ss > 0) group_row$sumsq / total_ss else NA_real_
  omega <- if (is.finite(total_ss + mse) && total_ss + mse > 0) {
    max(0, (group_row$sumsq - group_row$df * mse) / (total_ss + mse))
  } else NA_real_
  tibble::tibble(effect = c("Eta squared", "Omega squared"), estimate = c(eta, omega),
    note = c("SS(group) / SS(total).", "Bias-corrected estimate; negative estimates truncated to zero."))
}

welch_effect_sizes <- function(test_obj) {
  f <- unname(test_obj$statistic)
  df1 <- unname(test_obj$parameter[1]); df2 <- unname(test_obj$parameter[2])
  if (length(f) != 1 || !all(is.finite(c(f, df1, df2))) || f < 0 || df1 <= 0 || df2 <= 0) {
    stop("Finite Welch F statistic and positive degrees of freedom are required for effect sizes.")
  }
  # F-to-effect-size conversions; not SS-based variance-explained estimates.
  tibble::tibble(
    effect = c("Partial eta squared (approximate, Welch F)", "Partial omega squared (approximate, Welch F)"),
    estimate = c(f * df1 / (f * df1 + df2), max(0, (f - 1) * df1 / (f * df1 + df2 + 1))),
    note = c("F*df1 / (F*df1 + df2). Approximation using Welch degrees of freedom; not ordinary ANOVA variance explained.",
             "(F-1)*df1 / (F*df1 + df2 + 1). Approximation; negative estimates truncated to zero."))
}

kruskal_eta_squared_h <- function(test_obj, n, k) {
  estimate <- if (n > k) (unname(test_obj$statistic) - k + 1) / (n - k) else NA_real_
  if (is.finite(estimate)) estimate <- max(0, estimate)
  tibble::tibble(effect = "Eta squared (H-based)", estimate = estimate)
}

dunn_holm_test <- function(x, g) {
  keep <- !is.na(x) & !is.na(g)
  x <- x[keep]; g <- droplevels(factor(g[keep]))
  if (nlevels(g) < 2) return(tibble::tibble())
  n_total <- length(x)
  ranks <- rank(x, ties.method = "average")
  group_summary <- tibble::tibble(group = g, rank = ranks) %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(n = dplyr::n(), mean_rank = mean(rank), .groups = "drop")
  tie_sizes <- as.numeric(table(x))
  tie_term <- sum(tie_sizes^3 - tie_sizes)
  base_variance <- n_total * (n_total + 1) / 12
  if (n_total > 1) base_variance <- base_variance - tie_term / (12 * (n_total - 1))
  comparisons <- utils::combn(seq_len(nrow(group_summary)), 2, simplify = FALSE)
  result <- purrr::map_dfr(comparisons, function(idx) {
    i <- idx[1]; j <- idx[2]
    standard_error <- sqrt(base_variance * (1 / group_summary$n[i] + 1 / group_summary$n[j]))
    z_value <- if (is.finite(standard_error) && standard_error > 0) (group_summary$mean_rank[i] - group_summary$mean_rank[j]) / standard_error else NA_real_
    tibble::tibble(
      group1 = as.character(group_summary$group[i]),
      group2 = as.character(group_summary$group[j]),
      statistic = z_value,
      p.raw = if (is.na(z_value)) NA_real_ else 2 * stats::pnorm(-abs(z_value))
    )
  })
  result %>%
    dplyr::mutate(p.adj = stats::p.adjust(p.raw, method = "holm"), p.adj.signif = significance_code(p.adj)) %>%
    dplyr::select(group1, group2, statistic, p.adj, p.adj.signif)
}

make_sample_key <- function(source_file, sample_name, well_id) {
  paste(
    as.character(source_file),
    dplyr::coalesce(as.character(sample_name), "Unnamed sample"),
    dplyr::coalesce(as.character(well_id), "Unknown well"),
    sep = " | "
  )
}

significance_code <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p <= 0.001 ~ "***",
    p <= 0.01 ~ "**",
    p <= 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

split_tukey_contrast <- function(contrast, group_levels) {
  matches <- list()
  for (group2 in group_levels) {
    for (group1 in group_levels) {
      if (!identical(group1, group2) && identical(contrast, paste0(group2, "-", group1))) {
        matches[[length(matches) + 1L]] <- c(group1 = group1, group2 = group2)
      }
    }
  }
  if (length(matches) != 1L) stop("Tukey contrast could not be mapped uniquely: ", contrast,
    ". Rename groups to avoid ambiguous hyphenated pair names.")
  matches[[1]]
}

games_howell_test <- function(data, outcome, group_col = "group", conf.level = 0.95) {
  summary_df <- data %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarise(
      n = sum(!is.na(.data[[outcome]])),
      mean = mean(.data[[outcome]], na.rm = TRUE),
      variance = stats::var(.data[[outcome]], na.rm = TRUE),
      .groups = "drop"
    )

  invalid <- is.na(summary_df[[group_col]]) | summary_df$n < 2 |
    !is.finite(summary_df$mean) | !is.finite(summary_df$variance) | summary_df$variance <= 0
  if (nrow(summary_df) < 2 || any(invalid)) {
    stop("Games-Howell requires at least two groups, each with n >= 2 and positive finite variance; no groups were dropped.")
  }

  group_names <- as.character(summary_df[[group_col]])
  comparisons <- utils::combn(seq_len(nrow(summary_df)), 2, simplify = FALSE)
  k <- nrow(summary_df)
  alpha <- 1 - conf.level

  purrr::map_dfr(comparisons, function(idx) {
    i <- idx[1]; j <- idx[2]
    vi <- summary_df$variance[i] / summary_df$n[i]
    vj <- summary_df$variance[j] / summary_df$n[j]
    se_diff <- sqrt(vi + vj)
    estimate <- summary_df$mean[i] - summary_df$mean[j]

    if (!is.finite(se_diff) || se_diff <= 0) {
      return(tibble::tibble(
        group1 = group_names[i], group2 = group_names[j], estimate = estimate,
        conf.low = NA_real_, conf.high = NA_real_, statistic = NA_real_,
        df = NA_real_, p.adj = NA_real_, p.adj.signif = NA_character_
      ))
    }

    df_value <- (vi + vj)^2 / ((vi^2 / (summary_df$n[i] - 1)) + (vj^2 / (summary_df$n[j] - 1)))
    q_value <- sqrt(2) * abs(estimate) / se_diff
    p_value <- stats::ptukey(q_value, nmeans = k, df = df_value, lower.tail = FALSE)
    critical_q <- stats::qtukey(1 - alpha, nmeans = k, df = df_value)
    half_width <- critical_q * sqrt((vi + vj) / 2)

    tibble::tibble(
      group1 = group_names[i], group2 = group_names[j], estimate = estimate,
      conf.low = estimate - half_width, conf.high = estimate + half_width,
      statistic = q_value, df = df_value, p.adj = p_value,
      p.adj.signif = significance_code(p_value)
    )
  })
}

# Within-group diagnostics for the selected parameter, after missing-value filtering.
assess_group_normality <- function(data, metric) {
  purrr::map_dfr(levels(droplevels(factor(data$group))), function(group_name) {
    values <- data[[metric]][as.character(data$group) == group_name]
    values <- values[is.finite(values)]
    n <- length(values)
    reason <- if (n < 3) "Unavailable: fewer than 3 values" else if (n > 5000) {
      "Unavailable: more than 5000 values"
    } else if (length(unique(values)) < 2) "Unavailable: all values identical" else ""
    p <- if (nzchar(reason)) NA_real_ else tryCatch(
      stats::shapiro.test(values)$p.value, error = function(e) NA_real_)
    if (!is.finite(p) && !nzchar(reason)) reason <- "Unavailable: Shapiro-Wilk could not be computed"
    tibble::tibble(group = group_name, n = n, shapiro_p = p,
      diagnostic = if (nzchar(reason)) reason else if (p <= 0.05) "Departure detected" else "No departure detected")
  })
}

# Transparent screening rules; this is not a guarantee of test suitability.
statistical_test_choices <- function(n_groups = NULL) {
  choices <- list("Auto: select using diagnostics" = "Auto")
  means <- ranks <- character()
  if (is.null(n_groups) || n_groups == 2L) {
    means <- c(means, "Welch t-test (2 groups)" = "t.test")
    ranks <- c(ranks, "Wilcoxon rank-sum (2 groups)" = "wilcox.test")
  }
  if (is.null(n_groups) || n_groups > 2L) {
    means <- c(means, "ANOVA (equal variances)" = "aov", "Welch ANOVA" = "welch_aov")
    ranks <- c(ranks, "Kruskal-Wallis" = "kruskal.test")
  }
  if (length(means)) choices[["Manual: compare means"]] <- means
  if (length(ranks)) choices[["Manual: compare ranks"]] <- ranks
  choices
}

select_statistical_test <- function(choice, num_groups, min_n, shapiro_p, levene_p) {
  labels <- c("t.test" = "Welch t-test", "wilcox.test" = "Wilcoxon rank-sum test",
              "aov" = "One-way ANOVA", "welch_aov" = "Welch ANOVA",
              "kruskal.test" = "Kruskal-Wallis")
  automatic <- identical(choice, "Auto")
  flagged <- is.finite(shapiro_p) & shapiro_p <= 0.05
  normality_flag <- any(flagged)
  incomplete <- length(shapiro_p) != num_groups || any(!is.finite(shapiro_p))
  equal_variance_screen <- is.finite(levene_p) && levene_p > 0.05
  error <- NULL
  if (!automatic) {
    actual <- choice
  } else if (normality_flag) {
    actual <- if (num_groups == 2) "wilcox.test" else "kruskal.test"
  } else if (incomplete) {
    actual <- NA_character_
    error <- "Auto could not complete the within-group normality checks. Review the diagnostic table in full statistical details, then select a test manually for your study design. No test was run."
  } else if (num_groups == 2) {
    actual <- "t.test"
  } else {
    actual <- if (equal_variance_screen) "aov" else "welch_aov"
  }
  label <- if (is.na(actual)) "Not selected" else unname(labels[actual]) %||% actual
  group_labels <- names(shapiro_p)
  if (is.null(group_labels)) group_labels <- paste("Group", seq_along(shapiro_p))
  normality_detail <- if (normality_flag) {
    paste0("Within-group Shapiro-Wilk detected a departure from normality (p <= 0.05) in: ",
           paste(group_labels[flagged], collapse = ", "), ".")
  } else if (incomplete) {
    "At least one within-group Shapiro-Wilk check was unavailable; absence of a detected departure cannot be established for every group."
  } else {
    "All within-group Shapiro-Wilk p-values were above 0.05. No departure was detected; this does not establish normality."
  }
  if (normality_flag && incomplete) normality_detail <- paste(normality_detail,
    "Other normality checks were unavailable; selection used the detected departure in an assessable group.")
  variance_detail <- if (!is.finite(levene_p)) {
    "Median-centred Levene's test was unavailable. Equal variances were not assumed."
  } else if (equal_variance_screen) {
    sprintf("Median-centred Levene's test did not detect unequal variances (p = %.4g > 0.05); equality is not established.", levene_p)
  } else sprintf("Median-centred Levene's test flagged unequal variances (p = %.4g <= 0.05).", levene_p)

  reason <- if (!automatic) {
    paste0(label, " was selected manually; diagnostic checks did not override that choice.")
  } else if (!is.null(error)) {
    "Automatic selection stopped because normality screening was incomplete and no assessable group triggered the rank-test rule."
  } else if (normality_flag) {
    paste0(label, " was selected for ", num_groups, " independent groups because at least one within-group normality check flagged a departure.")
  } else if (num_groups == 2) {
    "Welch t-test was selected for two independent groups because neither normality check flagged a departure. It compares group means without assuming equal variances."
  } else if (equal_variance_screen) {
    "One-way ANOVA was selected because the within-group normality checks and variance check did not flag a departure. It compares group means assuming equal variances."
  } else {
    "Welch ANOVA was selected because no within-group normality check flagged a departure, while the variance check flagged unequal variances or was unavailable. It compares group means without assuming equal variances."
  }
  target_note <- if (actual %in% c("wilcox.test", "kruskal.test")) {
    "This is a rank-based comparison, not a test of arithmetic means. A location/median interpretation requires comparable distribution shapes."
  } else "Mean comparisons remain sensitive to skewness and outliers, particularly with small samples."
  size_note <- sprintf("Minimum group size for this parameter: n = %d independent samples after removing missing values.", min_n)
  if (min_n < 5) size_note <- paste(size_note,
    "Normality checks have limited sensitivity with very small groups. Small sample size alone did not trigger a rank test.")
  gate_note <- if (automatic) {
    "Auto is a screening heuristic using unadjusted diagnostic p-values, not proof of the correct test. Selecting a test from the same data can alter its error rate; reported p-values do not account for that selection."
  } else "The analysis requires independent samples; diagnostic tests do not verify independence."
  posthoc <- switch(actual %||% "none", "t.test" = "No additional post-hoc test (two groups)",
    "wilcox.test" = "No additional post-hoc test (two groups)", "aov" = "Tukey HSD",
    "welch_aov" = "Games-Howell", "kruskal.test" = "Dunn (Holm adjustment)", "None")
  list(test_name = actual, test_label = label, error = error, post_hoc_method = posthoc,
       selection_heading = if (automatic) "Auto-selection rationale" else "Manual-selection rationale",
       normality_status = if (normality_flag) "Departure detected" else if (incomplete) "Incomplete" else "No departure detected",
       selection_reason = paste(reason, normality_detail, variance_detail, size_note, target_note, gate_note))
}

# Constants
QC_LOG_FILE_BASENAME <- "nta_qc_log.csv"
DEFAULT_QC_TARGET_SIZE <- 100; DEFAULT_QC_TOLERANCE_PERCENT <- 10
DEFAULT_INTRA_CV_LIMIT <- 20; DEFAULT_INTER_CV_LIMIT <- 25
REPORT_TEMPLATE_FILENAME <- file.path("templates", "report_template.Rmd")
EXAMPLE_SUMMARY_FILENAME <- file.path("data", "examples", "example_nta_summary.csv")
DEFAULT_PLOT_TYPE <- "points"; DEFAULT_USE_LOG_SCALE <- FALSE; DEFAULT_COLOR_BY_SAMPLE <- TRUE
MAX_UPLOAD_FILES <- 20
MAX_UPLOAD_MB_PER_FILE <- 25
GITHUB_REPO_URL <- "" # Add your repository URL to show a GitHub link in the app.

# --- Parsing ---
# Read individual CSV rows with quoting intact (sample descriptions may contain commas).
nta_csv_cells <- function(line) {
  trimws(as.character(unlist(utils::read.csv(
    text = line, header = FALSE, colClasses = "character", check.names = FALSE,
    na.strings = "", strip.white = TRUE, comment.char = ""
  ), use.names = FALSE)))
}

# A manual-loading ExperimentSummary describes ONE sample with repeat captures.
# Use the instrument's exported Average/Standard Error, never treat captures as
# independent samples or recompute rounded instrument results.
parse_nta_manual_summary <- function(raw_lines, metadata) {
  labels <- trimws(sub(",.*$", "", raw_lines))
  section <- function(label) {
    start <- which(labels == label)
    if (length(start) != 1L) stop("Expected one ", label, " section in the manual-loading summary.")
    following <- which(startsWith(labels, "[") & seq_along(labels) > start)
    end <- if (length(following)) following[1] - 1L else length(labels)
    seq.int(start, end)
  }
  row_index <- function(label, rows, required = TRUE) {
    idx <- rows[labels[rows] == label]
    if (length(idx) > 1L) stop("Ambiguous row in manual-loading summary: ", label)
    if (!length(idx)) {
      if (required) stop("Missing row in manual-loading summary: ", label)
      return(NA_integer_)
    }
    idx
  }
  number <- function(value) {
    suppressWarnings(readr::parse_double(value, na = c("", "NA", "#DIV/0!", "Not recorded")))
  }
  main_error <- NULL
  distribution_error <- NULL
  size_rows <- NULL
  sample_name <- NA_character_
  well_id <- "Not applicable"

  tidy_data <- tryCatch({
    result_rows <- section("[Results]")
    size_rows <- section("[Size Data]")
    weighting <- nta_csv_cells(raw_lines[row_index("Weighting", size_rows)])[2]
    if (is.na(weighting) || tolower(weighting) != "number") {
      stop("Manual-loading import currently requires Number-weighted size data.")
    }
    metadata$analysis_method <- nta_csv_cells(raw_lines[row_index("Analysis Method", size_rows)])[2]
    metadata$weighting <- weighting

    description_row <- row_index("Sample Description", seq_along(labels), required = FALSE)
    if (!is.na(description_row)) sample_name <- nta_csv_cells(raw_lines[description_row])[2]
    if (is.na(sample_name) || !nzchar(sample_name)) sample_name <- metadata$experiment_name
    if (is.na(sample_name) || !nzchar(sample_name)) stop("Sample Description and Experiment Name are both missing.")

    result_header <- nta_csv_cells(raw_lines[result_rows[1]])
    size_header <- nta_csv_cells(raw_lines[row_index("Filename", size_rows)])
    average_columns <- function(header) {
      avg <- which(tolower(header) == "average")
      se <- which(tolower(header) == "standard error")
      if (length(avg) != 1L || length(se) != 1L) stop("Average/Standard Error columns missing or ambiguous.")
      c(average = avg, standard_error = se)
    }
    result_cols <- average_columns(result_header)
    size_cols <- average_columns(size_header)
    capture_cols <- setdiff(seq_along(size_header), c(1L, size_cols))
    capture_cols <- capture_cols[!is.na(size_header[capture_cols]) & nzchar(size_header[capture_cols])]
    if (!length(capture_cols)) stop("No capture filenames in the manual-loading summary.")
    metadata$capture_count <- as.character(length(capture_cols))
    metadata$import_format <- "Manual loading (one sample per file; exported averages)"

    parsed <- tibble::tibble(well_id = well_id, sample_name = sample_name)
    add_metric <- function(label, rows, columns, prefix, suffix = "", required = TRUE) {
      idx <- row_index(label, rows, required)
      values <- if (is.na(idx)) rep(NA_real_, 2L) else number(nta_csv_cells(raw_lines[idx])[columns])
      parsed[[paste0(prefix, "_average", suffix)]] <<- values[1]
      parsed[[paste0(prefix, "_standard_error", suffix)]] <<- values[2]
    }
    add_metric("Concentration (Particles / ml)", result_rows, result_cols, "concentration", "_particles_ml")
    add_metric("Particles per frame", result_rows, result_cols, "particles_per_frame", required = FALSE)
    add_metric("Centres per frame", result_rows, result_cols, "centres_per_frame", required = FALSE)
    for (metric in c("Mean", "Mode", "SD", "D10", "D50", "D90")) {
      add_metric(metric, size_rows, size_cols, tolower(metric), required = metric %in% c("Mean", "Mode"))
    }
    if (!any(is.finite(unlist(parsed[c("concentration_average_particles_ml", "mean_average", "mode_average")])))) {
      stop("No numeric concentration, mean or mode averages in the manual-loading summary.")
    }
    parsed$time_captured <- metadata$time_captured
    parsed
  }, error = function(e) {
    main_error <<- conditionMessage(e)
    NULL
  })

  distribution_data <- NULL
  if (!is.null(tidy_data)) {
    distribution_data <- tryCatch({
      header_idx <- row_index("Bin centre (nm)", size_rows, required = FALSE)
      if (is.na(header_idx)) {
        if (any(labels[size_rows] == "Graph Data")) stop("Graph Data is present but its bin header is missing.")
        NULL
      } else {
        header <- nta_csv_cells(raw_lines[header_idx])
        avg_col <- which(tolower(header) == "concentration average")
        se_col <- which(tolower(header) == "standard error")
        if (length(avg_col) != 1L || length(se_col) != 1L) stop("Distribution average/error columns missing or ambiguous.")
        rows <- size_rows[size_rows > header_idx]
        numeric_bin <- grepl("^[0-9]+(\\.[0-9]+)?$", labels[rows])
        first_other <- which(!numeric_bin)
        if (length(first_other)) rows <- head(rows, first_other[1] - 1L)
        if (!length(rows)) stop("No numeric distribution rows found.")
        dist <- readr::read_csv(
          I(paste(raw_lines[c(header_idx, rows)], collapse = "\n")),
          col_types = readr::cols(.default = "c"), name_repair = "minimal",
          show_col_types = FALSE, progress = FALSE
        )
        bins <- number(dist[[1]])
        averages <- number(dist[[avg_col]])
        if (!all(is.finite(bins)) || !any(is.finite(averages))) stop("No usable distribution values.")
        tibble::tibble(
          WellID_dist = well_id, SampleName_dist = sample_name,
          BinCentre_nm = bins, ConcentrationAverage_dist = averages,
          StandardError_dist = number(dist[[se_col]])
        )
      }
    }, error = function(e) {
      distribution_error <<- conditionMessage(e)
      NULL
    })
  }
  list(data = tidy_data, metadata = metadata, distribution_data = distribution_data,
       main_error = main_error, distribution_error = distribution_error,
       status = if (!is.null(tidy_data)) "success" else "failed")
}

parse_nta_summary <- function(file_path) {
  message("Reading summary file: ", basename(file_path))
  if (!file.exists(file_path)) stop("File not found: ", file_path)

  raw_lines <- readr::read_lines(file_path)
  if (!length(raw_lines)) stop("File is empty: ", basename(file_path))

  extract_metadata_value <- function(lines, pattern) {
    line_index <- which(stringr::str_detect(lines, pattern))
    if (!length(line_index)) return(NA_character_)
    line <- lines[line_index[1]]
    extracted <- if (grepl(",", line, fixed = TRUE)) nta_csv_cells(line) else strsplit(line, "\t")[[1]]
    if (length(extracted) >= 2) return(trimws(extracted[2]))
    value_part <- sub(pattern, "", line, ignore.case = TRUE)
    value_part <- trimws(sub("^[,\t]+", "", value_part))
    if (nzchar(value_part)) value_part else NA_character_
  }

  metadata <- list(
    software_version = extract_metadata_value(raw_lines, "^Software Version"),
    experiment_name = extract_metadata_value(raw_lines, "^Experiment Name"),
    operator_name = extract_metadata_value(raw_lines, "^Operator Name"),
    time_captured = extract_metadata_value(raw_lines, "^Time Captured")
  )

  has_wells <- any(grepl("^Well ID,", raw_lines))
  has_autosampler <- any(grepl("^\\[Auto.?sampler", raw_lines, ignore.case = TRUE))
  if (!has_wells && !has_autosampler) {
    return(parse_nta_manual_summary(raw_lines, metadata))
  }

  main_error <- NULL
  distribution_error <- NULL

  tidy_data <- tryCatch({
    well_id_line_index <- which(stringr::str_detect(raw_lines, "^Well ID"))
    if (!length(well_id_line_index)) stop("Well ID missing.")
    well_id_line_index <- well_id_line_index[1]

    sample_name_line_index <- which(stringr::str_detect(raw_lines, "^Sample Name"))
    if (!length(sample_name_line_index)) stop("Sample Name missing.")
    sample_name_line_index <- sample_name_line_index[1]

    possible_starts <- c("^Concentration Average", "^Mean Average", "^Mode Average", "^FTLA Size", "^Particles per Frame")
    metrics_start_indices <- which(stringr::str_detect(raw_lines, paste(possible_starts, collapse = "|")))
    if (!length(metrics_start_indices)) stop("Metrics start missing.")
    metrics_start_index <- min(metrics_start_indices)

    last_metric_index <- which(stringr::str_detect(raw_lines, "^D90 Standard Error"))
    if (!length(last_metric_index)) {
      potential_metric_indices <- which(stringr::str_detect(raw_lines, "^[A-Za-z0-9() /%-].*,[ 0-9.,\\+eE-]+"))
      potential_metric_indices <- potential_metric_indices[potential_metric_indices > well_id_line_index]
      if (!length(potential_metric_indices)) stop("Fallback for final metric row failed.")
      last_metric_index <- max(potential_metric_indices)
    } else {
      last_metric_index <- last_metric_index[length(last_metric_index)]
    }

    graph_data_loc <- which(stringr::str_detect(raw_lines, stringr::fixed("[Auto-sampler Graph Data]")))
    size_data_loc <- which(stringr::str_detect(raw_lines, stringr::fixed("[Size Data]")))
    section_starts <- c(graph_data_loc, size_data_loc)
    if (length(section_starts)) {
      earliest_section_start <- min(section_starts)
      if (last_metric_index >= earliest_section_start) {
        last_metric_index <- earliest_section_start - 1
        if (last_metric_index < metrics_start_index) stop("Conflicting section headers.")
      }
    }

    read_start_line <- well_id_line_index - 1
    read_n_lines <- last_metric_index - well_id_line_index + 1
    if (read_n_lines <= 0) stop("Invalid main-summary block length.")

    tf_main <- tempfile(fileext = ".csv")
    on.exit(unlink(tf_main), add = TRUE)
    writeLines(raw_lines[(read_start_line + 1):(read_start_line + read_n_lines)], tf_main)
    data_block <- readr::read_csv(
      tf_main, col_names = FALSE, col_types = cols(.default = "c"),
      guess_max = read_n_lines + 5, progress = FALSE, show_col_types = FALSE
    )

    well_ids <- unlist(data_block[1, -1], use.names = FALSE)
    valid_cols <- which(!is.na(well_ids) & nzchar(trimws(well_ids)))
    well_ids <- trimws(well_ids[valid_cols])
    if (!length(well_ids)) stop("No valid Well IDs in main summary block.")

    sample_name_row <- sample_name_line_index - read_start_line
    if (sample_name_row < 1 || sample_name_row > nrow(data_block)) stop("Sample Name row is outside the parsed block.")
    sample_names <- unlist(data_block[sample_name_row, -1], use.names = FALSE)[valid_cols]
    sample_names <- trimws(sample_names)
    if (length(sample_names) != length(well_ids)) stop("Well ID/Sample Name mismatch.")

    metric_start_row <- metrics_start_index - read_start_line
    metric_end_row <- last_metric_index - read_start_line
    row_numbers <- seq_len(nrow(data_block))
    metric_rows <- which(
      !is.na(data_block$X1) & nzchar(data_block$X1) & !startsWith(data_block$X1, "[") &
        row_numbers >= metric_start_row & row_numbers <= metric_end_row
    )
    if (!length(metric_rows)) stop("No metric rows found in main summary block.")

    metric_names <- data_block$X1[metric_rows]
    metric_values <- data_block[metric_rows, valid_cols + 1, drop = FALSE]
    if (ncol(metric_values) != length(well_ids)) stop("Metric-column count does not match sample count.")

    metric_values_t <- as.data.frame(t(metric_values), stringsAsFactors = FALSE)
    colnames(metric_values_t) <- make.unique(metric_names)

    parsed <- dplyr::bind_cols(
      tibble::tibble(WellID = well_ids, SampleName = sample_names),
      metric_values_t
    ) %>%
      clean_names_df()

    metric_columns <- setdiff(names(parsed), c("well_id", "sample_name"))
    parsed <- parsed %>%
      dplyr::mutate(
        dplyr::across(
          dplyr::all_of(metric_columns),
          ~ readr::parse_number(as.character(.x), na = c("", "NA", "#DIV/0!"))
        ),
        time_captured = metadata$time_captured %||% NA_character_
      )

    parsed
  }, error = function(e) {
    main_error <<- conditionMessage(e)
    NULL
  })

  graph_section <- which(stringr::str_detect(raw_lines, stringr::fixed("[Auto-sampler Graph Data]")))
  tidy_distribution_data <- NULL

  if (length(graph_section)) {
    tidy_distribution_data <- tryCatch({
      graph_start <- graph_section[1]
      header_idx <- NA_integer_
      for (offset in 1:10) {
        idx <- graph_start + offset
        if (idx > length(raw_lines)) break
        if (grepl("^Bin centre \\(nm\\),", raw_lines[idx], ignore.case = TRUE)) {
          header_idx <- idx
          break
        }
      }
      if (is.na(header_idx)) stop("Could not find 'Bin centre (nm)' header.")

      well_line <- NA_integer_; sample_line <- NA_integer_
      search_up_limit <- max(1, header_idx - 8)
      for (idx in seq(header_idx - 1, search_up_limit, by = -1)) {
        if (is.na(sample_line) && str_starts_compat(raw_lines[idx], "Sample Name,")) sample_line <- idx
        if (is.na(well_line) && str_starts_compat(raw_lines[idx], "Well Position,")) well_line <- idx
        if (!is.na(sample_line) && !is.na(well_line)) break
      }
      if (is.na(well_line) || is.na(sample_line)) stop("Could not locate distribution Well Position/Sample Name rows.")

      well_cells <- trimws(strsplit(raw_lines[well_line], ",", fixed = TRUE)[[1]][-1])
      sample_cells <- trimws(strsplit(raw_lines[sample_line], ",", fixed = TRUE)[[1]][-1])
      header_positions <- which(nzchar(well_cells))
      header_positions <- header_positions[header_positions <= length(sample_cells)]
      if (!length(header_positions)) stop("No valid distribution sample headers.")
      well_ids <- well_cells[header_positions]
      sample_names <- sample_cells[header_positions]
      sample_names[is.na(sample_names) | !nzchar(sample_names)] <- "Unnamed sample"

      first_data_row <- header_idx + 1
      last_data_row <- NA_integer_
      current_row <- first_data_row
      while (current_row <= length(raw_lines)) {
        current_line <- trimws(raw_lines[current_row])
        first_element <- trimws(strsplit(current_line, ",", fixed = TRUE)[[1]][1])
        if (grepl("^[0-9]+(\\.[0-9]+)?$", first_element)) {
          last_data_row <- current_row
        } else if (nzchar(current_line) && !startsWith(current_line, ",")) {
          if (!is.na(last_data_row)) break
        } else if (!nzchar(current_line) && !is.na(last_data_row) && current_row > last_data_row + 2) {
          break
        }
        current_row <- current_row + 1
      }
      if (is.na(last_data_row)) stop("No numeric distribution rows found.")

      tf_dist <- tempfile(fileext = ".csv")
      on.exit(unlink(tf_dist), add = TRUE)
      writeLines(raw_lines[header_idx:last_data_row], tf_dist)
      raw_dist <- readr::read_csv(
        tf_dist, col_types = cols(.default = "c"), progress = FALSE,
        show_col_types = FALSE, name_repair = "minimal", na = c("", "NA", "#DIV/0!")
      )

      valid_header <- header_positions[(1 + header_positions + 1) <= ncol(raw_dist)]
      if (!length(valid_header)) stop("Distribution block has no complete concentration/error column pairs.")
      keep <- match(valid_header, header_positions)
      header_positions <- valid_header
      well_ids <- well_ids[keep]
      sample_names <- sample_names[keep]
      n_samples <- length(header_positions)

      bin_centres <- readr::parse_number(raw_dist[[1]], na = c("", "NA", "#DIV/0!"))
      purrr::map_dfr(seq_len(n_samples), function(i) {
        concentration_col <- 1 + header_positions[i]
        error_col <- concentration_col + 1
        tibble::tibble(
          WellID_dist = well_ids[i],
          SampleName_dist = sample_names[i],
          BinCentre_nm = bin_centres,
          ConcentrationAverage_dist = readr::parse_number(raw_dist[[concentration_col]], na = c("", "NA", "#DIV/0!")),
          StandardError_dist = readr::parse_number(raw_dist[[error_col]], na = c("", "NA", "#DIV/0!"))
        )
      }) %>%
        dplyr::filter(!is.na(BinCentre_nm))
    }, error = function(e) {
      distribution_error <<- conditionMessage(e)
      NULL
    })
  }

  list(
    data = tidy_data,
    metadata = metadata,
    distribution_data = tidy_distribution_data,
    main_error = main_error,
    distribution_error = distribution_error,
    status = if (!is.null(tidy_data) && nrow(tidy_data) > 0) "success" else "failed"
  )
}

# Titles
format_col_name <- function(colnames_vec){
  sapply(colnames_vec, function(colname){
    if (is.na(colname) || colname == "none") return(colname)
    name <- stringr::str_replace_all(colname, "_", " ") |> stringr::str_to_title()
    if (grepl("Size|D10|D50|D90|Mode|Mean|Median|Sd", name, ignore.case = TRUE)) name <- paste0(name, " (nm)")
    else if (grepl("Concentration", name, ignore.case = TRUE)) {
      if (grepl("Particles Per Frame", name, ignore.case = TRUE)) name <- paste0(name, " (particles/frame)")
      else name <- paste0(name, " (particles/mL)")
    }
    name
  }, USE.NAMES = FALSE)
}

# --- UI ---
# Shorten display labels only; preserve imported names and sample identity.
comparison_display_names <- function(x) {
  labels <- trimws(as.character(x))
  labels <- sub("[ _]+[0-9]{4}-[0-9]{2}-[0-9]{2}([ T_]+[0-9]{2}[-:][0-9]{2}[-:][0-9]{2})?$", "", labels, perl = TRUE)
  labels[is.na(labels) | !nzchar(labels)] <- "Unnamed sample"
  labels
}

group_name_input <- function(id, label, value) {
  div(class = "form-group shiny-input-container", style = "width:100%;",
      tags$label(class = "sr-only", `for` = id, label),
      tags$input(id = id, type = "text", class = "form-control", value = value, list = "groupNameSuggestions"))
}

preserve_metric <- function(current, choices, fallback) {
  if (length(current) == 1 && current %in% choices) current else fallback
}

APP_VERSION <- "2.02"

app_logo_ui <- if (file.exists(file.path("www", "logo.PNG"))) {
  tags$img(src = "logo.PNG", height = "46px", alt = "shinyNTApro logo", class = "app-logo")
} else {
  div(class = "app-logo-fallback", "NTA")
}

github_link_ui <- if (nzchar(GITHUB_REPO_URL)) {
  tags$a(
    class = "btn btn-default header-button",
    icon("github"), "GitHub",
    href = GITHUB_REPO_URL, target = "_blank", rel = "noopener noreferrer"
  )
} else NULL

workflow_section <- function(title, icon_name, ..., open = FALSE) {
  details_tag <- tags$details(
    class = "workflow-section",
    tags$summary(icon(icon_name), tags$span(title)),
    div(class = "workflow-content", ...)
  )
  if (isTRUE(open)) details_tag$attribs$open <- "open"
  details_tag
}

ui <- fluidPage(
  shinyjs::useShinyjs(),
  theme = shinythemes::shinytheme("flatly"),
  tags$head(
    tags$title(paste0("shinyNTApro UI/UX v", APP_VERSION)),
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("
      :root {
        --nta-primary: #2c6e7f;
        --nta-primary-dark: #1f5260;
        --nta-accent: #3d8fa3;
        --nta-bg: #f4f7f9;
        --nta-card: #ffffff;
        --nta-border: #dfe7eb;
        --nta-text-muted: #60717a;
        --nta-success: #2e7d32;
        --nta-warning: #a86400;
        --nta-danger: #b3261e;
      }
      body { background: var(--nta-bg); }
      .container-fluid { padding-left: 20px; padding-right: 20px; }
      .app-header {
        background: #fff; border-bottom: 1px solid var(--nta-border);
        margin: 0 -20px 18px; padding: 12px 20px;
        box-shadow: 0 2px 8px rgba(32, 61, 72, 0.06);
      }
      .app-brand { display: flex; align-items: center; gap: 12px; }
      .app-logo { max-width: 150px; object-fit: contain; }
      .app-logo-fallback {
        width: 46px; height: 46px; border-radius: 12px; display: flex;
        align-items: center; justify-content: center; color: #fff;
        background: var(--nta-primary); font-weight: 700;
      }
      .app-title { margin: 0; font-size: 24px; font-weight: 700; color: #26383f; }
      .app-subtitle { color: var(--nta-text-muted); margin-top: 2px; font-size: 13px; }
      .header-actions { text-align: right; padding-top: 3px; }
      .header-actions .btn { margin-left: 7px; margin-bottom: 4px; }
      .header-button, .header-actions .shiny-download-link { border-radius: 8px; }
      .sidebar-shell {
        background: #fff; border: 1px solid var(--nta-border); border-radius: 12px;
        padding: 12px; box-shadow: 0 3px 12px rgba(32, 61, 72, 0.05);
      }
      .workflow-section { border-bottom: 1px solid var(--nta-border); padding: 2px 0; }
      .workflow-section:last-of-type { border-bottom: 0; }
      .workflow-section > summary {
        cursor: pointer; list-style: none; padding: 12px 4px; font-weight: 650;
        color: #31454d; display: flex; align-items: center; gap: 8px;
      }
      .workflow-section > summary::-webkit-details-marker { display: none; }
      .workflow-section > summary:after { content: '\\25B8'; margin-left: auto; transition: transform .18s; }
      .workflow-section[open] > summary:after { transform: rotate(90deg); }
      .workflow-content { padding: 2px 4px 13px; }
      .workflow-help { color: var(--nta-text-muted); font-size: 12px; margin-bottom: 10px; }
      .sidebar-footer { margin-top: 12px; display: grid; gap: 8px; }
      .sidebar-footer .btn, .sidebar-footer .shiny-download-link { width: 100%; }
      .btn-primary { background-color: var(--nta-primary); border-color: var(--nta-primary); }
      .btn-primary:hover, .btn-primary:focus { background-color: var(--nta-primary-dark); border-color: var(--nta-primary-dark); }
      .btn-danger-outline { color: var(--nta-danger); background: #fff; border: 1px solid #e5b8b4; }
      .nta-card {
        background: var(--nta-card); border: 1px solid var(--nta-border); border-radius: 12px;
        padding: 18px; margin-bottom: 16px; box-shadow: 0 3px 12px rgba(32, 61, 72, 0.05);
      }
      .nta-card h3, .nta-card h4 { margin-top: 0; }
      .card-header-row { display: flex; justify-content: space-between; align-items: center; gap: 10px; margin-bottom: 12px; }
      .card-header-row h3, .card-header-row h4, .card-header-row h5 { margin: 0; }
      .card-actions .btn, .card-actions .shiny-download-link { margin-left: 6px; }
      .metric-grid { display: grid; grid-template-columns: repeat(4, minmax(150px, 1fr)); gap: 12px; margin-bottom: 16px; }
      .metric-card {
        background: #fff; border: 1px solid var(--nta-border); border-radius: 12px;
        padding: 16px; min-height: 108px; box-shadow: 0 2px 8px rgba(32, 61, 72, 0.04);
      }
      .metric-label { color: var(--nta-text-muted); font-size: 12px; text-transform: uppercase; letter-spacing: .04em; }
      .metric-value { font-size: 28px; font-weight: 700; margin: 7px 0 3px; color: #26383f; line-height: 1.05; }
      .metric-detail { color: var(--nta-text-muted); font-size: 12px; }
      .status-pass { color: var(--nta-success); }
      .status-fail { color: var(--nta-danger); }
      .status-warning { color: var(--nta-warning); }
      .status-neutral { color: #425a64; }
      .status-badge {
        display: inline-block; border-radius: 999px; padding: 4px 9px; font-size: 12px;
        font-weight: 650; background: #edf3f5; color: #425a64;
      }
      .status-badge.pass { background: #e7f4e8; color: var(--nta-success); }
      .status-badge.fail { background: #fdecea; color: var(--nta-danger); }
      .status-badge.warning { background: #fff3dc; color: var(--nta-warning); }
      .hero-card {
        background: linear-gradient(135deg, #ffffff 0%, #edf6f8 100%);
        padding: 28px; border-radius: 16px; border: 1px solid #d6e7eb; margin-bottom: 18px;
      }
      .hero-card h2 { margin-top: 0; font-weight: 700; }
      .hero-actions .btn { margin-right: 8px; margin-top: 8px; }
      .step-grid { display: grid; grid-template-columns: repeat(4, minmax(140px, 1fr)); gap: 12px; margin-top: 18px; }
      .step-card { background: rgba(255,255,255,.85); border: 1px solid #dce9ec; border-radius: 10px; padding: 14px; }
      .step-number { width: 25px; height: 25px; border-radius: 50%; background: var(--nta-primary); color: #fff; display: inline-flex; align-items: center; justify-content: center; font-weight: 700; margin-right: 6px; }
      .next-action { border-left: 4px solid var(--nta-accent); background: #f1f8fa; padding: 13px 15px; border-radius: 8px; }
      .plot-card .shiny-plot-output { min-height: 320px; }
      .control-strip { background: #f8fafb; border: 1px solid var(--nta-border); border-radius: 10px; padding: 12px 14px; margin-bottom: 14px; }
      .control-strip .form-group { margin-bottom: 8px; }
      .details-panel { border: 1px solid var(--nta-border); border-radius: 9px; padding: 0 13px; background: #fbfcfd; margin-top: 12px; }
      .details-panel > summary { cursor: pointer; padding: 12px 0; font-weight: 650; }
      .details-panel[open] { padding-bottom: 12px; }
      .empty-state { text-align: center; padding: 42px 18px; color: var(--nta-text-muted); }
      .empty-state .fa { font-size: 38px; margin-bottom: 12px; color: #94a9b2; }
      .data-count { color: var(--nta-text-muted); margin-top: 7px; }
      .tab-content { padding-top: 16px; }
      .nav-tabs > li > a { font-weight: 600; padding: 10px 9px; }
      .nav-tabs { display: flex; flex-wrap: wrap; }
      .nav-tabs > li.active > a, .nav-tabs > li.active > a:hover { color: var(--nta-primary-dark); border-top: 3px solid var(--nta-primary); }
      .progress-list { margin: 0; padding: 0; list-style: none; }
      .progress-list li { padding: 9px 0; border-bottom: 1px solid #edf1f3; }
      .progress-list li:last-child { border-bottom: 0; }
      .small-muted { color: var(--nta-text-muted); font-size: 12px; }
      .warning-banner { background: #fff8e8; border: 1px solid #efd69b; border-radius: 9px; padding: 11px 13px; color: #7d5200; margin-bottom: 12px; }
      .info-banner { background: #eff7fa; border: 1px solid #cbe0e6; border-radius: 9px; padding: 11px 13px; color: #355b67; margin-bottom: 12px; }
      .form-control, .selectize-input, .btn, .shiny-download-link { border-radius: 7px; }
      @media (max-width: 1100px) {
        .metric-grid { grid-template-columns: repeat(2, minmax(150px, 1fr)); }
        .step-grid { grid-template-columns: repeat(2, minmax(140px, 1fr)); }
      }
      @media (max-width: 767px) {
        .container-fluid { padding-left: 10px; padding-right: 10px; }
        .app-header { margin-left: -10px; margin-right: -10px; }
        .header-actions { text-align: left; margin-top: 10px; }
        .header-actions .btn { margin-left: 0; margin-right: 6px; }
        .metric-grid, .step-grid { grid-template-columns: 1fr; }
        .sidebar-shell { margin-bottom: 14px; }
        .card-header-row { align-items: flex-start; flex-direction: column; }
        .card-actions .btn, .card-actions .shiny-download-link { margin-left: 0; margin-right: 5px; margin-top: 5px; }
      }
    ")),
    tags$script(HTML("
      Shiny.addCustomMessageHandler('nta-click-upload', function(message) {
        $('#summaryFile input[type=file]').trigger('click');
      });
    "))
  ),

  div(
    class = "app-header",
    fluidRow(
      column(
        7,
        div(
          class = "app-brand",
          app_logo_ui,
          div(
            h1(class = "app-title", "shinyNTApro"),
            div(class = "app-subtitle", "NTA quality control, visualisation and statistical analysis"),
            div(class = "status-badge", paste0("UI/UX v", APP_VERSION))
          )
        )
      ),
      column(
        5,
        div(
          class = "header-actions",
          actionButton("helpButton", "Help", icon = icon("question-circle"), class = "btn-default header-button"),
          github_link_ui,
          downloadButton("generateReport", "Full report", class = "btn-primary header-button")
        )
      )
    )
  ),

  sidebarLayout(
    sidebarPanel(
      width = 3,
      div(
        class = "sidebar-shell",
        div(
          id = "sidebar_inputs",
          workflow_section(
            "Data upload", "upload", open = TRUE,
            fileInput(
              "summaryFile", "Upload NTA summary files",
              accept = c("text/csv", ".csv"),
              placeholder = "Select NTA summary CSV files...", multiple = TRUE
            ),
            div(class = "workflow-help", paste0("Up to ", MAX_UPLOAD_FILES, " files; maximum ", MAX_UPLOAD_MB_PER_FILE, " MB per file.")),
            div(class = "workflow-help", "Supports autosampler and manual-loading NS300 summaries. Manual-loading files use one sample per file, with exported averages and Well ID shown as Not applicable."),
            actionButton("loadExampleSummary", "Load example data", icon = icon("play-circle"), class = "btn-info btn-sm"),
            uiOutput("sidebarUploadStatus")
          ),
          workflow_section(
            "Quality control", "check-circle",
            fileInput("qcFile", "Upload QC summary (optional)", accept = c("text/csv", ".csv"), placeholder = "QC summary CSV..."),
            numericInput("qcTargetSize", "QC target size (nm)", DEFAULT_QC_TARGET_SIZE, min = 1),
            numericInput("qcTolerance", "QC tolerance (%)", DEFAULT_QC_TOLERANCE_PERCENT, 0, 100, 1),
            div(class = "workflow-help", "Review the target and tolerance, then save the current assessment to QC history."),
            actionButton("saveQcResult", "Save QC result", icon = icon("save")),
            uiOutput("qcSaveStatus")
          ),
          workflow_section(
            "Filters", "filter",
            div(class = "workflow-help", "Filters apply to sample plots, statistics, data exports and group editing. QC and QC History use their own controls."),
            selectizeInput("filterSourceFile", "Source file", NULL, multiple = TRUE),
            selectizeInput("filterSampleName", "Sample name", NULL, multiple = TRUE),
            selectizeInput("filterWellID", "Well ID", NULL, multiple = TRUE),
            actionButton("clearFilters", "Clear filters", icon = icon("eraser"), class = "btn-default btn-sm")
          ),
          workflow_section(
            "Statistical groups", "users",
            div(class = "workflow-help", "Assign groups here, or upload a completed group CSV."),
            actionButton("editGroups", "Assign groups", icon = icon("edit")),
            tags$br(), tags$br(),
            div(class = "workflow-help", "Alternatively, use a CSV template."),
            downloadButton("downloadGroupTemplate", "Download group template", class = "btn-default"),
            tags$br(), tags$br(),
            fileInput("groupFile", "Upload completed group CSV", accept = c("text/csv", ".csv"), placeholder = "Completed group CSV..."),
            uiOutput("sidebarGroupStatus")
          ),
          workflow_section(
            "Advanced settings", "cog",
            textInput("qcLogDirPath", "QC log directory path", "."),
            div(class = "workflow-help", "For local use, '.' writes to the app directory. Hosted users can re-upload their QC log each session."),
            checkboxInput("showTechnicalColumns", "Show technical identity columns in the data table", value = FALSE)
          )
        ),
        div(
          class = "sidebar-footer",
          actionButton("resetAll", "Reset analysis", icon = icon("refresh"), class = "btn-danger-outline")
        )
      )
    ),

    mainPanel(
      width = 9,
      tabsetPanel(
        type = "tabs", id = "mainTabs", selected = "overviewTab",

        tabPanel(
          title = tagList(icon("home"), "Overview"), value = "overviewTab",
          uiOutput("overviewContent")
        ),

        tabPanel(
          title = tagList(icon("check-circle"), "QC"), value = "qcTab",
          div(
            class = "nta-card",
            div(class = "card-header-row", h3("Quality control")),
            uiOutput("qcStatusDisplay")
          ),
          div(
            class = "nta-card",
            div(class = "card-header-row", h4("Individual QC results")),
            DT::dataTableOutput("qcDetailsTable")
          )
        ),

        tabPanel(
          title = tagList(icon("bar-chart"), "Comparison"), value = "plotTab",
          div(
            class = "nta-card",
            div(class = "card-header-row", h3("Comparison plots")),
            div(
              class = "control-strip",
              fluidRow(
                column(3, selectInput("comparisonPlotType", "Display", c("Points" = "points", "Box plot" = "box", "Violin" = "violin"))),
                column(3, selectInput("comparisonXAxis", "X-axis", c("Sample" = "sample_key", "Source file" = "source_file"))),
                column(3, selectInput("comparisonColorBy", "Colour by", c("Sample" = "sample_key", "Source file" = "source_file", "None" = "none"))),
                column(3, radioButtons("comparisonLayout", "Layout", c("Side by side" = "side", "Stacked" = "stack"), inline = TRUE))
              ),
              fluidRow(
                column(3, selectInput("plotSize", "Plot size", c("Compact" = "compact", "Standard" = "standard", "Large" = "large"), selected = "standard")),
                column(3, selectInput("comparisonSort", "Sample order", c("Original" = "original", "Alphabetical" = "alphabetical", "By value" = "value"))),
                column(3, checkboxInput("comparisonLogY", "Use Log10 Y-axis", value = FALSE)),
                column(3, checkboxInput("comparisonShowLabels", "Show X-axis labels", value = TRUE))
              )
            ),
            fluidRow(
              column(6, selectInput("metricSelect1", "Plot 1 metric", NULL)),
              column(6, selectInput("metricSelect2", "Plot 2 metric", NULL))
            )
          ),
          actionButton("showSampleNames", "View full sample names", icon = icon("list")),
          div(class = "small-muted", "Sample names are shown without trailing dates and times. View full sample names for the original names and source files."),
          uiOutput("comparisonPlotsUI")
        ),

        tabPanel(
          title = tagList(icon("line-chart"), "Statistics"), value = "statsTab",
          div(
            class = "nta-card",
            div(class = "card-header-row", h3("Statistical comparison")),
            uiOutput("statsPrerequisiteMessage"),
            fluidRow(
              column(6, selectInput("statsMetric", "Particle parameter", NULL, width = "100%")),
              column(6, selectInput("statsTestType", "Test selection", statistical_test_choices(),
                                    selected = "Auto", selectize = FALSE, width = "100%"))
            ),
            selectInput("statsDesign", "Sample relationship", c(
              "Independent samples in each group" = "independent",
              "Paired / repeated measurements (not supported)" = "paired"
            ), selected = "independent", width = "100%"),
            uiOutput("statsSelectionHelp"),
            tags$details(
              tags$summary("Sample independence and multiple parameters"),
              div(class = "small-muted", "Each row must represent an independent sample; repeat recordings or aliquots of the same sample do not increase n. Paired or repeated-measures analyses are not supported. Each parameter is tested separately; p-values are not adjusted across different parameters.")
            )
          ),
          uiOutput("statsSummaryCards"),
          div(
            class = "nta-card plot-card",
            div(
              class = "card-header-row",
              h4("Group comparison plot"),
              div(class = "card-actions", downloadButton("downloadStatsPlot", "Download PNG", class = "btn-default btn-sm"))
            ),
            checkboxInput("statsShowNS", "Show non-significant comparisons (NS)", value = FALSE),
            plotOutput("statsPlot", height = 560)
          ),
          div(
            class = "nta-card",
            div(class = "card-header-row", h4("Group summary")),
            DT::dataTableOutput("statsGroupSummary")
          ),
          tags$details(
            class = "details-panel",
            tags$summary("Show full statistical details"),
            downloadButton("downloadStatsResults", "Download results (.txt)", class = "btn-default btn-sm"),
            verbatimTextOutput("statsOutput")
          )
        ),

        tabPanel(
          title = tagList(icon("area-chart"), "Distributions"), value = "distPlotTab",
          div(
            class = "nta-card",
            div(class = "card-header-row", h3("Particle size distributions")),
            div(
              class = "control-strip",
              fluidRow(
                column(6, selectizeInput(
                  "distSampleSelect", "Samples", NULL, multiple = TRUE, width = "100%",
                  options = list(
                    maxItems = 5000,
                    maxOptions = 5000,
                    closeAfterSelect = FALSE,
                    plugins = list("remove_button"),
                    placeholder = "Click or type to add samples..."
                  )
                )),
                column(3, selectizeInput("distSourceSelect", "Select by source file", NULL, multiple = TRUE)),
                column(3, tags$div(style = "padding-top:25px;", actionButton("distApplySource", "Apply source", class = "btn-default btn-sm")))
              ),
              fluidRow(
                column(3, actionButton("distSelectAll", "Select all", class = "btn-default btn-sm"), actionButton("distClearAll", "Clear", class = "btn-default btn-sm")),
                column(3, checkboxInput("distShowError", "Show standard error", value = TRUE)),
                column(3, checkboxInput("distLogY", "Log10 concentration", value = FALSE)),
                column(3, checkboxInput("distAverageByGroup", "Average by statistical group", value = FALSE))
              ),
              fluidRow(
                column(6, selectInput("distHighlight", "Highlight one curve", choices = c("None" = ""))),
                column(6, selectInput("distPlotSize", "Plot size", c("Compact" = "compact", "Standard" = "standard", "Large" = "large"), selected = "standard"))
              ),
              sliderInput("distBinRange", "Bin centre range (nm)", min = 0, max = 1000, value = c(0, 1000), step = 1, width = "100%"),
              uiOutput("distSelectionStatus"),
              div(class = "small-muted", "Three curves are selected initially for readability. Click the Samples field, use Select all, or add samples by source file to display more.")
            )
          ),
          div(
            class = "nta-card plot-card",
            div(
              class = "card-header-row",
              h4("Distribution plot"),
              div(class = "card-actions", downloadButton("downloadDistPlot", "Download PNG", class = "btn-default btn-sm"))
            ),
            uiOutput("distributionPlotContainer")
          )
        ),

        tabPanel(
          title = tagList(icon("history"), "QC History"), value = "qcLogAnalysisTab",
          div(
            class = "nta-card",
            div(class = "card-header-row", h3("QC history and precision")),
            fluidRow(
              column(4, fileInput("qcLogUpload", "Upload QC log", accept = c("text/csv", ".csv"), placeholder = "QC log CSV...")),
              column(4, selectInput("qcAnalysisMetric", "Metric", c("Mean size" = "MeanAverage", "Mode size" = "ModeAverage", "Concentration" = "ConcentrationAverage"))),
              column(4, dateRangeInput("qcDateRange", "Date range", start = Sys.Date() - 30, end = Sys.Date()))
            ),
            div(class = "small-muted", paste0("Required columns: TimeCaptured, Status, and at least one supported metric (MeanAverage, ModeAverage or ConcentrationAverage). Available metrics and the date range are detected automatically. Expected file: ", QC_LOG_FILE_BASENAME, ".")),
            tags$details(
              class = "details-panel",
              tags$summary("Advanced QC limits and baseline"),
              fluidRow(
                column(6, numericInput("intraDayCVLimit", "Maximum intra-day CV (%)", DEFAULT_INTRA_CV_LIMIT, 0, 100, 1)),
                column(6, numericInput("interDayCVLimit", "Maximum inter-day CV (%)", DEFAULT_INTER_CV_LIMIT, 0, 100, 1))
              ),
              checkboxInput("useFixedQcBaseline", "Use a fixed baseline mean and SD", value = FALSE),
              conditionalPanel(
                condition = "input.useFixedQcBaseline",
                fluidRow(
                  column(6, numericInput("qcBaselineMean", "Baseline mean", value = 100)),
                  column(6, numericInput("qcBaselineSD", "Baseline SD", value = 5, min = 0))
                )
              )
            )
          ),
          uiOutput("qcHistoryCards"),
          div(
            class = "nta-card plot-card",
            div(class = "card-header-row", h4("Levey-Jennings chart"),
              div(class = "card-actions",
                downloadButton("downloadQcHistoryPNG", "Download PNG", class = "btn-default btn-sm"),
                downloadButton("downloadQcHistoryPDF", "Download PDF", class = "btn-default btn-sm")
              )
            ),
            plotOutput("leveyJenningsPlot", height = 540)
          ),
          div(
            class = "nta-card",
            div(class = "card-header-row", h4("Daily summaries")),
            DT::dataTableOutput("qcDailyTable")
          ),
          tags$details(class = "details-panel", tags$summary("Show precision calculation details"), verbatimTextOutput("qcAnalysisSummary"))
        ),

        tabPanel(
          title = tagList(icon("table"), "Data"), value = "dataTab",
          div(
            class = "nta-card",
            div(
              class = "card-header-row",
              div(h3("Parsed data"), uiOutput("dataRowCount")),
              div(class = "card-actions", downloadButton("downloadData", "Download filtered CSV", class = "btn-default btn-sm"))
            ),
            selectizeInput("visibleDataColumns", "Visible columns", choices = NULL, multiple = TRUE, width = "100%"),
            DT::dataTableOutput("parsedTable")
          )
        ),

        tabPanel(
          title = tagList(icon("wrench"), "Diagnostics"), value = "parsingInfoTab",
          div(
            class = "nta-card",
            div(class = "card-header-row", h3("File parsing status")),
            DT::dataTableOutput("parsingStatusTable")
          ),
          div(
            class = "nta-card",
            div(class = "card-header-row", h4("File metadata")),
            DT::dataTableOutput("fileMetadataDisplay")
          ),
          tags$details(
            class = "details-panel",
            tags$summary("Show technical parsing log"),
            verbatimTextOutput("parsingInfo")
          )
        )
      )
    )
  )
)

# --- Server ---

server <- function(input, output, session) {
  
  rv <- reactiveValues(
    parsed_result = NULL,
    grouped_data = NULL,
    group_map = NULL,
    qc_saved_signature = NULL,
    qc_saved_path = NULL,
    group_editor_samples = NULL,
    active_files_info = NULL,
    processed_files_info = NULL,
    parse_log = character(),
    active_qc_file_info = list(path = NULL, name = NULL)
  )

  observe({
    has_data <- !is.null(rv$parsed_result) && !is.null(rv$parsed_result$data) && nrow(rv$parsed_result$data) > 0
    shinyjs::toggleState("generateReport", condition = has_data)
    shinyjs::toggleState("downloadGroupTemplate", condition = has_data)
  })

  process_summary_files <- function(files_df) {
    req(files_df)
    if (nrow(files_df) > MAX_UPLOAD_FILES) {
      showNotification(paste0("Upload no more than ", MAX_UPLOAD_FILES, " files at once."), type = "error")
      return(NULL)
    }
    if (any(is.na(files_df$size)) || any(files_df$size > MAX_UPLOAD_MB_PER_FILE * 1024^2)) {
      showNotification(paste0("Each file must be ", MAX_UPLOAD_MB_PER_FILE, " MB or smaller."), type = "error")
      return(NULL)
    }

    files_df <- files_df %>%
      dplyr::mutate(
        name = basename(name),
        source_file = make.unique(name, sep = "_")
      )
    file_paths <- stats::setNames(files_df$datapath, files_df$source_file)

    withProgress(message = "Parsing files...", value = 0, {
      all_results <- purrr::map2(file_paths, names(file_paths), function(path, source_id) {
        incProgress(1 / length(file_paths), detail = source_id)
        result <- tryCatch(
          parse_nta_summary(path),
          error = function(e) list(
            data = NULL, metadata = list(), distribution_data = NULL,
            main_error = conditionMessage(e), distribution_error = NULL, status = "failed"
          )
        )
        result$source_file <- source_id
        result
      })
      names(all_results) <- names(file_paths)

      successful_results <- purrr::keep(all_results, ~ identical(.x$status, "success") && !is.null(.x$data) && nrow(.x$data) > 0)
      failed_results <- purrr::keep(all_results, ~ !identical(.x$status, "success") || is.null(.x$data) || nrow(.x$data) == 0)

      issue_messages <- character()
      if (length(failed_results)) {
        issue_messages <- c(issue_messages, purrr::imap_chr(failed_results, function(x, nm) {
          paste0(nm, ": main-data parsing failed - ", x$main_error %||% "unknown error")
        }))
      }
      distribution_warnings <- purrr::imap_chr(successful_results, function(x, nm) {
        if (!is.null(x$distribution_error)) paste0(nm, ": distribution parsing failed - ", x$distribution_error) else NA_character_
      })
      issue_messages <- c(issue_messages, stats::na.omit(distribution_warnings))

      if (length(issue_messages)) {
        safe_issues <- htmltools::htmlEscape(issue_messages)
        showNotification(
          HTML(paste("Parsing issues detected:<br>", paste(safe_issues, collapse = "<br>"))),
          type = "warning", duration = 15
        )
      }

      if (!length(successful_results)) {
        rv$parsed_result <- NULL
        rv$grouped_data <- NULL
        rv$processed_files_info <- files_df
        rv$parse_log <- c("No files were successfully parsed.", issue_messages)
        showNotification("No summary files could be successfully parsed.", type = "error", duration = 10)
        return(NULL)
      }

      combined_data <- purrr::imap_dfr(successful_results, function(x, source_id) {
        x$data %>%
          dplyr::mutate(
            source_file = source_id,
            sample_key = make_sample_key(source_file, sample_name, well_id),
            .before = 1
          )
      })

      combined_dist_data <- purrr::imap_dfr(successful_results, function(x, source_id) {
        if (is.null(x$distribution_data) || !nrow(x$distribution_data)) return(tibble::tibble())
        x$distribution_data %>%
          dplyr::mutate(
            source_file = source_id,
            sample_key = make_sample_key(source_file, SampleName_dist, WellID_dist),
            .before = 1
          )
      })

      combined_metadata <- purrr::imap_dfr(successful_results, function(x, source_id) {
        tibble::as_tibble(x$metadata) %>%
          dplyr::mutate(
            source_file = source_id,
            main_status = "Parsed",
            distribution_status = if (is.null(x$distribution_error)) {
              if (!is.null(x$distribution_data) && nrow(x$distribution_data)) "Parsed" else "Not present"
            } else {
              paste0("Failed: ", x$distribution_error)
            },
            .before = 1
          )
      })

      rv$parsed_result <- list(
        data = combined_data,
        distribution_data = if (nrow(combined_dist_data)) combined_dist_data else NULL,
        metadata = combined_metadata
      )
      rv$grouped_data <- NULL
      rv$processed_files_info <- files_df
      rv$parse_log <- c(
        paste0("Successfully parsed: ", paste(names(successful_results), collapse = ", ")),
        if (length(issue_messages)) c("", "Issues:", issue_messages) else NULL
      )

      showNotification(paste(length(successful_results), "file(s) parsed and combined."), type = "message")

      updateSelectizeInput(session, "filterSourceFile", choices = unique(combined_data$source_file), selected = character(), server = TRUE)
      updateSelectizeInput(session, "filterSampleName", choices = sort(unique(stats::na.omit(combined_data$sample_name))), selected = character(), server = TRUE)
      updateSelectizeInput(session, "filterWellID", choices = sort(unique(stats::na.omit(combined_data$well_id))), selected = character(), server = TRUE)

      if (!is.null(rv$parsed_result$distribution_data)) {
        dist_data_all <- rv$parsed_result$distribution_data
        dist_keys <- unique(dist_data_all$sample_key)
        dist_choices <- stats::setNames(dist_keys, dist_keys)
        updateSelectizeInput(
          session, "distSampleSelect", choices = dist_choices,
          selected = head(dist_keys, 3), server = FALSE
        )
        updateSelectizeInput(
          session, "distSourceSelect", choices = sort(unique(dist_data_all$source_file)),
          selected = character(), server = TRUE
        )
        updateSelectInput(
          session, "distHighlight",
          choices = c("None" = "", stats::setNames(unique(dist_data_all$sample_key), unique(dist_data_all$sample_key))),
          selected = ""
        )
        min_bin <- floor(min(dist_data_all$BinCentre_nm, na.rm = TRUE))
        max_bin <- ceiling(max(dist_data_all$BinCentre_nm, na.rm = TRUE))
        if (is.finite(min_bin) && is.finite(max_bin) && min_bin < max_bin) {
          updateSliderInput(session, "distBinRange", min = min_bin, max = max_bin, value = c(min_bin, max_bin))
        }
      } else {
        updateSelectizeInput(session, "distSampleSelect", choices = character(), selected = character(), server = FALSE)
        updateSelectizeInput(session, "distSourceSelect", choices = character(), selected = character(), server = TRUE)
        updateSelectInput(session, "distHighlight", choices = c("None" = ""), selected = "")
      }
      updateTabsetPanel(session, "mainTabs", selected = "overviewTab")
    })
  }

  # MODIFIED: observeEvent for multiple files
  observeEvent(input$summaryFile, {
    req(input$summaryFile$datapath)
    shinyjs::reset("groupFile")
    rv$group_map <- NULL
    rv$group_editor_samples <- NULL
    rv$grouped_data <- NULL
    rv$active_files_info <- input$summaryFile
    if (!is.null(rv$active_qc_file_info$path)) rv$active_qc_file_info <- list(path = NULL, name = NULL)
  })
  
  load_example_data <- function() {
    example_summary_path <- EXAMPLE_SUMMARY_FILENAME
    if (!file.exists(example_summary_path)) {
      showNotification(paste("Example file not found:", example_summary_path), type = "error")
      return(invisible(NULL))
    }
    shinyjs::reset("groupFile")
    rv$group_map <- NULL
    rv$group_editor_samples <- NULL
    rv$grouped_data <- NULL
    rv$active_files_info <- tibble::tibble(
      name = basename(example_summary_path), size = file.size(example_summary_path),
      type = "text/csv", datapath = example_summary_path
    )
    shinyjs::reset("summaryFile")
    invisible(NULL)
  }

  observeEvent(input$loadExampleSummary, load_example_data())
  observeEvent(input$landingLoadExample, load_example_data())
  observeEvent(input$landingUpload, {
    session$sendCustomMessage("nta-click-upload", list())
  })
  
  observeEvent(rv$active_files_info, {
    if (!is.null(rv$active_files_info)) {
      process_summary_files(rv$active_files_info)
    } else {
      rv$parsed_result <- NULL
      updateSelectizeInput(session, "filterSourceFile", choices = character())
      updateSelectizeInput(session, "filterSampleName", choices = character())
      updateSelectizeInput(session, "filterWellID", choices = character())
      updateSelectizeInput(session, "distSampleSelect", choices = character())
      updateSelectizeInput(session, "distSourceSelect", choices = character())
      updateTabsetPanel(session, "mainTabs", selected = "overviewTab")
    }
  }, ignoreNULL = FALSE, priority = 10)
  
  observeEvent(input$clearFilters, {
    updateSelectizeInput(session, "filterSourceFile", selected = character())
    updateSelectizeInput(session, "filterSampleName", selected = character())
    updateSelectizeInput(session, "filterWellID", selected = character())
  })

  observeEvent(input$helpButton, {
    showModal(modalDialog(
      title = "Using shinyNTApro",
      tags$ol(
        tags$li(tags$strong("Upload data:"), " select one or more NTA summary CSV files."),
        tags$li(tags$strong("Review:"), " confirm parsing, QC and available measurements in Overview."),
        tags$li(tags$strong("Assign groups:"), " click Assign groups to edit them directly, or complete and upload the group CSV template."),
        tags$li(tags$strong("Analyse and export:"), " use the result-specific download buttons or generate the full report.")
      ),
      tags$p(class = "small-muted", paste("Application version", APP_VERSION)),
      easyClose = TRUE,
      footer = modalButton("Close")
    ))
  })

  observeEvent(input$resetAll, {
    showModal(modalDialog(
      title = "Reset the analysis?",
      "All uploaded data, filters, group assignments and current selections will be cleared.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirmReset", "Reset everything", class = "btn-danger")
      )
    ))
  })

  observeEvent(input$confirmReset, {
    removeModal()
    updateSelectInput(session, "statsDesign", selected = "independent")
    updateSelectInput(session, "statsTestType", selected = "Auto")
    updateCheckboxInput(session, "statsShowNS", value = FALSE)
    shinyjs::reset("sidebar_inputs")
    rv$qc_saved_signature <- NULL
    rv$qc_saved_path <- NULL
    rv$group_map <- NULL
    rv$group_editor_samples <- NULL
    rv$parsed_result <- NULL
    rv$grouped_data <- NULL
    rv$active_files_info <- NULL
    rv$processed_files_info <- NULL
    rv$parse_log <- character()
    rv$active_qc_file_info <- list(path = NULL, name = NULL)
    updateSelectizeInput(session, "filterSourceFile", choices = character(), selected = character())
    updateSelectizeInput(session, "filterSampleName", choices = character(), selected = character())
    updateSelectizeInput(session, "filterWellID", choices = character(), selected = character())
    updateSelectizeInput(session, "distSampleSelect", choices = character(), selected = character())
    updateSelectizeInput(session, "distSourceSelect", choices = character(), selected = character())
    updateSelectInput(session, "distHighlight", choices = c("None" = ""), selected = "")
    shinyjs::reset("qcLogUpload")
    updateSelectInput(session, "comparisonPlotType", selected = "points")
    updateSelectInput(session, "comparisonXAxis", selected = "sample_key")
    updateSelectInput(session, "comparisonColorBy", selected = "sample_key")
    updateRadioButtons(session, "comparisonLayout", selected = "side")
    updateSelectInput(session, "plotSize", selected = "standard")
    updateSelectInput(session, "comparisonSort", selected = "original")
    updateCheckboxInput(session, "comparisonLogY", value = FALSE)
    updateCheckboxInput(session, "comparisonShowLabels", value = TRUE)
    updateCheckboxInput(session, "distShowError", value = TRUE)
    updateCheckboxInput(session, "distLogY", value = FALSE)
    updateCheckboxInput(session, "distAverageByGroup", value = FALSE)
    updateSelectInput(session, "distPlotSize", selected = "standard")
    updateTabsetPanel(session, "mainTabs", selected = "overviewTab")
    showNotification("Analysis reset.", type = "message")
  })
  
  parsed_data_raw <- reactive({ if (!is.null(rv$parsed_result) && !is.null(rv$parsed_result$data)) rv$parsed_result$data else NULL })
  parsed_metadata <- reactive({ req(rv$parsed_result, rv$parsed_result$metadata); rv$parsed_result$metadata })
  parsed_distribution_data <- reactive({ if (!is.null(rv$parsed_result) && !is.null(rv$parsed_result$distribution_data)) rv$parsed_result$distribution_data else NULL })
  
  # MODIFIED: Display metadata as an interactive table
  output$fileMetadataDisplay <- DT::renderDataTable({
    meta_df <- parsed_metadata(); req(meta_df)
    datatable(
      meta_df,
      options = list(scrollX = TRUE, dom = "t", pageLength = -1),
      rownames = FALSE
    )
  })
  
  filtered_data <- reactive({
    df <- parsed_data_raw()
    req(df)
    df_filtered <- df
    if (length(input$filterSourceFile)) df_filtered <- df_filtered %>% dplyr::filter(source_file %in% input$filterSourceFile)
    if (length(input$filterSampleName)) df_filtered <- df_filtered %>% dplyr::filter(sample_name %in% input$filterSampleName)
    if (length(input$filterWellID)) df_filtered <- df_filtered %>% dplyr::filter(well_id %in% input$filterWellID)
    validate(need(nrow(df_filtered) > 0, "No data remain after filtering."))
    df_filtered
  })

  metric_card <- function(label, value, detail = NULL, status_class = "status-neutral") {
    div(
      class = "metric-card",
      div(class = "metric-label", label),
      div(class = paste("metric-value", status_class), value),
      if (!is.null(detail)) div(class = "metric-detail", detail)
    )
  }

  output$sidebarUploadStatus <- renderUI({
    if (is.null(rv$parsed_result) || is.null(rv$parsed_result$data)) return(NULL)
    n_files <- dplyr::n_distinct(rv$parsed_result$data$source_file)
    n_samples <- dplyr::n_distinct(rv$parsed_result$data$sample_key)
    div(
      style = "margin-top:10px;",
      tags$span(class = "status-badge pass", icon("check"), paste(n_files, "file(s) parsed")),
      div(class = "small-muted", paste(n_samples, "sample records available"), style = "margin-top:5px;")
    )
  })

  output$sidebarGroupStatus <- renderUI({
    if (is.null(rv$grouped_data) || !nrow(rv$grouped_data)) {
      return(div(style = "margin-top:8px;", tags$span(class = "status-badge", "Groups not assigned")))
    }
    div(
      style = "margin-top:8px;",
      tags$span(class = "status-badge pass", icon("check"), paste(dplyr::n_distinct(rv$grouped_data$group), "groups applied")),
      div(class = "small-muted", paste(dplyr::n_distinct(rv$grouped_data$sample_key), "samples assigned"), style = "margin-top:5px;")
    )
  })

  output$overviewContent <- renderUI({
    if (is.null(rv$parsed_result) || is.null(rv$parsed_result$data)) {
      return(tagList(
        div(
          class = "hero-card",
          h2("Analyse NTA summary files with confidence"),
          p("Upload one or more NTA summary CSV files to review measurements, assess QC, compare samples, run group statistics and export publication-ready results."),
          div(
            class = "hero-actions",
            actionButton("landingUpload", "Upload NTA files", icon = icon("upload"), class = "btn-primary btn-lg"),
            actionButton("landingLoadExample", "Explore example data", icon = icon("play-circle"), class = "btn-default btn-lg")
          ),
          div(
            class = "step-grid",
            div(class = "step-card", tags$span(class = "step-number", "1"), tags$strong("Upload"), p(class = "small-muted", "Select up to 20 NTA summary files.")),
            div(class = "step-card", tags$span(class = "step-number", "2"), tags$strong("Review"), p(class = "small-muted", "Confirm parsing, QC and available metrics.")),
            div(class = "step-card", tags$span(class = "step-number", "3"), tags$strong("Group"), p(class = "small-muted", "Assign experimental groups for statistics.")),
            div(class = "step-card", tags$span(class = "step-number", "4"), tags$strong("Export"), p(class = "small-muted", "Download data, plots and a full report."))
          )
        ),
        div(
          class = "nta-card",
          h4("Before uploading"),
          p("Use CSV summary exports produced by the supported NTA software format. Distribution plots are available when the file contains an [Auto-sampler Graph Data] section."),
          div(class = "info-banner", icon("lock"), " Uploaded files are processed within the active app session. Avoid identifiable or confidential sample names when using a public deployment."),
          div(class = "small-muted", paste("shinyNTApro version", APP_VERSION))
        )
      ))
    }

    df_all <- rv$parsed_result$data
    df_filtered <- tryCatch(filtered_data(), error = function(e) df_all)
    n_files <- dplyr::n_distinct(df_all$source_file)
    n_samples <- dplyr::n_distinct(df_all$sample_key)
    n_filtered <- dplyr::n_distinct(df_filtered$sample_key)
    dist_df <- rv$parsed_result$distribution_data
    n_dist <- if (is.null(dist_df)) 0L else dplyr::n_distinct(dist_df$sample_key)
    group_count <- if (is.null(rv$grouped_data)) 0L else dplyr::n_distinct(rv$grouped_data$group)
    assigned_count <- if (is.null(rv$grouped_data)) 0L else dplyr::n_distinct(rv$grouped_data$sample_key)

    qc_status <- "Not uploaded"
    qc_class <- "status-neutral"
    if (!is.null(input$qcFile$datapath)) {
      qc_obj <- tryCatch(qc_results(), error = function(e) NULL)
      if (!is.null(qc_obj)) {
        qc_status <- qc_obj$status %||% "ERROR"
        qc_class <- switch(qc_status, PASS = "status-pass", FAIL = "status-fail", ERROR = "status-warning", "status-neutral")
      }
    }

    issues <- rv$parse_log[grepl("failed|issue|missing|not present", rv$parse_log, ignore.case = TRUE)]
    next_action <- if (is.null(rv$grouped_data)) {
      "Compare samples now, or click Assign groups to prepare statistical comparisons. QC is optional."
    } else {
      "Groups are assigned. Open Statistics to review usable sample counts and analysis requirements for your selected parameter."
    }

    file_items <- lapply(seq_len(nrow(rv$parsed_result$metadata)), function(i) {
      row <- rv$parsed_result$metadata[i, , drop = FALSE]
      dist_status <- as.character(row$distribution_status %||% "Not present")
      dist_badge <- if (grepl("^Parsed", dist_status)) "pass" else if (grepl("^Failed", dist_status)) "fail" else "warning"
      tags$li(
        tags$strong(as.character(row$source_file)),
        tags$span(class = "status-badge pass", "Summary parsed"),
        tags$span(class = paste("status-badge", dist_badge), paste("Distribution:", dist_status))
      )
    })

    tagList(
      div(
        class = "metric-grid",
        metric_card("Files loaded", n_files, paste(n_samples, "total samples")),
        metric_card("Active samples", n_filtered, if (n_filtered < n_samples) paste("Filtered from", n_samples) else "No filters applied"),
        metric_card("Distribution data", if (n_dist > 0) "Available" else "Unavailable", paste(n_dist, "samples with curves"), if (n_dist > 0) "status-pass" else "status-warning"),
        metric_card("QC status", qc_status, if (is.null(input$qcFile$datapath)) "Optional QC file" else basename(input$qcFile$name), qc_class)
      ),
      div(
        class = "nta-card",
        div(class = "card-header-row", h3("Analysis overview"), div(class = "card-actions", downloadButton("downloadDataOverview", "Download filtered data", class = "btn-default btn-sm"))),
        fluidRow(
          column(
            7,
            h4("Files"),
            do.call(tags$ul, c(list(class = "progress-list"), file_items))
          ),
          column(
            5,
            h4("Workflow status"),
            tags$ul(
              class = "progress-list",
              tags$li(icon("check", class = "status-pass"), paste(n_samples, "samples parsed")),
              tags$li(if (n_dist > 0) icon("check", class = "status-pass") else icon("exclamation-triangle", class = "status-warning"), paste(n_dist, "samples with distribution data")),
              tags$li(if (group_count > 0) icon("check", class = "status-pass") else icon("circle", class = "status-neutral"), if (group_count > 0) paste(group_count, "groups covering", assigned_count, "samples") else "Statistical groups not assigned")
            )
          )
        ),
        if (length(issues)) div(class = "warning-banner", icon("exclamation-triangle"), paste(unique(issues), collapse = " ")),
        div(class = "next-action", tags$strong("Suggested next action: "), next_action),
        div(
          class = "hero-actions", style = "margin-top:12px;",
          actionButton("goComparison", "Open comparison plots", icon = icon("bar-chart"), class = "btn-default"),
          actionButton("goStatistics", "Open statistics", icon = icon("line-chart"), class = "btn-default"),
          actionButton("goDiagnostics", "View diagnostics", icon = icon("wrench"), class = "btn-default")
        )
      )
    )
  })

  observeEvent(input$goComparison, updateTabsetPanel(session, "mainTabs", selected = "plotTab"))
  observeEvent(input$goStatistics, updateTabsetPanel(session, "mainTabs", selected = "statsTab"))
  observeEvent(input$goDiagnostics, updateTabsetPanel(session, "mainTabs", selected = "parsingInfoTab"))

  output$dataRowCount <- renderUI({
    df <- tryCatch(filtered_data(), error = function(e) NULL)
    all_df <- parsed_data_raw()
    if (is.null(df) || is.null(all_df)) return(NULL)
    div(class = "data-count", paste("Showing", nrow(df), "of", nrow(all_df), "rows"))
  })

  output$parsingStatusTable <- DT::renderDataTable({
    meta <- parsed_metadata()
    req(meta)
    display <- meta %>%
      dplyr::transmute(
        `Source file` = source_file,
        `Summary data` = paste0("<span class='status-badge pass'>", htmltools::htmlEscape(main_status), "</span>"),
        `Distribution data` = dplyr::case_when(
          grepl("^Parsed", distribution_status) ~ paste0("<span class='status-badge pass'>", htmltools::htmlEscape(distribution_status), "</span>"),
          grepl("^Failed", distribution_status) ~ paste0("<span class='status-badge fail'>", htmltools::htmlEscape(distribution_status), "</span>"),
          TRUE ~ paste0("<span class='status-badge warning'>", htmltools::htmlEscape(distribution_status), "</span>")
        )
      )
    DT::datatable(display, escape = FALSE, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  parsing_log <- reactive({
    req(rv$active_files_info)
    rv$parse_log %||% "No parsing information available."
  })
  output$parsingInfo <- renderPrint({ cat(paste(parsing_log(), collapse = "\n")) })

  sampleIdCol_active <- reactive({
    df <- tryCatch(filtered_data(), error = function(e) NULL)
    if (!is.null(df) && "sample_key" %in% names(df)) "sample_key" else "sample_name"
  })

  observe({
    df <- tryCatch(filtered_data(), error = function(e) NULL)
    req(df)
    numeric_metrics <- names(df %>% dplyr::select(where(is.numeric)))
    if (!length(numeric_metrics)) return()
    main_numeric_metrics <- numeric_metrics[!grepl("standard_error$|_se$|_standard_error_particles_ml$", numeric_metrics, ignore.case = TRUE)]
    if (!length(main_numeric_metrics)) main_numeric_metrics <- numeric_metrics
    metric_choices <- stats::setNames(main_numeric_metrics, format_col_name(main_numeric_metrics))
    preferred_m1 <- intersect(c("mean_average", "concentration_average_particles_ml"), main_numeric_metrics)
    sel_m1 <- if (length(preferred_m1)) preferred_m1[1] else main_numeric_metrics[1]
    sel_m2_pref <- intersect("mode_average", main_numeric_metrics)[1]
    sel_m2 <- if (!is.null(sel_m2_pref) && !is.na(sel_m2_pref) && sel_m2_pref != sel_m1) {
      sel_m2_pref
    } else if (length(main_numeric_metrics) > 1) {
      setdiff(main_numeric_metrics, sel_m1)[1] %||% main_numeric_metrics[1]
    } else {
      main_numeric_metrics[1]
    }
    updateSelectInput(session, "metricSelect1", choices = metric_choices, selected = preserve_metric(isolate(input$metricSelect1), main_numeric_metrics, sel_m1))
    updateSelectInput(session, "metricSelect2", choices = metric_choices, selected = preserve_metric(isolate(input$metricSelect2), main_numeric_metrics, sel_m2))
    updateSelectInput(session, "statsMetric", choices = metric_choices, selected = preserve_metric(isolate(input$statsMetric), main_numeric_metrics, sel_m1))
  }) |> bindEvent(filtered_data())

  output$downloadGroupTemplate <- downloadHandler(
    filename = function() "group_assignment_template.csv",
    content = function(file) {
      df <- filtered_data()
      req(df)
      validate(need(all(c("source_file", "sample_key", "sample_name", "well_id") %in% names(df)), "Sample identity columns are missing."))
      template_df <- df %>%
        dplyr::select(source_file, sample_key, sample_name, well_id) %>%
        dplyr::distinct() %>%
        dplyr::mutate(Group = "")
      if (!is.null(rv$group_map)) {
        template_df <- template_df %>% dplyr::select(-Group) %>%
          dplyr::left_join(rv$group_map, by = c("source_file", "sample_key")) %>%
          dplyr::mutate(Group = dplyr::coalesce(group, "")) %>% dplyr::select(-group)
      }
      readr::write_csv(template_df, file)
      showNotification("Group template generated.", type = "message")
    }
  )

  group_assignments_from_file <- reactive({
    req(input$groupFile)
    tryCatch({
      group_df <- readr::read_csv(input$groupFile$datapath, col_types = cols(.default = "c"), na = "", show_col_types = FALSE)
      required <- c("source_file", "sample_key", "Group")
      if (!all(required %in% names(group_df))) {
        stop("Group file must contain: ", paste(required, collapse = ", "))
      }

      group_df <- group_df %>%
        dplyr::transmute(
          source_file = trimws(source_file),
          sample_key = trimws(sample_key),
          group = trimws(Group)
        ) %>%
        dplyr::filter(!is.na(source_file), !is.na(sample_key), !is.na(group),
                      nzchar(source_file), nzchar(sample_key), nzchar(group))

      
      duplicate_keys <- group_df %>% dplyr::count(source_file, sample_key) %>% dplyr::filter(n > 1)
      if (nrow(duplicate_keys)) stop("Duplicate sample keys were found in the group file. Each sample may appear only once.")
      group_df
    }, error = function(e) {
      showNotification(paste("Error reading group file:", conditionMessage(e)), type = "error")
      NULL
    })
  })

  observeEvent(input$groupFile, {
    imported <- group_assignments_from_file()
    if (!is.null(imported)) rv$group_map <- imported
  })

  observeEvent(input$editGroups, {
    df <- tryCatch(filtered_data(), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) {
      showNotification("Load samples before assigning groups.", type = "warning")
      return()
    }
    samples <- df %>% dplyr::distinct(source_file, sample_key, sample_name, well_id)
    if (!is.null(rv$group_map)) samples <- dplyr::left_join(samples, rv$group_map, by = c("source_file", "sample_key"))
    if (!"group" %in% names(samples)) samples$group <- ""
    samples$group[is.na(samples$group)] <- ""
    rv$group_editor_samples <- samples
    rows <- lapply(seq_len(nrow(samples)), function(i) tags$tr(
      tags$td(samples$sample_name[i], style = "overflow-wrap:anywhere;"),
      tags$td(samples$well_id[i]),
      tags$td(samples$source_file[i], style = "overflow-wrap:anywhere;"),
      tags$td(group_name_input(paste0("sampleGroup_", i), paste("Group for", samples$sample_name[i], "well", samples$well_id[i]), samples$group[i]))
    ))
    showModal(modalDialog(
      title = "Assign statistical groups", size = "l", easyClose = FALSE,
      tags$datalist(id = "groupNameSuggestions", lapply(unique(samples$group[nzchar(samples$group)]), function(x) tags$option(value = x))),
      tags$p("Enter the same group name for samples that belong together. Blank fields leave samples unassigned. Only samples in the current filters are edited."),
      div(style = "max-height:60vh;overflow-y:auto;",
        tags$table(class = "table table-striped", style = "table-layout:fixed;width:100%;",
          tags$thead(tags$tr(tags$th("Sample"), tags$th("Well"), tags$th("Source file"), tags$th("Group"))),
          tags$tbody(rows))),
      footer = tagList(modalButton("Cancel"), actionButton("saveGroups", "Apply groups", class = "btn-primary"))
    ))
  })

  observeEvent(input$saveGroups, {
    samples <- rv$group_editor_samples
    req(samples)
    edited <- samples %>% dplyr::select(source_file, sample_key)
    edited$group <- vapply(seq_len(nrow(samples)), function(i) trimws(input[[paste0("sampleGroup_", i)]] %||% samples$group[i]), character(1))
    previous <- rv$group_map
    if (!is.null(previous)) previous <- dplyr::anti_join(previous, edited, by = c("source_file", "sample_key"))
    rv$group_map <- dplyr::bind_rows(previous, dplyr::filter(edited, nzchar(group)))
    rv$group_editor_samples <- NULL
    removeModal()
  })

  observe({
    df <- tryCatch(filtered_data(), error = function(e) NULL)
    group_map <- rv$group_map
    if (is.null(df) || is.null(group_map) || !nrow(group_map)) {
      rv$grouped_data <- NULL
      return()
    }

    join_cols <- c("source_file", "sample_key")
    df_grouped <- tryCatch(
      df %>%
        dplyr::inner_join(group_map, by = join_cols) %>%
        dplyr::mutate(group = droplevels(factor(group))),
      error = function(e) {
        showNotification(paste("Error joining groups:", conditionMessage(e)), type = "error")
        NULL
      }
    )

    if (is.null(df_grouped)) {
      rv$grouped_data <- NULL
      return()
    }
    if (nrow(df_grouped) > nrow(df)) {
      showNotification("Group assignment duplicated observations; grouping was rejected.", type = "error")
      rv$grouped_data <- NULL
      return()
    }
    if (!nrow(df_grouped)) {
      showNotification("No samples matched the group assignments.", type = "error")
      rv$grouped_data <- NULL
      return()
    }
    if (nrow(df_grouped) < nrow(df)) {
      showNotification("Some filtered samples were not assigned a group and were excluded from statistics.", type = "warning")
    }
    rv$grouped_data <- df_grouped
    showNotification("Groups applied.", type = "message")
  }) |> bindEvent(filtered_data(), rv$group_map)

  data_for_display <- reactive({
    df <- filtered_data()
    if (!is.null(rv$group_map)) df <- dplyr::left_join(df, rv$group_map, by = c("source_file", "sample_key"))
    req(df)
    priority <- intersect(c("source_file", "sample_name", "well_id", "group", "sample_key"), names(df))
    df %>% dplyr::relocate(dplyr::all_of(priority))
  })

  observe({
    df <- data_for_display()
    choices <- names(df)
    default_columns <- choices
    if (!isTRUE(input$showTechnicalColumns)) default_columns <- setdiff(default_columns, "sample_key")
    current <- isolate(input$visibleDataColumns)
    selected <- if (length(current)) intersect(current, choices) else default_columns
    if (!length(selected)) selected <- head(default_columns, min(8, length(default_columns)))
    updateSelectizeInput(session, "visibleDataColumns", choices = choices, selected = selected, server = TRUE)
  }) |> bindEvent(data_for_display(), input$showTechnicalColumns)

  output$parsedTable <- DT::renderDataTable({
    df_display <- data_for_display()
    req(df_display)
    selected_columns <- input$visibleDataColumns
    if (length(selected_columns)) {
      selected_columns <- intersect(selected_columns, names(df_display))
      if (length(selected_columns)) df_display <- df_display %>% dplyr::select(dplyr::all_of(selected_columns))
    }

    numeric_cols <- names(df_display)[vapply(df_display, is.numeric, logical(1))]
    concentration_cols <- numeric_cols[grepl("concentration|particles_ml", numeric_cols, ignore.case = TRUE)]
    standard_numeric <- setdiff(numeric_cols, concentration_cols)

    table_obj <- DT::datatable(
      df_display,
      extensions = c("Buttons", "FixedColumns"),
      options = list(
        scrollX = TRUE, pageLength = 12, dom = "Bfrtip",
        buttons = c("copy"), fixedColumns = list(leftColumns = min(3, ncol(df_display))),
        columnDefs = list(list(className = "dt-nowrap", targets = "_all"))
      ),
      rownames = FALSE, filter = "top"
    )
    if (length(concentration_cols)) table_obj <- DT::formatSignif(table_obj, columns = concentration_cols, digits = 3)
    if (length(standard_numeric)) table_obj <- DT::formatRound(table_obj, columns = standard_numeric, digits = 3)
    table_obj
  })
  
  qc_results <- reactive({
    req(input$qcFile$datapath, input$qcTargetSize, input$qcTolerance)
    parsed <- parse_nta_summary(input$qcFile$datapath)
    if (is.null(parsed$data) || !nrow(parsed$data)) {
      return(list(status = "ERROR", message = paste("QC parsing failed:", parsed$main_error %||% "unknown error"), source_file = input$qcFile$name))
    }

    qc_df <- parsed$data
    required_metrics <- c("mean_average", "mode_average")
    missing_metrics <- setdiff(required_metrics, names(qc_df))
    if (length(missing_metrics)) {
      return(list(status = "ERROR", message = paste("QC file is missing:", paste(missing_metrics, collapse = ", ")), source_file = input$qcFile$name))
    }

    validate(need(is.finite(input$qcTargetSize) && input$qcTargetSize > 0, "Enter a positive QC target."),
             need(is.finite(input$qcTolerance) && input$qcTolerance >= 0 && input$qcTolerance <= 100, "Enter a QC tolerance from 0 to 100%."))
    target <- input$qcTargetSize
    tolerance_percent <- input$qcTolerance
    tolerance_fraction <- tolerance_percent / 100
    lower <- target * (1 - tolerance_fraction)
    upper <- target * (1 + tolerance_fraction)

    details <- qc_df %>%
      dplyr::transmute(
        sample_name = as.character(sample_name),
        well_id = as.character(well_id),
        mean_value = mean_average,
        mode_value = mode_average,
        mean_se = if ("mean_standard_error" %in% names(qc_df)) mean_standard_error else NA_real_,
        mode_se = if ("mode_standard_error" %in% names(qc_df)) mode_standard_error else NA_real_,
        concentration_value = if ("concentration_average_particles_ml" %in% names(qc_df)) concentration_average_particles_ml else NA_real_,
        concentration_se = if ("concentration_standard_error_particles_ml" %in% names(qc_df)) concentration_standard_error_particles_ml else NA_real_,
        mean_check = !is.na(mean_value) & mean_value >= lower & mean_value <= upper,
        mode_check = !is.na(mode_value) & mode_value >= lower & mode_value <= upper,
        row_status = dplyr::case_when(
          is.na(mean_value) | is.na(mode_value) ~ "ERROR",
          mean_check & mode_check ~ "PASS",
          TRUE ~ "FAIL"
        )
      )

    overall_status <- if (any(details$row_status == "ERROR")) "ERROR" else if (all(details$row_status == "PASS")) "PASS" else "FAIL"
    first_row <- details[1, , drop = FALSE]

    list(
      status = overall_status,
      message = if (overall_status == "ERROR") "One or more QC wells contain missing Mean or Mode values." else NULL,
      details = details,
      n_qc = nrow(details),
      mean_value = round(first_row$mean_value, 2), mean_check = first_row$mean_check,
      mode_value = round(first_row$mode_value, 2), mode_check = first_row$mode_check,
      target = target, lower_bound = lower, upper_bound = upper,
      qc_sample_name = first_row$sample_name, qc_well_id = first_row$well_id,
      tolerance = tolerance_percent, mean_se = first_row$mean_se, mode_se = first_row$mode_se,
      time_captured = parsed$metadata$time_captured %||% NA_character_,
      source_file = basename(input$qcFile$name)
    )
  })

  output$qcStatusDisplay <- renderUI({
    if (is.null(input$qcFile$datapath)) {
      return(div(
        class = "empty-state",
        icon("check-circle"),
        h4("No QC summary uploaded"),
        p("Upload an optional QC summary in the Quality control section to evaluate mean and mode size against the target range.")
      ))
    }

    results <- tryCatch(
      qc_results(),
      error = function(e) list(status = "ERROR", message = conditionMessage(e), source_file = basename(input$qcFile$name))
    )
    req(results)
    if (identical(results$status, "ERROR") && is.null(results$details)) {
      return(div(class = "warning-banner", icon("exclamation-triangle"), tags$strong("QC error: "), results$message))
    }

    details <- results$details
    mean_pass <- sum(details$mean_check %in% TRUE, na.rm = TRUE)
    mode_pass <- sum(details$mode_check %in% TRUE, na.rm = TRUE)
    overall_class <- switch(results$status, PASS = "status-pass", FAIL = "status-fail", ERROR = "status-warning", "status-neutral")
    overall_symbol <- switch(results$status, PASS = "✓", FAIL = "✕", ERROR = "!", "—")

    tagList(
      div(
        class = "metric-grid",
        metric_card("Overall QC", paste(overall_symbol, results$status), paste(results$n_qc, "QC wells assessed"), overall_class),
        metric_card("Mean size", paste0(mean_pass, "/", results$n_qc), "Wells within target range", if (mean_pass == results$n_qc) "status-pass" else "status-fail"),
        metric_card("Mode size", paste0(mode_pass, "/", results$n_qc), "Wells within target range", if (mode_pass == results$n_qc) "status-pass" else "status-fail"),
        metric_card("Accepted range", paste0(round(results$lower_bound, 1), "–", round(results$upper_bound, 1)), paste0("Target ", results$target, " nm ±", results$tolerance, "%"))
      ),
      div(
        class = "info-banner",
        icon("file-lines"),
        tags$strong(results$source_file),
        paste0(" — captured: ", results$time_captured %||% "not reported")
      ),
      if (!is.null(results$message)) div(class = "warning-banner", icon("exclamation-triangle"), results$message)
    )
  })

  output$qcDetailsTable <- DT::renderDataTable({
    req(input$qcFile$datapath)
    results <- qc_results()
    req(results$details)
    display <- results$details %>%
      dplyr::transmute(
        Sample = sample_name, Well = well_id,
        `Mean (nm)` = round(mean_value, 2),
        `Mean check` = dplyr::case_when(is.na(mean_value) ~ "ERROR", mean_check ~ "PASS", TRUE ~ "FAIL"),
        `Mode (nm)` = round(mode_value, 2),
        `Mode check` = dplyr::case_when(is.na(mode_value) ~ "ERROR", mode_check ~ "PASS", TRUE ~ "FAIL"),
        Status = row_status
      )
    DT::datatable(display, options = list(dom = "t", scrollX = TRUE), rownames = FALSE) %>%
      DT::formatStyle(
        c("Mean check", "Mode check", "Status"),
        color = DT::styleEqual(c("PASS", "FAIL", "ERROR"), c("#2e7d32", "#b3261e", "#a86400")),
        fontWeight = "bold"
      )
  })

  qc_current_signature <- reactive({
    req(input$qcFile$datapath)
    list(file = input$qcFile$datapath, target = input$qcTargetSize,
         tolerance = input$qcTolerance, directory = input$qcLogDirPath)
  })
  output$qcSaveStatus <- renderUI({
    if (is.null(input$qcFile)) return(NULL)
    if (identical(rv$qc_saved_signature, qc_current_signature())) {
      div(class = "small-muted", paste("Current assessment saved to", rv$qc_saved_path))
    } else div(class = "small-muted", "Current assessment is not saved. Review the criteria, then click Save QC result.")
  })
  observeEvent(input$saveQcResult, {
    req(input$qcFile$datapath)
    if (identical(rv$qc_saved_signature, qc_current_signature())) {
      showNotification("This assessment is already saved.", type = "message")
      return()
    }
    parsed <- parse_nta_summary(input$qcFile$datapath)
    if (is.null(parsed$data) || !nrow(parsed$data)) {
      showNotification("Failed to parse the uploaded QC file for logging.", type = "warning")
      return()
    }

    results <- tryCatch(qc_results(), error = function(e) NULL)
    if (is.null(results) || is.null(results$details)) {
      showNotification("QC results were unavailable for logging.", type = "warning")
      return()
    }

    user_dir_path <- trimws(input$qcLogDirPath %||% ".")
    if (!nzchar(user_dir_path) || identical(user_dir_path, ".")) {
      log_dir <- getwd()
    } else {
      log_dir <- normalizePath(user_dir_path, mustWork = FALSE)
      if (!dir.exists(log_dir)) {
        showNotification("QC Log Directory is invalid.", type = "error")
        return()
      }
    }

    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    server_user <- unname(Sys.info()["user"])
    if (is.na(server_user) || !nzchar(server_user)) server_user <- "Unknown"
    qc_filename <- basename(input$qcFile$name %||% "N/A")

    log_entry <- results$details %>%
      dplyr::transmute(
        Timestamp = timestamp,
        User = server_user,
        QcFile = qc_filename,
        QcSample = sample_name,
        QcWell = well_id,
        MeanAverage = mean_value,
        ModeAverage = mode_value,
        MeanSE = mean_se,
        ModeSE = mode_se,
        ConcentrationAverage = concentration_value,
        ConcentrationSE = concentration_se,
        TimeCaptured = results$time_captured,
        Status = row_status,
        TargetSize = results$target,
        TolerancePercent = results$tolerance,
        AppVersion = APP_VERSION
      )

    tryCatch({
      log_file_path <- file.path(log_dir, QC_LOG_FILE_BASENAME)
      if (file.exists(log_file_path) && file.size(log_file_path) > 0) {
        existing_header <- names(readr::read_csv(log_file_path, n_max = 0, show_col_types = FALSE))
        if (!identical(existing_header, names(log_entry))) {
          log_file_path <- file.path(log_dir, "nta_qc_log_v2.csv")
          if (file.exists(log_file_path) && file.size(log_file_path) > 0 &&
              !identical(names(readr::read_csv(log_file_path, n_max = 0, show_col_types = FALSE)), names(log_entry))) {
            stop("QC log format differs. Choose a new log directory to save this assessment.")
          }
        }
      }
      write_headers <- !file.exists(log_file_path) || file.size(log_file_path) == 0
      readr::write_csv(log_entry, log_file_path, append = !write_headers, col_names = write_headers, na = "")
      rv$qc_saved_signature <- qc_current_signature()
      rv$qc_saved_path <- log_file_path
      showNotification(paste(nrow(log_entry), "QC row(s) saved to", basename(log_file_path)), type = "message", duration = 7)
    }, error = function(e) {
      showNotification(paste("Error writing QC log:", conditionMessage(e)), type = "error")
    })
  }, ignoreNULL = TRUE, ignoreInit = TRUE)

  comparison_data <- reactive({
    df <- filtered_data()
    req(df)
    if (!is.null(rv$grouped_data) && nrow(rv$grouped_data)) {
      group_map <- rv$grouped_data %>% dplyr::distinct(source_file, sample_key, group)
      if (!"group" %in% names(df)) df <- df %>% dplyr::left_join(group_map, by = c("source_file", "sample_key"))
    }
    df
  })

  observe({
    df <- tryCatch(comparison_data(), error = function(e) NULL)
    has_groups <- !is.null(df) && "group" %in% names(df) && any(!is.na(df$group))
    x_choices <- c("Sample" = "sample_key", "Source file" = "source_file")
    color_choices <- c("Sample" = "sample_key", "Source file" = "source_file", "None" = "none")
    if (has_groups) {
      x_choices <- c(x_choices, "Statistical group" = "group")
      color_choices <- c(color_choices[1:2], "Statistical group" = "group", "None" = "none")
    }
    current_x <- isolate(input$comparisonXAxis) %||% "sample_key"
    current_color <- isolate(input$comparisonColorBy) %||% "sample_key"
    if (!current_x %in% unname(x_choices)) current_x <- "sample_key"
    if (!current_color %in% unname(color_choices)) current_color <- "sample_key"
    updateSelectInput(session, "comparisonXAxis", choices = x_choices, selected = current_x)
    updateSelectInput(session, "comparisonColorBy", choices = color_choices, selected = current_color)
  }) |> bindEvent(comparison_data(), rv$grouped_data)

  comparison_sample_labels <- reactive({
    # Derive IDs before filtering/sorting so they remain stable across both plots.
    df <- parsed_data_raw() %>% dplyr::distinct(sample_key, sample_name, well_id, source_file)
    df$plot_label <- comparison_display_names(df$sample_name)
    duplicate_names <- duplicated(df$plot_label) | duplicated(df$plot_label, fromLast = TRUE)
    df$plot_label[duplicate_names] <- paste0(df$plot_label[duplicate_names], " [", df$well_id[duplicate_names], "]")
    df$plot_label <- make.unique(df$plot_label, sep = " #")
    df
  })

  observeEvent(input$showSampleNames, {
    req(comparison_sample_labels())
    showModal(modalDialog(title = "Comparison sample names", size = "l",
      tags$p("Samples in the current filters. Samples without values for a selected plot metric do not appear in that plot."),
      DTOutput("comparisonSampleNames"), easyClose = TRUE, footer = modalButton("Close")))
  })
  output$comparisonSampleNames <- renderDT({
    df <- comparison_sample_labels() %>% dplyr::filter(sample_key %in% filtered_data()$sample_key)
    names(df) <- c("Sample key", "Full sample name", "Well", "Source file", "Plot label")
    DT::datatable(df[c("Plot label", "Full sample name", "Well", "Source file", "Sample key")],
      rownames = FALSE, options = list(pageLength = 10, scrollX = TRUE))
  })

  comparison_plot_height <- reactive({
    switch(input$plotSize %||% "standard", compact = 380, large = 760, 560)
  })

  generate_comparison_plot <- reactive({
    function(df, metric_col, x_col, color_col, plot_type, use_log, sort_mode, show_labels) {
      req(df, metric_col, x_col, plot_type)
      if (!metric_col %in% names(df) || !x_col %in% names(df)) return(NULL)
      if (!is.numeric(df[[metric_col]])) return(NULL)
      if (!identical(color_col, "none") && !color_col %in% names(df)) color_col <- "none"

      plot_df <- df %>% dplyr::filter(!is.na(.data[[metric_col]]), !is.na(.data[[x_col]]))
      if (!nrow(plot_df)) return(NULL)
      plot_df[[x_col]] <- as.character(plot_df[[x_col]])

      x_levels <- unique(plot_df[[x_col]])
      if (identical(sort_mode, "alphabetical")) {
        x_levels <- sort(x_levels)
      } else if (identical(sort_mode, "value")) {
        order_df <- plot_df %>%
          dplyr::group_by(.data[[x_col]]) %>%
          dplyr::summarise(order_value = mean(.data[[metric_col]], na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(order_value)
        x_levels <- as.character(order_df[[x_col]])
      }
      plot_df[[x_col]] <- factor(plot_df[[x_col]], levels = x_levels)
      if (!identical(color_col, "none")) plot_df[[color_col]] <- factor(plot_df[[color_col]])

      derived_error_col <- NULL
      if (identical(plot_type, "points")) {
        candidates <- if (metric_col == "concentration_average_particles_ml") {
          "concentration_standard_error_particles_ml"
        } else if (stringr::str_detect(metric_col, "_average$")) {
          stringr::str_replace(metric_col, "_average$", "_standard_error")
        } else if (stringr::str_detect(metric_col, "_(mean|mode)$")) {
          stringr::str_replace(metric_col, "_(mean|mode)$", "_standard_error")
        } else character()
        if (length(candidates) && candidates[1] %in% names(plot_df) && is.numeric(plot_df[[candidates[1]]])) derived_error_col <- candidates[1]
      }

      base_mapping <- if (identical(color_col, "none")) {
        ggplot2::aes(x = .data[[x_col]], y = .data[[metric_col]])
      } else {
        ggplot2::aes(x = .data[[x_col]], y = .data[[metric_col]], color = .data[[color_col]])
      }
      p <- ggplot2::ggplot(plot_df, base_mapping)

      if (identical(plot_type, "points")) {
        p <- p + ggplot2::geom_point(size = 3, alpha = 0.88, na.rm = TRUE)
        if (!is.null(derived_error_col)) {
          error_mapping <- if (identical(color_col, "none")) {
            ggplot2::aes(ymin = .data[[metric_col]] - .data[[derived_error_col]], ymax = .data[[metric_col]] + .data[[derived_error_col]])
          } else {
            ggplot2::aes(ymin = .data[[metric_col]] - .data[[derived_error_col]], ymax = .data[[metric_col]] + .data[[derived_error_col]], color = .data[[color_col]])
          }
          p <- p + ggplot2::geom_errorbar(mapping = error_mapping, width = 0.22, linewidth = 0.7, na.rm = TRUE, show.legend = FALSE)
        }
      } else if (identical(plot_type, "box")) {
        fill_mapping <- if (identical(color_col, "none")) NULL else ggplot2::aes(fill = .data[[color_col]])
        p <- p +
          ggplot2::geom_boxplot(mapping = fill_mapping, alpha = 0.55, outlier.shape = NA, na.rm = TRUE, show.legend = FALSE) +
          ggplot2::geom_jitter(width = 0.14, height = 0, alpha = 0.72, na.rm = TRUE)
      } else if (identical(plot_type, "violin")) {
        fill_mapping <- if (identical(color_col, "none")) NULL else ggplot2::aes(fill = .data[[color_col]])
        p <- p +
          ggplot2::geom_violin(mapping = fill_mapping, alpha = 0.48, scale = "width", trim = FALSE, na.rm = TRUE, show.legend = FALSE) +
          ggplot2::geom_boxplot(width = 0.11, fill = "white", alpha = 0.65, outlier.shape = NA, na.rm = TRUE) +
          ggplot2::geom_jitter(width = 0.08, height = 0, alpha = 0.65, na.rm = TRUE)
      }

      if (!identical(color_col, "none")) {
        p <- p + viridis::scale_color_viridis(discrete = TRUE, option = "plasma", begin = 0.08, end = 0.92,
          labels = if (identical(color_col, "sample_key")) {
            lookup <- comparison_sample_labels()
            function(x) lookup$plot_label[match(x, lookup$sample_key)]
          } else ggplot2::waiver())
        if (plot_type %in% c("box", "violin")) p <- p + viridis::scale_fill_viridis(discrete = TRUE, option = "plasma", begin = 0.08, end = 0.92)
      }

      subtitle <- if (!is.null(derived_error_col)) paste("Error bars: ±", format_col_name(derived_error_col)) else NULL
      p <- p +
        ggplot2::labs(
          subtitle = subtitle,
          x = if (identical(x_col, "sample_key")) "Sample" else format_col_name(x_col), y = format_col_name(metric_col),
          color = if (identical(color_col, "none")) NULL else format_col_name(color_col)
        ) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          axis.text.x = if (isTRUE(show_labels)) ggplot2::element_text(size = 11, colour = "grey15",
            angle = 45, vjust = 1, hjust = 1) else ggplot2::element_blank(),
          legend.text = ggplot2::element_text(size = 13, colour = "grey15"),
          legend.title = ggplot2::element_text(size = 13),
          plot.margin = ggplot2::margin(10, 15, 10, 40),
          axis.ticks.x = if (isTRUE(show_labels)) ggplot2::element_line() else ggplot2::element_blank(),
          plot.subtitle = ggplot2::element_text(size = ggplot2::rel(0.9), colour = "grey40"),
          legend.position = if (!identical(color_col, "none") && !identical(color_col, x_col)) "bottom" else "none"
        ) +
        ggplot2::scale_x_discrete(labels = function(x) {
          if (identical(x_col, "sample_key")) {
            lookup <- comparison_sample_labels()
            lookup$plot_label[match(x, lookup$sample_key)]
          } else stringr::str_wrap(x, width = 24)
        })


      if (isTRUE(use_log)) {
        low_check <- plot_df[[metric_col]]
        if (!is.null(derived_error_col)) low_check <- low_check - plot_df[[derived_error_col]]
        if (any(low_check <= 0, na.rm = TRUE)) {
          showNotification(paste(format_col_name(metric_col), "contains zero or negative values; Log10 was not applied."), type = "warning")
        } else {
          p <- p + ggplot2::scale_y_log10() + ggplot2::labs(y = paste(format_col_name(metric_col), "(Log10)"))
        }
      }
      p
    }
  })

  plot_object1 <- reactive({
    req(comparison_data(), input$metricSelect1, input$comparisonXAxis)
    generate_comparison_plot()(
      comparison_data(), input$metricSelect1, input$comparisonXAxis,
      input$comparisonColorBy %||% "none", input$comparisonPlotType %||% "points",
      isTRUE(input$comparisonLogY), input$comparisonSort %||% "original",
      isTRUE(input$comparisonShowLabels)
    )
  })

  plot_object2 <- reactive({
    req(comparison_data(), input$metricSelect2, input$comparisonXAxis)
    generate_comparison_plot()(
      comparison_data(), input$metricSelect2, input$comparisonXAxis,
      input$comparisonColorBy %||% "none", input$comparisonPlotType %||% "points",
      isTRUE(input$comparisonLogY), input$comparisonSort %||% "original",
      isTRUE(input$comparisonShowLabels)
    )
  })

  output$plotTitle1 <- renderText({ req(input$metricSelect1); format_col_name(input$metricSelect1) })
  output$plotTitle2 <- renderText({ req(input$metricSelect2); format_col_name(input$metricSelect2) })
  output$comparisonPlot1 <- renderPlot({ plot_object1() }, height = function() comparison_plot_height())
  output$comparisonPlot2 <- renderPlot({ plot_object2() }, height = function() comparison_plot_height())

  output$comparisonPlotsUI <- renderUI({
    if (is.null(rv$parsed_result) || is.null(rv$parsed_result$data)) {
      return(div(class = "nta-card empty-state", icon("bar-chart"), h4("Upload NTA summary files to create comparison plots.")))
    }
    warning_ui <- NULL
    if ((input$comparisonPlotType %||% "points") %in% c("box", "violin") && identical(input$comparisonXAxis, "sample_key")) {
      warning_ui <- div(class = "warning-banner", icon("exclamation-triangle"), "Box and violin plots are most informative when the X-axis is Source file or Statistical group, because individual samples often contain only one value.")
    }
    card1 <- div(
      class = "nta-card plot-card",
      div(class = "card-header-row", h4(textOutput("plotTitle1", inline = TRUE)), div(class = "card-actions", downloadButton("downloadPlot1", "Download PNG", class = "btn-default btn-sm"))),
      plotOutput("comparisonPlot1")
    )
    card2 <- div(
      class = "nta-card plot-card",
      div(class = "card-header-row", h4(textOutput("plotTitle2", inline = TRUE)), div(class = "card-actions", downloadButton("downloadPlot2", "Download PNG", class = "btn-default btn-sm"))),
      plotOutput("comparisonPlot2")
    )
    if (identical(input$comparisonLayout, "stack")) {
      tagList(warning_ui, card1, card2)
    } else {
      tagList(warning_ui, fluidRow(column(6, card1), column(6, card2)))
    }
  })


  # --- Statistics ---
  stats_available_counts <- reactive({
    grouped <- rv$grouped_data
    metric <- input$statsMetric
    if (is.null(grouped) || !nrow(grouped) || is.null(metric) ||
        !metric %in% names(grouped) || !is.numeric(grouped[[metric]])) return(NULL)
    grouped %>% dplyr::filter(!is.na(group), nzchar(trimws(as.character(group)))) %>% dplyr::group_by(group) %>%
      dplyr::summarise(n = sum(!is.na(.data[[metric]])), missing = sum(is.na(.data[[metric]])), .groups = "drop")
  })

  output$statsSelectionHelp <- renderUI({
    choice <- input$statsTestType %||% "Auto"
    message <- switch(choice,
      Auto = "Auto uses within-group normality and variance checks. The executed test appears below; the selection reason is in full statistical details. Screening is guidance, not a guarantee.",
      aov = "Manual mean comparison: assumes equal variances. Your chosen test will run without automatic switching.",
      welch_aov = "Manual mean comparison: allows unequal variances and group sizes. Your chosen test will run without automatic switching.",
      t.test = "Manual comparison of two group means: allows unequal variances and group sizes.",
      "Manual rank comparison: does not test equality of arithmetic means. A median/location interpretation requires comparable distribution shapes.")
    div(class = "small-muted", style = "margin-bottom:10px;", role = "status", message)
  })

  output$statsPrerequisiteMessage <- renderUI({
    if (is.null(rv$parsed_result) || is.null(rv$parsed_result$data)) {
      return(div(class = "info-banner", icon("info-circle"), "Upload NTA summary files before configuring statistics."))
    }
    if (is.null(rv$grouped_data) || !nrow(rv$grouped_data)) {
      return(div(
        class = "warning-banner",
        icon("users"), tags$strong("Groups are required. "),
        "Download the group template from the sidebar, complete the Group column and upload it to enable statistical comparisons."
      ))
    }
    group_counts <- stats_available_counts()
    if (is.null(group_counts)) return(div(class = "info-banner", "Select a particle parameter to see the usable group sizes."))
    usable <- group_counts$n[group_counts$n > 0]
    tagList(
      div(class = "info-banner", role = "status",
          tags$strong("Usable samples for this parameter: "),
          paste(paste0(group_counts$group, ": n = ", group_counts$n), collapse = "; "),
          if (length(unique(usable)) > 1L) div("Unequal group sizes are supported; all usable samples are retained."),
          if (sum(group_counts$missing) > 0L) div(paste(sum(group_counts$missing), "missing values excluded."))),
      if (any(group_counts$n == 0L)) div(class = "warning-banner", "Groups with n = 0 have no usable values and are excluded from this parameter's comparison."),
      if (length(usable) < 2L) div(class = "warning-banner", "At least two groups with usable values are required."),
      if (length(usable) >= 2L && any(usable < 2L)) div(class = "warning-banner", "Each included group needs at least two usable observations. No test can run yet."),
      if (any(usable >= 2L & usable < 5L)) div(class = "warning-banner", icon("exclamation-triangle"), "Some groups have fewer than five usable samples. Normality checks and estimates are uncertain; a small n alone does not determine the test.")
    )
  })

  observe({
    counts <- stats_available_counts()
    choices <- statistical_test_choices(if (is.null(counts)) NULL else sum(counts$n > 0L))
    selected <- isolate(input$statsTestType) %||% "Auto"
    if (!selected %in% unname(unlist(choices))) {
      selected <- "Auto"
      showNotification("The number of usable groups changed. Test selection has been reset to Auto.", type = "message")
    }
    updateSelectInput(session, "statsTestType", choices = choices, selected = selected)
  }) |> bindEvent(stats_available_counts())

  stats_data_prep <- reactive({
    data_grouped <- rv$grouped_data
    metric <- input$statsMetric
    req(data_grouped, metric, nrow(data_grouped) > 0)
    validate(
      need(identical(input$statsDesign, "independent"),
           if (identical(input$statsDesign, "paired"))
             "Paired or repeated measurements need a paired/repeated-measures analysis, which this version does not provide."
           else "Select the sample relationship before running statistics."),
      need(metric %in% names(data_grouped), "Selected metric is missing."),
      need(is.numeric(data_grouped[[metric]]), "Selected metric is not numeric.")
    )

    validate(need(!any(is.infinite(data_grouped[[metric]])),
                  "The selected parameter contains infinite values. Correct them before analysis."))
    data_analysis <- data_grouped %>% dplyr::filter(!is.na(group), nzchar(trimws(as.character(group))), !is.na(.data[[metric]]))
    validate(need(nrow(data_analysis) > 0, "No non-missing values are available for statistics."))
    data_analysis$group <- droplevels(factor(data_analysis$group))
    group_counts <- data_analysis %>% dplyr::count(group, name = "n")
    num_groups <- nrow(group_counts)
    min_n <- min(group_counts$n)
    validate(
      need(num_groups >= 2, "At least two groups are required."),
      need(min_n >= 2, "Each group requires at least two non-missing observations.")
    )
    list(data = data_analysis, metric = metric, num_groups = num_groups, min_n_per_group = min_n)
  })

  run_statistical_test <- reactive({
    prep <- stats_data_prep()
    req(prep)
    data_analysis <- prep$data
    metric <- prep$metric
    num_groups <- prep$num_groups
    min_n <- prep$min_n_per_group
    formula <- stats::as.formula(paste0("`", metric, "` ~ group"))

    normality <- assess_group_normality(data_analysis, metric)
    normality_p <- stats::setNames(normality$shapiro_p, normality$group)
    minimum_shapiro_p <- if (any(is.finite(normality_p))) min(normality_p[is.finite(normality_p)]) else NA_real_
    levene_p <- tryCatch(car::leveneTest(formula, data = data_analysis, center = stats::median)$`Pr(>F)`[1], error = function(e) NA_real_)

    test_choice <- input$statsTestType %||% "Auto"
    selection <- select_statistical_test(test_choice, num_groups, min_n, normality_p, levene_p)
    actual_test <- selection$test_name
    selection_reason <- selection$selection_reason

    results <- list(
      test_name = actual_test, test_label = selection$test_label,
      selection_heading = selection$selection_heading,
      selection_reason = selection_reason, error = selection$error,
      post_hoc_method = selection$post_hoc_method,
      study_design = "Independent samples (selected design)",
      inference_note = "Each parameter is tested separately; p-values are not adjusted across different parameters.",
      assumptions = list(shapiro_p = minimum_shapiro_p, levene_p = levene_p,
                         normality_by_group = normality, normality_status = selection$normality_status),
      main_result = NULL, main_result_print = NULL,
      post_hoc_result = NULL, effect_size = NULL, data = data_analysis
    )

    if (!is.null(results$error)) return(results)

    if (actual_test %in% c("t.test", "wilcox.test") && num_groups != 2) {
      results$error <- "The selected test requires exactly two groups."
      return(results)
    }
    if (actual_test %in% c("aov", "welch_aov", "kruskal.test") && num_groups <= 2) {
      results$error <- "The selected test requires more than two groups."
      return(results)
    }

    if (actual_test == "welch_aov") {
      variances <- tapply(data_analysis[[metric]], data_analysis$group, stats::var)
      if (any(!is.finite(variances) | variances <= 0)) {
        results$error <- "Welch ANOVA requires a positive, finite variance in every group. Check groups with identical values; no alternative test has been substituted."
        return(results)
      }
    }

    tryCatch({
      if (actual_test == "t.test") {
        test_obj <- stats::t.test(formula, data = data_analysis, var.equal = FALSE)
        results$main_result <- broom::tidy(test_obj)
        results$effect_size <- tryCatch(cohens_d_unpooled(data_analysis, metric), error = function(e) NULL)
        group_levels <- levels(data_analysis$group)
        results$post_hoc_result <- tibble::tibble(
          group1 = group_levels[1], group2 = group_levels[2],
          p.value = test_obj$p.value, p.signif = significance_code(test_obj$p.value),
          adjustment = "None (single comparison)"
        )
      } else if (actual_test == "wilcox.test") {
        group_sizes <- table(data_analysis$group)
        use_exact <- all(group_sizes < 50) && !anyDuplicated(data_analysis[[metric]])
        test_obj <- stats::wilcox.test(formula, data = data_analysis, exact = use_exact)
        results$inference_note <- paste(results$inference_note,
          if (use_exact) "Wilcoxon uses an exact two-sided p-value (small samples without ties)."
          else "Wilcoxon uses a normal approximation with continuity correction because of ties or sample size; the approximation may be unreliable for very small groups.")
        results$main_result <- broom::tidy(test_obj)
        results$effect_size <- tryCatch(rank_biserial_from_wilcox(test_obj, data_analysis), error = function(e) NULL)
        group_levels <- levels(data_analysis$group)
        results$post_hoc_result <- tibble::tibble(
          group1 = group_levels[1], group2 = group_levels[2],
          p.value = test_obj$p.value, p.signif = significance_code(test_obj$p.value),
          adjustment = "None (single comparison)"
        )
      } else if (actual_test == "aov") {
        aov_obj <- stats::aov(formula, data = data_analysis)
        results$main_result <- broom::tidy(aov_obj)
        group_row <- results$main_result %>% dplyr::filter(term == "group")
        if (nrow(group_row) == 1) {
          results$main_result_print <- sprintf(
            "ANOVA F(%d, %d) = %.3f, p = %.4g",
            group_row$df, stats::df.residual(aov_obj), group_row$statistic, group_row$p.value
          )
          if (!is.na(group_row$p.value) && group_row$p.value < 0.05) {
            tukey_df <- as.data.frame(stats::TukeyHSD(aov_obj)$group)
            tukey_df$contrast <- rownames(tukey_df)
            split_groups <- t(vapply(tukey_df$contrast, split_tukey_contrast, FUN.VALUE = c(group1 = "", group2 = ""), group_levels = levels(data_analysis$group)))
            results$post_hoc_result <- tibble::as_tibble(tukey_df) %>%
              dplyr::mutate(
                group1 = split_groups[, "group1"], group2 = split_groups[, "group2"],
                p.adj = `p adj`, p.adj.signif = significance_code(p.adj),
                estimate = diff, conf.low = lwr, conf.high = upr
              ) %>%
              dplyr::filter(!is.na(group1), !is.na(group2)) %>%
              dplyr::select(group1, group2, p.adj, p.adj.signif, estimate, conf.low, conf.high)
          }
        }
        results$effect_size <- tryCatch(anova_effect_sizes(aov_obj), error = function(e) NULL)
      } else if (actual_test == "welch_aov") {
        welch_obj <- stats::oneway.test(formula, data = data_analysis, var.equal = FALSE)
        results$main_result <- broom::tidy(welch_obj)
        results$effect_size <- welch_effect_sizes(welch_obj)
        results$main_result_print <- sprintf(
          "Welch ANOVA F(%.2f, %.2f) = %.3f, p = %.4g",
          unname(welch_obj$parameter[1]), unname(welch_obj$parameter[2]),
          unname(welch_obj$statistic), welch_obj$p.value
        )
        if (!is.na(welch_obj$p.value) && welch_obj$p.value < 0.05) {
          results$post_hoc_result <- games_howell_test(data_analysis, metric)
        }
      } else if (actual_test == "kruskal.test") {
        test_obj <- stats::kruskal.test(formula, data = data_analysis)
        results$inference_note <- paste(results$inference_note,
          "Kruskal-Wallis uses a chi-squared approximation; its p-value may be unreliable for very small groups.")
        results$main_result <- broom::tidy(test_obj)
        if (!is.na(test_obj$p.value) && test_obj$p.value < 0.05) {
          results$post_hoc_result <- tryCatch(
            dunn_holm_test(data_analysis[[metric]], data_analysis$group),
            error = function(e) paste("Dunn post-hoc test failed:", conditionMessage(e))
          )
        }
        results$effect_size <- tryCatch(kruskal_eta_squared_h(test_obj, nrow(data_analysis), num_groups), error = function(e) NULL)
      }
    }, error = function(e) {
      results$error <<- paste("Statistical test error:", conditionMessage(e))
    })
    main_p <- if (actual_test == "aov") results$main_result$p.value[results$main_result$term == "group"] else results$main_result$p.value
    if (is.null(results$error) && !is.null(results$main_result) &&
        (!length(main_p) || any(!is.finite(main_p)))) {
      results$error <- "The selected test could not produce a finite p-value. Check for constant or degenerate data."
    }
    if (is.null(results$effect_size)) {
      results$effect_size <- tibble::tibble(effect = "Not available", estimate = NA_real_,
        note = if (!is.null(results$error)) results$error else "An effect-size estimate could not be computed for these data.")
    }
    results
  })

  stats_results <- reactive({
    prep <- stats_data_prep()
    stat_results <- run_statistical_test()
    req(prep, stat_results)

    results_capture <- capture.output({
      cat("===========================================\nStatistical Comparison Results\n===========================================\n\n")
      cat("Metric:", format_col_name(prep$metric), "\n")
      cat("Groups:", paste(levels(prep$data$group), collapse = ", "), "\n")
      cat("Study design:", stat_results$study_design, "\n")
      cat("Requested test:", input$statsTestType, "\nExecuted test:", stat_results$test_label, "\n\n")
      cat(stat_results$selection_heading, ":\n", sep = "")
      cat(paste(strwrap(stat_results$selection_reason, width = 88), collapse = "\n"), "\n\n")
      cat(paste(strwrap(stat_results$inference_note, width = 88), collapse = "\n"), "\n\n")
      cat("--- Within-group Normality Diagnostics ---\n")
      print(as.data.frame(stat_results$assumptions$normality_by_group), row.names = FALSE)
      cat("\n")

      if (!is.null(stat_results$error)) {
        cat("ERROR:", stat_results$error, "\n")
      } else {
        cat("--- Assumption Checks ---\n")
        shapiro_p <- stat_results$assumptions$shapiro_p
        levene_p <- stat_results$assumptions$levene_p
        cat("Smallest available within-group Shapiro p:", if (is.na(shapiro_p)) "N/A" else sprintf("p = %.4g%s", shapiro_p, if (shapiro_p <= 0.05) " (warning)" else ""), "\n")
        cat("Median-centred Levene variance test:", if (is.na(levene_p)) "N/A" else sprintf("p = %.4g%s", levene_p, if (levene_p <= 0.05) " (unequal variances)" else ""), "\n\n")

        cat("--- Main Test Results ---\n")
        main_res <- stat_results$main_result
        if (!is.null(stat_results$main_result_print)) {
          cat(stat_results$main_result_print, "\n")
        } else if (!is.null(main_res) && stat_results$test_name == "t.test") {
          cat(sprintf("Welch t(%.2f) = %.3f, p = %.4g\n", main_res$parameter, main_res$statistic, main_res$p.value))
        } else if (!is.null(main_res) && stat_results$test_name == "wilcox.test") {
          cat(sprintf("Wilcoxon W = %.1f, p = %.4g\n", main_res$statistic, main_res$p.value))
        } else if (!is.null(main_res) && stat_results$test_name == "kruskal.test") {
          cat(sprintf("Kruskal-Wallis chi-squared = %.3f, df = %d, p = %.4g\n", main_res$statistic, main_res$parameter, main_res$p.value))
        } else if (!is.null(main_res)) {
          print(main_res)
        } else {
          cat("Main results unavailable.\n")
        }

        cat("\n--- Effect Sizes ---\n")
        if (!is.null(stat_results$effect_size)) {
          print(as.data.frame(stat_results$effect_size), row.names = FALSE)
        } else cat("Not available for this analysis.\n")

        post_hoc <- stat_results$post_hoc_result
        if (!is.null(post_hoc)) {
          if (is.character(post_hoc)) {
            cat("\n--- Post-Hoc Notes ---\n", post_hoc, "\n")
          } else if (inherits(post_hoc, "data.frame") && nrow(post_hoc)) {
            cat("\n--- ", if (stat_results$test_name %in% c("t.test", "wilcox.test")) "Two-group comparison (unadjusted p-value)" else paste("Post-hoc results:", stat_results$post_hoc_method), " ---\n", sep = "")
            print(as.data.frame(post_hoc), row.names = FALSE)
            cat("*** p ≤ 0.001, ** p ≤ 0.01, * p ≤ 0.05\n")
          }
        }

        cat("\n--- Group Summaries ---\n")
        group_summary <- prep$data %>%
          dplyr::group_by(group) %>%
          dplyr::summarise(
            N = dplyr::n(),
            Mean = mean(.data[[prep$metric]], na.rm = TRUE),
            SD = stats::sd(.data[[prep$metric]], na.rm = TRUE),
            Median = stats::median(.data[[prep$metric]], na.rm = TRUE),
            Q1 = stats::quantile(.data[[prep$metric]], 0.25, na.rm = TRUE),
            Q3 = stats::quantile(.data[[prep$metric]], 0.75, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
        print(as.data.frame(group_summary), row.names = FALSE)
      }
    })
    paste(results_capture, collapse = "\n")
  })

  output$statsOutput <- renderPrint({
    validate(need(rv$grouped_data, "Assign groups in the popup or upload a completed group CSV."), need(input$statsMetric, "Select a metric."))
    withProgress(message = "Running statistics...", value = 0.5, {
      results_text <- stats_results()
      incProgress(0.5)
      cat(results_text)
    })
  })

  output$statsSummaryCards <- renderUI({
    if (is.null(rv$grouped_data) || !nrow(rv$grouped_data) || is.null(input$statsMetric)) return(NULL)
    result <- tryCatch(run_statistical_test(), error = function(e) NULL)
    prep <- tryCatch(stats_data_prep(), error = function(e) NULL)
    if (is.null(result) || is.null(prep)) {
      if (!identical(input$statsDesign, "independent")) return(div(class = "info-banner",
        if (identical(input$statsDesign, "paired")) "Paired or repeated measurements require a different analysis. This version supports independent groups only."
        else "Select the sample relationship to enable statistical analysis."))
      return(NULL)
    }
    if (!is.null(result$error)) {
      return(div(class = "warning-banner", icon("exclamation-triangle"), tags$strong("Statistical analysis error: "), result$error))
    }

    test_label <- result$test_label

    main_res <- result$main_result
    p_value <- NA_real_
    statistic_text <- "Available in details"
    if (!is.null(main_res)) {
      if (result$test_name == "aov") {
        group_row <- main_res %>% dplyr::filter(term == "group")
        if (nrow(group_row)) {
          p_value <- group_row$p.value[1]
          statistic_text <- sprintf("F = %.3f", group_row$statistic[1])
        }
      } else {
        p_value <- main_res$p.value[1] %||% NA_real_
        if ("statistic" %in% names(main_res)) statistic_text <- sprintf("Statistic = %.3f", main_res$statistic[1])
      }
    }
    p_text <- if (is.na(p_value)) "N/A" else if (p_value < 0.001) "<0.001" else sprintf("%.4f", p_value)
    p_class <- if (!is.na(p_value) && p_value <= 0.05) "status-pass" else "status-neutral"

    effect_text <- "Not available"
    effect_detail <- ""
    if (!is.null(result$effect_size) && inherits(result$effect_size, "data.frame") && nrow(result$effect_size)) {
      est_col <- intersect(c("estimate", "Eta2", "Epsilon2"), names(result$effect_size))[1]
      if (!is.na(est_col)) effect_text <- paste(ifelse(is.finite(result$effect_size[[est_col]]), sprintf("%.3f", result$effect_size[[est_col]]), "N/A"), collapse = " / ")
      effect_detail <- paste(result$effect_size$effect, collapse = " / ")
    }

    shapiro_p <- result$assumptions$shapiro_p
    levene_p <- result$assumptions$levene_p
    assumption_text <- result$assumptions$normality_status
    assumption_class <- if (identical(assumption_text, "No departure detected")) "status-neutral" else "status-warning"
    variance_detail <- if (is.na(levene_p)) "Variance test unavailable" else if (levene_p > 0.05) "No unequal variance detected" else "Unequal variance indicated"

    tagList(
      div(
        class = "metric-grid",
        metric_card("Executed test", test_label, paste(if (identical(input$statsTestType, "Auto")) "Selected by Auto" else "Selected manually",
                  "|", prep$num_groups, "groups;", nrow(prep$data), "observations")),
        metric_card("P-value", p_text, statistic_text, p_class),
        metric_card("Effect size", effect_text, effect_detail),
        metric_card("Normality screen", assumption_text, variance_detail, assumption_class)
      )
    )
  })

  output$statsGroupSummary <- DT::renderDataTable({
    prep <- stats_data_prep()
    req(prep)
    summary_df <- prep$data %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        N = dplyr::n(),
        Mean = mean(.data[[prep$metric]], na.rm = TRUE),
        SD = stats::sd(.data[[prep$metric]], na.rm = TRUE),
        Median = stats::median(.data[[prep$metric]], na.rm = TRUE),
        IQR = stats::IQR(.data[[prep$metric]], na.rm = TRUE),
        Minimum = min(.data[[prep$metric]], na.rm = TRUE),
        Maximum = max(.data[[prep$metric]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))
    DT::datatable(summary_df, rownames = FALSE, options = list(dom = "t", scrollX = TRUE))
  })

  stats_plot_object <- reactive({
    prep <- stats_data_prep()
    stat_results <- run_statistical_test()
    req(prep, stat_results)
    data_analysis <- prep$data
    metric <- prep$metric

    p <- ggplot(data_analysis, aes(x = group, y = .data[[metric]])) +
      geom_boxplot(aes(fill = group), outlier.shape = NA, alpha = 0.5, show.legend = FALSE) +
      geom_jitter(aes(color = group), width = 0.2, height = 0, alpha = 0.7, show.legend = FALSE) +
      scale_color_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
      scale_fill_viridis_d(option = "plasma", begin = 0.1, end = 0.9) +
      labs(x = "Group", y = format_col_name(metric)) +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(size = 14, colour = "grey15", angle = 0, vjust = 0.5, hjust = 0.5))

    pvals <- stat_results$post_hoc_result
    # Normalize only the plotting copy; exported two-group p-values stay unadjusted.
    if (inherits(pvals, "data.frame") && "p.value" %in% names(pvals)) {
      pvals <- dplyr::rename(pvals, p.adj = p.value, p.adj.signif = p.signif)
    }
    if (inherits(pvals, "data.frame") && nrow(pvals) && all(c("group1", "group2", "p.adj") %in% names(pvals))) {
      pvals <- pvals %>% dplyr::filter(!is.na(group1), !is.na(group2), !is.na(p.adj))
      # Filter before stacking brackets so hidden comparisons leave no empty space.
      # This only changes the plot; full statistical results retain every comparison.
      if (!isTRUE(input$statsShowNS)) {
        pvals <- pvals %>% dplyr::filter(p.adj <= 0.05)
      }
      if (nrow(pvals)) {
        max_y <- max(data_analysis[[metric]], na.rm = TRUE)
        min_y <- min(data_analysis[[metric]], na.rm = TRUE)
        range_y <- max_y - min_y
        if (!is.finite(range_y) || range_y == 0) range_y <- max(abs(max_y) * 0.1, 1e-6)
        step <- range_y * 0.1
        pvals <- pvals %>%
          dplyr::arrange(abs(match(group1, levels(data_analysis$group)) - match(group2, levels(data_analysis$group)))) %>%
          dplyr::mutate(y.position = max_y + seq_len(dplyr::n()) * step)
        p <- p + ggpubr::stat_pvalue_manual(
          pvals, label = "p.adj.signif", tip.length = 0.01,
          bracket.size = 0.5, size = 4, vjust = 0.5, hide.ns = FALSE
        ) +
          coord_cartesian(ylim = c(min_y - 0.05 * range_y, max(pvals$y.position) + 0.1 * range_y), clip = "off")
      } else {
        p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
      }
    } else {
      p <- p + scale_y_continuous(expand = expansion(mult = c(0.05, 0.08)))
    }
    p
  })

  output$statsPlot <- renderPlot({
    validate(need(rv$grouped_data, "Assign groups in the popup or upload a completed group CSV."), need(input$statsMetric, "Select a metric."))
    stats_plot_object()
  })

  # Size distributions
  observeEvent(parsed_distribution_data(), {
    dist_data <- parsed_distribution_data()
    if (!is.null(dist_data) && nrow(dist_data) > 0 && "BinCentre_nm" %in% names(dist_data)) {
      min_val <- floor(min(dist_data$BinCentre_nm, na.rm = TRUE))
      max_val <- ceiling(max(dist_data$BinCentre_nm, na.rm = TRUE))
      if (is.finite(min_val) && is.finite(max_val) && min_val < max_val) {
        updateSliderInput(session, "distBinRange", min = min_val, max = max_val, value = c(min_val, max_val))
      }
    }
  }, ignoreNULL = FALSE)

  observeEvent(input$distSelectAll, {
    dist_data <- parsed_distribution_data()
    req(dist_data)
    updateSelectizeInput(session, "distSampleSelect", selected = unique(dist_data$sample_key), server = FALSE)
  })

  observeEvent(input$distClearAll, {
    updateSelectizeInput(session, "distSampleSelect", selected = character(), server = FALSE)
    updateSelectInput(session, "distHighlight", selected = "")
  })

  observeEvent(input$distApplySource, {
    dist_data <- parsed_distribution_data()
    req(dist_data)
    if (!length(input$distSourceSelect)) {
      showNotification("Select at least one source file.", type = "warning")
      return()
    }
    source_samples <- dist_data %>%
      dplyr::filter(source_file %in% input$distSourceSelect) %>%
      dplyr::pull(sample_key) %>% unique()
    selected_samples <- unique(c(input$distSampleSelect %||% character(), source_samples))
    updateSelectizeInput(session, "distSampleSelect", selected = selected_samples, server = FALSE)
  })

  observe({
    if (isTRUE(input$distAverageByGroup) && (is.null(rv$grouped_data) || !nrow(rv$grouped_data))) {
      updateCheckboxInput(session, "distAverageByGroup", value = FALSE)
      showNotification("Assign statistical groups before averaging distributions by group.", type = "warning")
    }
  }) |> bindEvent(input$distAverageByGroup, rv$grouped_data)

  filtered_dist_data <- reactive({
    dist_data <- parsed_distribution_data()
    req(dist_data, nrow(dist_data) > 0)
    if (!length(input$distSampleSelect)) return(dist_data[0, , drop = FALSE])
    req(input$distBinRange)

    selected <- dist_data %>%
      dplyr::filter(
        sample_key %in% input$distSampleSelect,
        sample_key %in% filtered_data()$sample_key,
        BinCentre_nm >= input$distBinRange[1],
        BinCentre_nm <= input$distBinRange[2]
      )

    if (isTRUE(input$distAverageByGroup)) {
      req(rv$grouped_data)
      group_map <- rv$grouped_data %>% dplyr::distinct(source_file, sample_key, group)
      selected <- selected %>%
        dplyr::left_join(group_map, by = c("source_file", "sample_key")) %>%
        dplyr::filter(!is.na(group))
      validate(need(nrow(selected) > 0, "None of the selected distribution samples have statistical group assignments."))
      selected %>%
        dplyr::group_by(group, BinCentre_nm) %>%
        dplyr::summarise(
          mean_concentration = mean(ConcentrationAverage_dist, na.rm = TRUE),
          se_concentration = if (sum(!is.na(ConcentrationAverage_dist)) > 1) stats::sd(ConcentrationAverage_dist, na.rm = TRUE) / sqrt(sum(!is.na(ConcentrationAverage_dist))) else NA_real_,
          n_samples = dplyr::n_distinct(sample_key),
          .groups = "drop"
        ) %>%
        dplyr::transmute(
          group, BinCentre_nm,
          ConcentrationAverage_dist = mean_concentration,
          StandardError_dist = se_concentration,
          n_samples,
          series_key = paste0("Group: ", as.character(group)),
          source_file = "Grouped average"
        )
    } else {
      selected %>% dplyr::mutate(series_key = sample_key, n_samples = 1L)
    }
  })

  observe({
    plot_data <- tryCatch(filtered_dist_data(), error = function(e) NULL)
    if (is.null(plot_data) || !nrow(plot_data)) {
      updateSelectInput(session, "distHighlight", choices = c("None" = ""), selected = "")
      return()
    }
    series <- unique(as.character(plot_data$series_key))
    current <- isolate(input$distHighlight) %||% ""
    if (!current %in% series) current <- ""
    updateSelectInput(session, "distHighlight", choices = c("None" = "", stats::setNames(series, series)), selected = current)
  }) |> bindEvent(filtered_dist_data(), input$distAverageByGroup)

  output$distSelectionStatus <- renderUI({
    dist_data <- parsed_distribution_data()
    if (is.null(dist_data)) return(div(class = "info-banner", icon("info-circle"), "No size-distribution section was found in the uploaded data."))
    selected_n <- length(unique(input$distSampleSelect %||% character()))
    total_n <- dplyr::n_distinct(dist_data$sample_key)
    messages <- list(
      tags$span(class = "status-badge", paste(selected_n, "of", total_n, "samples selected"))
    )
    if (selected_n > 12 && !isTRUE(input$distAverageByGroup)) {
      messages <- c(messages, list(div(class = "warning-banner", style = "margin-top:8px;", icon("exclamation-triangle"), "Many individual curves are selected. Use source filtering, highlight one curve, or average by group for a clearer plot.")))
    }
    if (isTRUE(input$distAverageByGroup)) {
      messages <- c(messages, list(div(class = "info-banner", style = "margin-top:8px;", icon("users"), "Selected samples are averaged within their assigned statistical groups. Ribbons show the standard error across samples.")))
    }
    do.call(tagList, messages)
  })

  dist_legend_layout <- reactive({
    data <- tryCatch(filtered_dist_data(), error = function(e) NULL)
    series <- if (is.null(data) || !nrow(data) || !"series_key" %in% names(data)) character() else sort(unique(as.character(data$series_key)))
    highlight <- input$distHighlight %||% ""
    if (nzchar(highlight) && highlight %in% series) series <- highlight
    grouped <- isTRUE(input$distAverageByGroup)
    columns <- if (grouped && length(series) > 1L) 2L else 1L
    wrap_width <- if (grouped) 24L else 42L
    labels <- function(x) {
      if (grouped) x <- sub("^Group: ", "", x)
      vapply(x, function(value) {
        # Wrap words, then split long unbroken identifiers without truncating them.
        lines <- strsplit(stringr::str_wrap(value, width = wrap_width), "\n", fixed = TRUE)[[1]]
        lines <- unlist(lapply(lines, function(line) {
          starts <- seq.int(1L, max(1L, nchar(line)), by = wrap_width)
          substring(line, starts, starts + wrap_width - 1L)
        }), use.names = FALSE)
        paste(lines, collapse = "\n")
      }, character(1), USE.NAMES = FALSE)
    }
    wrapped <- labels(series)
    line_counts <- vapply(strsplit(wrapped, "\n", fixed = TRUE), length, integer(1))
    if (length(line_counts)) {
      # The guide fills rows in the same order; reserve space for the tallest label.
      padded <- c(line_counts, rep(0L, (-length(line_counts)) %% columns))
      row_lines <- apply(matrix(padded, ncol = columns, byrow = TRUE), 1, max)
      legend_height <- 40 + sum(row_lines * 24 + 12)
    } else legend_height <- 0
    list(columns = columns, labels = labels, extra_height = max(0, legend_height - 112))
  })

  dist_plot_height <- reactive({
    base_height <- switch(input$distPlotSize %||% "standard", compact = 400, large = 780, 600)
    base_height + dist_legend_layout()$extra_height
  })

  output$distributionPlotContainer <- renderUI({
    plotOutput("sizeDistributionPlot", height = paste0(dist_plot_height(), "px"))
  })

  size_dist_plot_object <- reactive({
    plot_data <- filtered_dist_data()
    req(plot_data, nrow(plot_data) > 0, input$distBinRange)

    use_log <- isTRUE(input$distLogY)
    has_nonpositive <- any(plot_data$ConcentrationAverage_dist <= 0, na.rm = TRUE)
    if (use_log && has_nonpositive) {
      showNotification("Log10 scale was not applied because the distribution contains zero or negative concentrations.", type = "warning")
      use_log <- FALSE
    }

    highlight <- input$distHighlight %||% ""
    highlight_active <- nzchar(highlight) && highlight %in% plot_data$series_key

    legend_layout <- dist_legend_layout()
    legend_labels <- legend_layout$labels

    if (highlight_active) {
      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = BinCentre_nm, y = ConcentrationAverage_dist, group = series_key)) +
        ggplot2::geom_line(color = "grey75", linewidth = 0.7, alpha = 0.7, na.rm = TRUE) +
        ggplot2::geom_line(
          data = plot_data %>% dplyr::filter(series_key == highlight),
          ggplot2::aes(color = series_key), linewidth = 1.35, na.rm = TRUE
        ) +
        viridis::scale_color_viridis(discrete = TRUE, option = "plasma", labels = legend_labels)
    } else {
      p <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = BinCentre_nm, y = ConcentrationAverage_dist, color = series_key, group = series_key)
      ) +
        ggplot2::geom_line(linewidth = 1, na.rm = TRUE) +
        viridis::scale_color_viridis(discrete = TRUE, option = "plasma", labels = legend_labels)
    }

    if (isTRUE(input$distShowError) && "StandardError_dist" %in% names(plot_data)) {
      ribbon_data <- plot_data %>%
        dplyr::mutate(
          valid_error = !is.na(StandardError_dist) & StandardError_dist >= 0,
          ymin_ribbon = ifelse(valid_error, ConcentrationAverage_dist - StandardError_dist, ConcentrationAverage_dist),
          ymax_ribbon = ifelse(valid_error, ConcentrationAverage_dist + StandardError_dist, ConcentrationAverage_dist),
          ymin_ribbon = if (use_log) ifelse(ymin_ribbon > 0, ymin_ribbon, NA_real_) else pmax(0, ymin_ribbon)
        )
      if (highlight_active) ribbon_data <- ribbon_data %>% dplyr::filter(series_key == highlight)
      p <- p +
        ggplot2::geom_ribbon(
          data = ribbon_data,
          ggplot2::aes(x = BinCentre_nm, ymin = ymin_ribbon, ymax = ymax_ribbon, fill = series_key, group = series_key),
          inherit.aes = FALSE, alpha = 0.18, linetype = "blank", show.legend = FALSE
        ) +
        viridis::scale_fill_viridis(discrete = TRUE, option = "plasma")
    }

    p <- p +
      ggplot2::labs(
        title = if (isTRUE(input$distAverageByGroup)) "Group-Averaged Particle Size Distributions" else "Particle Size Distributions",
        subtitle = if (highlight_active) paste("Highlighted:", highlight) else NULL,
        x = "Bin Centre (nm)",
        y = if (use_log) "Concentration (particles/mL, Log10)" else "Concentration (particles/mL)",
        color = if (isTRUE(input$distAverageByGroup)) "Group" else "Source | Sample | Well"
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        legend.position = "bottom",
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
        legend.text = ggplot2::element_text(size = 13, colour = "grey15"),
        legend.title = ggplot2::element_text(size = 13),
        legend.key.width = grid::unit(1.1, "lines")
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = legend_layout$columns, byrow = TRUE, title.position = "top")) +
      ggplot2::coord_cartesian(xlim = input$distBinRange)

    if (use_log) p <- p + ggplot2::scale_y_log10()
    p
  })

  output$sizeDistributionPlot <- renderPlot({
    req(parsed_distribution_data())
    validate(
      need(nrow(parsed_distribution_data()) > 0, "No size-distribution data were found."),
      need(length(input$distSampleSelect) > 0, "Select at least one sample.")
    )
    size_dist_plot_object()
  }, height = function() dist_plot_height())


  # --- QC Log Analysis ---
  observeEvent(input$qcLogUpload, {
    req(input$qcLogUpload$datapath)
    preview <- tryCatch(
      readr::read_csv(input$qcLogUpload$datapath, show_col_types = FALSE, col_types = cols(.default = "c")),
      error = function(e) NULL
    )
    if (is.null(preview) || !nrow(preview)) return()

    metric_choices_all <- c(
      "Mean size" = "MeanAverage",
      "Mode size" = "ModeAverage",
      "Concentration" = "ConcentrationAverage"
    )
    available_metric_choices <- metric_choices_all[unname(metric_choices_all) %in% names(preview)]

    if (!length(available_metric_choices)) {
      showNotification(
        "The QC log does not contain MeanAverage, ModeAverage or ConcentrationAverage.",
        type = "error", duration = 10
      )
      return()
    }

    selected_metric <- input$qcAnalysisMetric
    if (is.null(selected_metric) || !selected_metric %in% unname(available_metric_choices)) {
      selected_metric <- unname(available_metric_choices)[1]
    }
    updateSelectInput(
      session, "qcAnalysisMetric",
      choices = available_metric_choices,
      selected = selected_metric
    )

    if ("TimeCaptured" %in% names(preview)) {
      parsed_times <- parse_qc_datetime(preview$TimeCaptured)
      valid_dates <- as.Date(parsed_times[!is.na(parsed_times)], tz = DEFAULT_QC_TIMEZONE)
      if (length(valid_dates)) {
        date_min <- min(valid_dates)
        date_max <- max(valid_dates)
        updateDateRangeInput(
          session, "qcDateRange",
          start = date_min, end = date_max,
          min = date_min, max = date_max
        )
      }
    }
  }, ignoreInit = TRUE)

  qc_log_data_filtered <- reactive({
    req(input$qcLogUpload, input$qcDateRange)
    df <- tryCatch(
      readr::read_csv(input$qcLogUpload$datapath, show_col_types = FALSE, col_types = cols(.default = "c")),
      error = function(e) {
        showNotification(paste("QC log read error:", conditionMessage(e)), type = "error")
        NULL
      }
    )
    req(df)
    validate(need(nrow(df) > 0, "QC log is empty."))

    required_base <- c("TimeCaptured", "Status")
    missing_base <- setdiff(required_base, names(df))
    validate(need(!length(missing_base), paste("QC log is missing:", paste(missing_base, collapse = ", "))))

    supported_metrics <- c("MeanAverage", "ModeAverage", "ConcentrationAverage")
    available_metrics <- intersect(supported_metrics, names(df))
    validate(
      need(length(available_metrics) > 0, "QC log must contain at least one supported metric: MeanAverage, ModeAverage or ConcentrationAverage."),
      need(input$qcAnalysisMetric %in% available_metrics, "The selected QC metric is not present in this file.")
    )

    se_choice <- intersect(c("ConcentrationSE", "concentration_se", "concentration_standard_error_particles_ml"), names(df))[1] %||% NA_character_
    if (!is.na(se_choice) && se_choice != "ConcentrationSE") df <- dplyr::rename(df, ConcentrationSE = dplyr::all_of(se_choice))

    df <- df %>%
      dplyr::filter(!is.na(TimeCaptured), nzchar(trimws(TimeCaptured))) %>%
      dplyr::mutate(
        Timestamp_parsed = parse_qc_datetime(TimeCaptured)
      )
    validate(need(any(!is.na(df$Timestamp_parsed)), "The TimeCaptured values could not be parsed."))

    numeric_columns <- intersect(c("MeanAverage", "ModeAverage", "MeanSE", "ModeSE", "ConcentrationAverage", "ConcentrationSE"), names(df))
    df <- df %>%
      dplyr::mutate(
        dplyr::across(dplyr::all_of(numeric_columns), ~ readr::parse_number(as.character(.x))),
        Status = toupper(trimws(Status)),
        Status = ifelse(Status %in% c("PASS", "FAIL", "ERROR"), Status, "N/A"),
        Date = as.Date(Timestamp_parsed, tz = DEFAULT_QC_TIMEZONE)
      ) %>%
      dplyr::filter(
        !is.na(Timestamp_parsed), !is.na(Date),
        Date >= input$qcDateRange[1], Date <= input$qcDateRange[2]
      ) %>%
      dplyr::arrange(Timestamp_parsed) %>%
      dplyr::mutate(Sequence = dplyr::row_number())

    validate(need(nrow(df) > 0, "No QC data are present in the selected date range."))
    validate(need(input$qcAnalysisMetric %in% names(df), "The selected QC metric is missing."))
    df <- df %>%
      dplyr::filter(!is.na(.data[[input$qcAnalysisMetric]])) %>%
      dplyr::mutate(Sequence = dplyr::row_number())
    validate(need(nrow(df) > 0, "The selected QC metric has no numeric values."))
    df
  })

  qc_analysis_results <- reactive({
    df_filtered <- qc_log_data_filtered()
    req(df_filtered)
    metric_col <- input$qcAnalysisMetric
    validate(need(nrow(df_filtered) > 1, "At least two QC measurements are required."))

    daily_stats <- df_filtered %>%
      dplyr::group_by(Date) %>%
      dplyr::summarise(
        N = sum(!is.na(.data[[metric_col]])),
        Mean = mean(.data[[metric_col]], na.rm = TRUE),
        SD = stats::sd(.data[[metric_col]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        CV_Intra = ifelse(N > 1 & Mean != 0 & !is.na(SD), SD / Mean * 100, NA_real_),
        Intra_Day_Pass = ifelse(!is.na(CV_Intra), CV_Intra <= input$intraDayCVLimit, NA)
      )

    n_days <- dplyr::n_distinct(df_filtered$Date)
    point_mean <- mean(df_filtered[[metric_col]], na.rm = TRUE)
    point_sd <- stats::sd(df_filtered[[metric_col]], na.rm = TRUE)
    mean_of_daily_means <- mean(daily_stats$Mean, na.rm = TRUE)
    sd_of_means <- if (n_days > 1) stats::sd(daily_stats$Mean, na.rm = TRUE) else NA_real_
    cv_inter <- if (n_days > 1 && mean_of_daily_means != 0 && !is.na(sd_of_means)) sd_of_means / mean_of_daily_means * 100 else NA_real_
    inter_day_pass <- if (!is.na(cv_inter)) cv_inter <= input$interDayCVLimit else NA

    list(
      metric_selected = metric_col, n_total_points = nrow(df_filtered), n_days = n_days,
      overall_mean = point_mean, overall_sd = point_sd,
      mean_of_daily_means = mean_of_daily_means, sd_between_days = sd_of_means,
      cv_inter = cv_inter, inter_day_pass = inter_day_pass,
      inter_day_limit = input$interDayCVLimit, daily_stats = daily_stats,
      intra_day_limit = input$intraDayCVLimit
    )
  })

  output$qcHistoryCards <- renderUI({
    if (is.null(input$qcLogUpload)) {
      return(div(class = "nta-card empty-state", icon("history"), h4("Upload a QC log to review historical precision.")))
    }
    results <- tryCatch(qc_analysis_results(), error = function(e) NULL)
    if (is.null(results)) return(NULL)
    failing_days <- sum(results$daily_stats$Intra_Day_Pass == FALSE, na.rm = TRUE)
    inter_text <- if (is.na(results$cv_inter)) "N/A" else sprintf("%.2f%%", results$cv_inter)
    inter_class <- if (is.na(results$inter_day_pass)) "status-neutral" else if (isTRUE(results$inter_day_pass)) "status-pass" else "status-fail"
    mean_text <- if (is.finite(results$overall_mean)) format(signif(results$overall_mean, 5), scientific = abs(results$overall_mean) >= 1e5) else "N/A"
    tagList(
      div(
        class = "metric-grid",
        metric_card("QC records", results$n_total_points, paste(results$n_days, "days in selected range")),
        metric_card("Overall mean", mean_text, format_col_name(results$metric_selected)),
        metric_card("Inter-day CV", inter_text, paste0("Limit: ", results$inter_day_limit, "%"), inter_class),
        metric_card("Intra-day failures", failing_days, paste0("Limit: ", results$intra_day_limit, "%"), if (failing_days == 0) "status-pass" else "status-fail")
      )
    )
  })

  output$qcAnalysisSummary <- renderPrint({
    validate(need(input$qcLogUpload, "Upload a QC log."))
    results <- qc_analysis_results()
    req(results)
    cat("--- Precision Summary (", format(input$qcDateRange[1]), " to ", format(input$qcDateRange[2]), ") ---\n", sep = "")
    cat("Metric:", format_col_name(results$metric_selected), "\n")
    cat("Valid points:", results$n_total_points, "\nDays:", results$n_days, "\n")
    cat("Overall mean (all records):", sprintf("%.3f", results$overall_mean), "\n")
    cat("Overall SD (all records):", if (is.na(results$overall_sd)) "N/A" else sprintf("%.3f", results$overall_sd), "\n")
    cat("Mean of daily means:", sprintf("%.3f", results$mean_of_daily_means), "\n")
    cat("SD between daily means:", if (is.na(results$sd_between_days)) "N/A (<2 days)" else sprintf("%.3f", results$sd_between_days), "\n")
    cat("Inter-day CV:", if (is.na(results$cv_inter)) "N/A (<2 days)" else sprintf("%.2f%%", results$cv_inter), "\n")
    cat("Inter-day limit:", sprintf("%.2f%%", results$inter_day_limit), "\n")
    cat("Inter-day status:", if (is.na(results$inter_day_pass)) "N/A (<2 days)" else ifelse(results$inter_day_pass, "PASS", "FAIL"), "\n\n")

    valid_intra <- results$daily_stats$CV_Intra[!is.na(results$daily_stats$CV_Intra)]
    avg_intra <- if (length(valid_intra)) mean(valid_intra) else NA_real_
    cat("--- Intra-day Summary ---\n")
    cat("Average intra-day CV:", if (is.na(avg_intra)) "N/A (no day has repeated measurements)" else sprintf("%.2f%%", avg_intra), "\n")
    cat("Intra-day limit:", sprintf("%.2f%%", results$intra_day_limit), "\n")
    cat("Days failing intra-day limit:", sum(results$daily_stats$Intra_Day_Pass == FALSE, na.rm = TRUE), "\n")
  })

  output$qcDailyTable <- DT::renderDataTable({
    validate(need(input$qcLogUpload, "Upload a QC log."))
    results <- qc_analysis_results()
    req(results)
    display <- results$daily_stats %>%
      dplyr::mutate(
        dplyr::across(where(is.numeric), ~ round(.x, 3)),
        Intra_Day_Pass = dplyr::case_when(
          N <= 1 ~ "N/A (N ≤ 1)",
          is.na(CV_Intra) ~ "N/A",
          Intra_Day_Pass ~ "PASS",
          TRUE ~ "FAIL"
        ),
        Date = format(Date, "%Y-%m-%d")
      ) %>%
      dplyr::rename(`Intra CV (%)` = CV_Intra, `Intra Status` = Intra_Day_Pass)
    DT::datatable(display, options = list(pageLength = 5, searching = FALSE, lengthChange = FALSE), rownames = FALSE)
  })

  qc_history_plot_object <- reactive({
    log_df <- qc_log_data_filtered()
    metric <- input$qcAnalysisMetric
    req(log_df, metric)
    plot_df <- log_df %>% dplyr::filter(!is.na(.data[[metric]]))
    validate(need(nrow(plot_df) > 1, "At least two valid measurements are required."))

    if (isTRUE(input$useFixedQcBaseline)) {
      validate(
        need(is.finite(input$qcBaselineMean), "Enter a valid baseline mean."),
        need(is.finite(input$qcBaselineSD) && input$qcBaselineSD > 0, "Baseline SD must be greater than zero.")
      )
      centre <- input$qcBaselineMean
      spread <- input$qcBaselineSD
      baseline_note <- "Fixed baseline"
    } else {
      centre <- mean(plot_df[[metric]], na.rm = TRUE)
      spread <- stats::sd(plot_df[[metric]], na.rm = TRUE)
      validate(need(is.finite(centre) && is.finite(spread) && spread > 0, "The selected range cannot provide a valid mean and SD."))
      baseline_note <- "Estimated from selected date range"
    }

    limits <- c(centre - 3 * spread, centre - 2 * spread, centre, centre + 2 * spread, centre + 3 * spread)
    timestamps <- sort(unique(plot_df$Timestamp_parsed))
    short_range <- as.numeric(difftime(max(timestamps), min(timestamps), units = "days")) <= 2
    if (short_range) {
      tick_candidates <- timestamps
      date_format <- "%d %b %Y\n%H:%M"
    } else {
      dates <- sort(unique(as.Date(timestamps, tz = DEFAULT_QC_TIMEZONE)))
      # Calendar ticks must be spaced by elapsed time, not by record count.
      # Otherwise clusters of measurement dates produce overlapping labels.
      day_span <- as.integer(max(dates) - min(dates))
      tick_dates <- min(dates) + unique(round(seq(0, day_span, length.out = min(5, day_span + 1))))
      tick_candidates <- as.POSIXct(as.character(tick_dates), tz = DEFAULT_QC_TIMEZONE)
      date_format <- "%d %b\n%Y"
    }
    tick_indices <- unique(round(seq(1, length(tick_candidates), length.out = min(6, length(tick_candidates)))))
    ggplot(plot_df, aes(x = Timestamp_parsed, y = .data[[metric]])) +
      geom_hline(yintercept = limits, linetype = c("dashed", "dotted", "solid", "dotted", "dashed"), color = c("red", "orange", "blue", "orange", "red")) +
      geom_point(aes(color = Status), size = 2.5, alpha = 0.8) +
      geom_line(color = "grey50", alpha = 0.6) +
      scale_color_manual(values = c(PASS = "darkgreen", FAIL = "red", ERROR = "purple", `N/A` = "grey"), na.value = "grey") +
      scale_x_datetime(breaks = tick_candidates[tick_indices], date_labels = date_format,
        timezone = DEFAULT_QC_TIMEZONE,
        limits = range(c(timestamps, tick_candidates[tick_indices]))) +
      labs(
        title = paste("Levey-Jennings Plot:", format_col_name(metric)),
        subtitle = sprintf("Mean = %.2f, SD = %.2f - %s", centre, spread, baseline_note),
        x = paste0("Measurement date (", DEFAULT_QC_TIMEZONE, ")"), y = format_col_name(metric), color = "QC Status"
      ) +
      theme_bw(base_size = 12) + theme(legend.position = "bottom")
  })

  output$leveyJenningsPlot <- renderPlot({ qc_history_plot_object() })

  qc_history_download_name <- function(extension) {
    paste0("NTA_QC_history_", safe_file_stub(input$qcAnalysisMetric %||% "metric"), "_",
      format(as.Date(input$qcDateRange[1]), "%Y%m%d"), "_",
      format(as.Date(input$qcDateRange[2]), "%Y%m%d"), ".", extension)
  }

  output$downloadQcHistoryPNG <- downloadHandler(
    filename = function() qc_history_download_name("png"),
    contentType = "image/png",
    content = function(file) {
      p <- qc_history_plot_object()
      ggplot2::ggsave(file, plot = p, device = "png", width = 9, height = 6,
        units = "in", dpi = 300, bg = "white")
    }
  )

  output$downloadQcHistoryPDF <- downloadHandler(
    filename = function() qc_history_download_name("pdf"),
    contentType = "application/pdf",
    content = function(file) {
      p <- qc_history_plot_object()
      ggplot2::ggsave(file, plot = p, device = grDevices::cairo_pdf,
        width = 9, height = 6, units = "in", bg = "white")
    }
  )

  # Report download
  output$generateReport <- downloadHandler(
    filename = function(){ paste0("NTA_Report_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".html") },
    content = function(file){
      showNotification("Generating report...", type="message", duration=NULL, id="report-progress")
      on.exit(removeNotification("report-progress"))
      if (!file.exists(REPORT_TEMPLATE_FILENAME)) stop(paste0("Report template '", REPORT_TEMPLATE_FILENAME, "' not found in app working directory."))
      req(rv$parsed_result, message = "Upload or load example summary first.")
      active_id_col <- sampleIdCol_active(); req(active_id_col)
      params_list <- list(
        main_data = tryCatch(data_for_display(), error=function(e) NULL),
        application_version = APP_VERSION,
        active_filters = list(source_file = input$filterSourceFile, sample_name = input$filterSampleName, well_id = input$filterWellID),
        assignment_map = tryCatch({
          identities <- parsed_data_raw() %>% dplyr::distinct(source_file, sample_key, sample_name, well_id)
          if (!is.null(rv$group_map)) identities <- dplyr::left_join(identities, rv$group_map, by = c("source_file", "sample_key"))
          if (!"group" %in% names(identities)) identities$group <- NA_character_
          identities$included_by_filters <- identities$sample_key %in% filtered_data()$sample_key
          identities$included_in_statistics <- identities$included_by_filters & !is.na(identities$group) &
            identities$sample_key %in% tryCatch(stats_data_prep()$data$sample_key, error = function(e) character())
          identities
        }, error=function(e) NULL),
        comparison_plot1 = tryCatch(plot_object1(), error=function(e) NULL),
        comparison_plot2 = tryCatch(plot_object2(), error=function(e) NULL),
        metadata = tryCatch(parsed_metadata(), error=function(e) NULL),
        qc_status = tryCatch(qc_results(), error=function(e) NULL),
        group_assignments_file = if (!is.null(rv$group_map) && nrow(rv$group_map)) "Current assignments (popup or CSV)" else "None",
        stat_test_results_object = tryCatch(run_statistical_test(), error=function(e) NULL),
        stat_plot_object = tryCatch(stats_plot_object(), error=function(e) NULL),
        precision_results = tryCatch(qc_analysis_results(), error=function(e) NULL),
        selected_id_col = active_id_col,
        selected_stat_metric = input$statsMetric,
        report_run_time = Sys.time(),
        input_file_name = paste((rv$processed_files_info %||% rv$active_files_info)$name, collapse = ", ") %||% "N/A",
        distribution_data_for_report = tryCatch(filtered_dist_data(), error=function(e) NULL),
        dist_plot_object_for_report = tryCatch(size_dist_plot_object(), error=function(e) NULL),
        distribution_plot_height = dist_plot_height() / 90
      )
      rmarkdown::render(input = REPORT_TEMPLATE_FILENAME, output_format = "html_document", output_file = basename(file), output_dir = dirname(file), params = params_list, envir = new.env(parent = globalenv()), quiet = TRUE)
    }
  )
  
  # Contextual downloads
  output$downloadData <- downloadHandler(
    filename = function() paste0("parsed_nta_data_", Sys.Date(), ".csv"),
    content = function(file) {
      data_out <- data_for_display()
      req(data_out)
      readr::write_csv(data_out, file)
    }
  )
  output$downloadDataOverview <- downloadHandler(
    filename = function() paste0("filtered_nta_data_", Sys.Date(), ".csv"),
    content = function(file) {
      data_out <- tryCatch(filtered_data(), error = function(e) NULL)
      req(data_out)
      readr::write_csv(data_out, file)
    }
  )
  output$downloadPlot1 <- downloadHandler(
    filename = function() paste0("NTA_plot1_", safe_file_stub(input$metricSelect1 %||% "plot"), "_", input$comparisonPlotType %||% "points", "_", Sys.Date(), ".png"),
    content = function(file) {
      p <- plot_object1()
      req(p)
      p_dl <- p + ggplot2::labs(title = format_col_name(input$metricSelect1)) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5))
      ggplot2::ggsave(file, plot = p_dl, device = "png", width = 8.5, height = comparison_plot_height() / 90, dpi = 300, units = "in", limitsize = FALSE)
    }
  )
  output$downloadPlot2 <- downloadHandler(
    filename = function() paste0("NTA_plot2_", safe_file_stub(input$metricSelect2 %||% "plot"), "_", input$comparisonPlotType %||% "points", "_", Sys.Date(), ".png"),
    content = function(file) {
      p <- plot_object2()
      req(p)
      p_dl <- p + ggplot2::labs(title = format_col_name(input$metricSelect2)) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5))
      ggplot2::ggsave(file, plot = p_dl, device = "png", width = 8.5, height = comparison_plot_height() / 90, dpi = 300, units = "in", limitsize = FALSE)
    }
  )
  output$downloadStatsPlot <- downloadHandler(
    filename = function() paste0("NTA_stats_plot_", safe_file_stub(input$statsMetric %||% "stats"), "_", Sys.Date(), ".png"),
    content = function(file) {
      validate(need(rv$grouped_data, "Groups needed."), need(stats_plot_object(), "Statistics plot unavailable."))
      p <- stats_plot_object()
      req(p)
      p_dl <- p + ggplot2::labs(title = paste("Comparison of", format_col_name(input$statsMetric), "by Group")) + ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0.5))
      ggplot2::ggsave(file, plot = p_dl, device = "png", width = 7.5, height = 7, dpi = 300, units = "in", limitsize = FALSE)
    }
  )

  output$downloadStatsResults <- downloadHandler(
    filename = function() paste0("NTA_statistical_results_", Sys.Date(), ".txt"),
    content = function(file) writeLines(stats_results(), file, useBytes = TRUE)
  )
  output$downloadDistPlot <- downloadHandler(
    filename = function() paste0("NTA_distribution_plot_", Sys.Date(), ".png"),
    content = function(file) {
      p_dist <- size_dist_plot_object()
      req(p_dist)
      ggplot2::ggsave(file, plot = p_dist, device = "png", width = 10, height = dist_plot_height() / 90, dpi = 300, units = "in", limitsize = FALSE)
    }
  )

}

# --- Run App ---
shinyApp(ui = ui, server = server)
