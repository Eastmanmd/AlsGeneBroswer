# ============================================================
#  ALS Differential Gene Expression Explorer — Shiny App
#  ALS vs Control: fixed TPM plots, selectable DEG model
#  ALS Only: one selected analysis, one tissue or all tissues
#  No JavaScript — pure Shiny reactive tab switching
# ============================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(DT)

# ── colour palette ──────────────────────────────────────────
PAL_ALS    <- "#C0392B"
PAL_CTRL   <- "#2471A3"
PAL_C9_YES <- "#8E44AD"
PAL_C9_NO  <- "#17A589"
PAL_BG     <- "#F4F6F9"
PAL_PANEL  <- "#FFFFFF"
PAL_BORDER <- "#D5D8DC"
PAL_TEXT   <- "#1C2833"
PAL_MUTED  <- "#717D7E"
PAL_ACCENT <- "#2980B9"

# ── helpers ─────────────────────────────────────────────────
# Convert identifiers before concatenating factors with character choices.
clean_identifier <- function(x) {
  value <- trimws(gsub("\u00a0", " ", as.character(x), fixed = TRUE))
  value[!is.na(value) & !nzchar(value)] <- NA_character_
  value
}

# ── shared ggplot theme ──────────────────────────────────────
als_theme <- function() {
  theme_classic(base_family = "sans") +
    theme(
      plot.background   = element_rect(fill = PAL_PANEL, color = NA),
      panel.background  = element_rect(fill = PAL_PANEL, color = NA),
      panel.grid.minor  = element_blank(),
      axis.text         = element_text(color = "black",  size = 13),
      axis.title        = element_text(color = "black", size = 14),
      axis.title.y      = element_text(color = "black", margin = margin(r = 10)),
      axis.title.x      = element_text(color = "black", margin = margin(t = 8)),
      legend.position   = "top",
      legend.text       = element_text(color = "black", size = 13),
      legend.background = element_rect(fill = PAL_PANEL, color = NA),
      legend.key        = element_rect(fill = "transparent"),
      plot.margin       = margin(12, 12, 12, 12)
    )
}

ANALYSES <- list(
  c9orf72 = list(label = "C9orf72 status", column = "c9orf72", type = "categorical",
                 groups = c("C9 -", "C9 +"), colors = c("C9 -" = PAL_C9_NO, "C9 +" = PAL_C9_YES),
                 effect = "log2FC", effect_cutoff = 1, direction = "C9+ vs C9-; confirm against the fitted model"),
  sod1 = list(label = "SOD1 status", column = "sod1", type = "categorical",
              groups = c("SOD1 -", "SOD1 +"), colors = c("SOD1 -" = "#D68910", "SOD1 +" = "#117864"),
              effect = "log2FC", effect_cutoff = 1, direction = "SOD1+ vs SOD1-; confirm against the fitted model"),
  age_at_death = list(label = "Age at death", column = "age_at_death", type = "continuous",
                      axis = "Age at death (years)", effect = "log2FC", effect_cutoff = 0.0025, direction = "Log2 fold change for age at death"),
  disease_duration = list(label = "Disease duration", column = "disease_duration", type = "continuous",
                          axis = "Disease duration (months)", effect = "log2FC", effect_cutoff = 0.0025, direction = "Log2 fold change for disease duration")
)
MODEL_LABELS <- c(existing = "Standard model", music = "MuSiC proportion adjusted")
TISSUE_LABELS <- c("Cerebellum" = "CBL", "Lumbar Spinal Cord" = "LSC",
                   "Cervical Spinal Cord" = "CSC", "Frontal Cortex" = "FCX", "Motor Cortex" = "MCX")

load_optional_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) {
    warning("Could not load ", path, ": ", conditionMessage(e)); NULL
  })
}

# Choose readable facet dimensions from the available plot width.
scatter_layout <- function(width = 1200, tissues = 5L) {
  if (length(width) != 1 || !is.finite(width) || width <= 0) width <- 1200
  columns <- if (width < 650) 1L else if (width < 1100) 2L else 3L
  columns <- min(columns, max(1L, as.integer(tissues)))
  rows    <- ceiling(max(1L, as.integer(tissues)) / columns)
  list(columns = columns, height = 80L + 260L * rows)
}

expression_plot <- function(d, analysis = NULL, all_tissues = TRUE, shared_y = FALSE, facet_columns = 3L) {
  if (is.null(analysis)) {
    column <- "condition"; groups <- c("Control", "ALS")
    colors <- c(Control = PAL_CTRL, ALS = PAL_ALS)
    label <- "Condition"; categorical <- TRUE
  } else {
    spec <- ANALYSES[[analysis]]
    column <- spec$column; groups <- spec$groups; colors <- spec$colors
    label <- spec$label; categorical <- spec$type == "categorical"
  }
  validate(need(column %in% names(d), paste(label, "metadata is not available yet.")))
  # Base-R subsetting is faster than dplyr::filter for small slices
  ok <- !is.na(d[[column]]) & is.finite(d$tpm) & !is.na(d$tissue)
  d  <- d[ok, , drop = FALSE]
  if (categorical) d <- d[d[[column]] %in% groups, , drop = FALSE]
  validate(need(nrow(d) > 0, paste("No usable expression observations for", label)))
  d$value <- d[[column]]
  if (categorical) {
    d$value <- factor(d$value, levels = groups)
    p <- ggplot(d, aes(y = log2(tpm + 1), fill = value, color = value))
    if (all_tissues) {
      p <- p + aes(x = tissue) +
        geom_boxplot(alpha = .25, outlier.shape = NA, width = .55,
                     position = position_dodge(.7), color = "black") +
        geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7, seed = 42),
                   size = 1.2, alpha = .65)
    } else {
      p <- p + aes(x = value) +
        geom_boxplot(alpha = .25, outlier.shape = NA, width = .45, color = "black") +
        geom_point(position = position_jitter(width = .18, height = 0, seed = 42), size = 1.5, alpha = .7)
    }
    p <- p +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_color_manual(values = colors, drop = FALSE) +
      labs(x = if (all_tissues) NULL else label, fill = NULL, color = NULL)
    if (all_tissues) p <- p + scale_x_discrete(labels = function(x) {
      ifelse(x %in% names(TISSUE_LABELS), TISSUE_LABELS[x], x)
    })
  } else {
    validate(need(is.numeric(d$value), paste(label, "must be numeric.")))
    d <- d[is.finite(d$value), , drop = FALSE]
    validate(need(nrow(d) > 0, "No finite clinical measurements available."))
    p <- ggplot(d, aes(x = value, y = log2(tpm + 1))) +
      geom_point(color = PAL_ALS, size = 1.15, alpha = .45)
    # Eligibility check via base R (avoids group_by overhead)
    if (all_tissues) {
      tissue_chr  <- as.character(d$tissue)
      n_obs       <- ave(rep(1L, nrow(d)), tissue_chr, FUN = sum)
      n_distinct  <- ave(d$value, tissue_chr, FUN = function(v) length(unique(v)))
      eligible    <- d[n_obs >= 3 & n_distinct >= 2, , drop = FALSE]
    } else {
      eligible <- if (nrow(d) >= 3 && length(unique(d$value)) >= 2) d else d[0L, , drop = FALSE]
    }
    if (nrow(eligible)) {
      p <- p + geom_smooth(data = eligible, method = "lm", formula = y ~ x,
                           color = PAL_ACCENT, fill = PAL_ACCENT, alpha = .15, linewidth = .8)
    }
    if (all_tissues) {
      counts        <- table(as.character(d$tissue))
      facet_labels  <- setNames(paste0(names(counts), " \u00b7 n = ", as.integer(counts)), names(counts))
      p <- p + facet_wrap(~tissue, ncol = facet_columns,
                          scales = if (shared_y) "fixed" else "free_y",
                          labeller = as_labeller(facet_labels))
    }
    p <- p + labs(x = spec$axis)
  }
  p <- p + labs(y = "log\u2082(TPM + 1)") + als_theme()
  if (!categorical && all_tissues) {
    p <- p + theme(
      strip.background  = element_rect(fill = "#F4F6F9", color = NA),
      strip.text        = element_text(face = "bold", size = 11, color = PAL_TEXT, margin = margin(8, 5, 8, 5)),
      panel.spacing     = grid::unit(1.2, "lines"),
      panel.grid.major.y = element_line(color = "#EDF0F3", linewidth = .3),
      axis.text         = element_text(size = 10),
      axis.title        = element_text(size = 12))
  }
  p
}

cvc_tissue_choices <- function(...) {
  tissues <- sort(unique(c(...)))
  tissues <- clean_identifier(tissues)
  tissues <- tissues[!is.na(tissues)]
  c("All tissues" = "__all__", setNames(tissues, tissues))
}

# Pre-compute column aliases once at load time
.DEG_ALIASES <- list(
  gene   = c("symbol", "gene"),
  tissue = "tissue",
  effect = c("effect", "log2FoldChange", "logFC", "LFC", "coefficient"),
  pvalue = c("pvalue", "PValue", "pval", "P.Value"),
  fdr    = c("fdr", "padj", "FDR", "adj.P.Val", "p.adjust")
)

normalize_deg <- function(x) {
  x <- as.data.frame(x)
  cn <- names(x)
  # Find each column once, up front
  col_idx <- lapply(.DEG_ALIASES, function(aliases) {
    found <- aliases[aliases %in% cn]
    if (!length(found)) stop("Missing DEG column: ", paste(aliases, collapse = "/"))
    found[1]
  })
  out <- x[, unlist(col_idx), drop = FALSE]
  names(out) <- names(.DEG_ALIASES)
  out$gene   <- clean_identifier(out$gene)
  out$tissue <- clean_identifier(out$tissue)
  missing_id <- is.na(out$gene) | is.na(out$tissue)
  if (any(missing_id)) {
    warning(sum(missing_id), " DEG row(s) excluded because the gene or tissue identifier is missing.", call. = FALSE)
    out <- out[!missing_id, , drop = FALSE]
  }
  for (nm in c("effect", "pvalue", "fdr")) {
    old        <- out[[nm]]
    out[[nm]]  <- suppressWarnings(as.numeric(as.character(old)))
    if (any(!is.na(old) & is.na(out[[nm]]))) stop("Non-numeric DEG values: ", nm)
  }
  out <- unique(out)
  for (nm in c("pvalue", "fdr")) {
    if (any(out[[nm]] < 0 | out[[nm]] > 1, na.rm = TRUE)) stop(nm, " outside [0, 1]")
  }
  out
}

# Fully vectorized — no vapply loops
format_deg_number <- function(x) {
  result <- character(length(x))
  is_bad  <- is.na(x) | !is.finite(x)
  is_zero <- !is_bad & x == 0
  rest    <- !is_bad & !is_zero
  result[is_bad]  <- "NA"
  result[is_zero] <- "0"
  if (any(rest)) {
    vals           <- x[rest]
    result[rest]   <- format(signif(vals, 3), scientific = abs(vals) < 0.01, trim = TRUE)
  }
  result
}

style_deg_value <- function(value, metric, effect_cutoff = 1) {
  text <- format_deg_number(value)
  if (metric == "log2FC") {
    color <- ifelse(is.na(value), PAL_MUTED,
                    ifelse(value < -effect_cutoff, PAL_CTRL,
                           ifelse(value > effect_cutoff, PAL_ALS, PAL_MUTED)))
    paste0('<span style="color:', color, '">', text, "</span>")
  } else if (metric %in% c("p-value", "FDR")) {
    ifelse(!is.na(value) & value < .05,
           paste0('<span style="font-weight:600">', text, "</span>"), text)
  } else {
    text
  }
}

prepare_cvc_table <- function(df, mobile = FALSE, effect_metric = "log2FC", effect_cutoff = 1) {
  if (mobile) {
    result <- data.frame(Tissue = df$tissue, check.names = FALSE)
    result$Tissue <- ifelse(result$Tissue %in% names(TISSUE_LABELS),
                            TISSUE_LABELS[result$Tissue], result$Tissue)
    result[[effect_metric]] <- style_deg_value(df$effect, effect_metric, effect_cutoff)
    result[["p-value"]]     <- style_deg_value(df$pvalue, "p-value")
    result$FDR              <- style_deg_value(df$fdr, "FDR")
    return(result)
  }
  result <- data.frame(Metric = c(effect_metric, "p-value", "FDR"), check.names = FALSE)
  for (i in seq_len(nrow(df))) {
    result[[as.character(df$tissue[i])]] <- c(
      style_deg_value(df$effect[i], effect_metric, effect_cutoff),
      style_deg_value(df$pvalue[i], "p-value"),
      style_deg_value(df$fdr[i], "FDR"))
  }
  result
}

# ── Load data ONCE ──────────────────────────────────────
tpm_long_data <- tryCatch(
  suppressMessages(readRDS("./data/tpm_long_data_factor_final.rds")),
  error = function(e) { message("TPM load error: ", e$message); NULL }
)
deg_all_data <- tryCatch(
  suppressMessages(readRDS("./data/deg_all_data_subset_final.rds")),
  error = function(e) { message("DEG load error: ", e$message); NULL }
)
case_only_data <- tryCatch(
  suppressMessages(readRDS("./data/case_only_data_subset_final.rds")),
  error = function(e) { message("Case-only load error: ", e$message); list() }
)

deg_cell_adjusted_data <- load_optional_rds("./data/deg_all_data_celltype_final.rds")
sod1_results           <- load_optional_rds("./data/deg_case_only_sod1.rds")
if (!is.null(sod1_results)) case_only_data$sod1 <- sod1_results

sod1_metadata <- load_optional_rds("./data/sod1_metadata.rds")
if (!is.null(tpm_long_data) && !is.null(sod1_metadata)) {
  if (all(c("sample", "sod1") %in% names(sod1_metadata)) &&
      !anyNA(sod1_metadata$sample) && !anyDuplicated(sod1_metadata$sample)) {
    added <- sod1_metadata$sod1[match(tpm_long_data$sample, sod1_metadata$sample)]
    if (!"sod1" %in% names(tpm_long_data)) {
      tpm_long_data$sod1 <- added
    } else {
      missing <- is.na(tpm_long_data$sod1)
      tpm_long_data$sod1 <- as.character(tpm_long_data$sod1)
      tpm_long_data$sod1[missing] <- as.character(added[missing])
    }
  } else {
    warning("sod1_metadata.rds requires unique sample IDs and sample/sod1 columns")
  }
}

# ALS-only subset — computed once at startup
tpm_als_only <- if (!is.null(tpm_long_data) && "condition" %in% colnames(tpm_long_data))
  tpm_long_data[tpm_long_data$condition == "ALS", , drop = FALSE] else NULL

# Pre-compute gene list once for autocomplete
all_genes <- if (!is.null(tpm_long_data) && "gene" %in% colnames(tpm_long_data))
  sort(unique(tpm_long_data$gene)) else character(0)

# Pre-normalize DEG tables once at startup so reactives never call normalize_deg() again.
# Errors are surfaced as character messages; valid results are stored as data frames.
.norm_safe <- function(x, label = "") {
  if (is.null(x)) return(list(data = NULL, message = paste(label, "results not available yet")))
  tryCatch(list(data = normalize_deg(x), message = NULL),
           error = function(e) list(data = NULL, message = conditionMessage(e)))
}
.deg_norm        <- .norm_safe(deg_all_data, "Standard model")
.deg_music_norm  <- .norm_safe(deg_cell_adjusted_data, "MuSiC-adjusted model")
.case_norm       <- setNames(lapply(names(ANALYSES), function(nm) .norm_safe(case_only_data[[nm]], ANALYSES[[nm]]$label)),
                             names(ANALYSES))

# ── CSS ──────────────────────────────────────────────────────
app_css <- sprintf("
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body {
    background: %s; color: %s;
    font-family: 'Inter', sans-serif; font-size: 16px; min-height: 100vh;
  }
  .app-header {
    background: %s; border-bottom: 2px solid %s;
    padding: 16px 32px; display: flex; align-items: center; gap: 16px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.07);
  }
  .app-header h1 {
    font-family: 'Open Sans', monospace; font-size: 16px; font-weight: 700;
    letter-spacing: 0.04em; color: %s;
  }
  .app-header .subtitle { font-size: 12px; color: %s; margin-top: 2px; }
  .navbar { display: none; }
  .tab-bar {
    display: flex; background: %s; border-bottom: 2px solid %s;
    padding: 0 0 0 260px;
  }
  .tab-btn {
    background: none; border: none; padding: 13px 24px;
    font-family: 'Open Sans', monospace; font-size: 11px; font-weight: 700;
    letter-spacing: 0.08em; text-transform: uppercase; color: %s;
    cursor: pointer; border-bottom: 3px solid transparent;
    transition: color .15s, border-color .15s; position: relative; bottom: -2px;
  }
  .tab-btn.active { color: %s; border-bottom-color: %s; }
  .tab-btn:hover:not(.active) { color: %s; }
  .main-layout {
    display: grid; grid-template-columns: 260px 1fr;
    gap: 0; height: calc(100vh - 113px);
  }
  .sidebar {
    background: %s; border-right: 1px solid %s;
    padding: 24px 18px; overflow-y: auto;
  }
  .sidebar-section { margin-bottom: 26px; }
  .sidebar-label {
    font-family: 'Open Sans', monospace; font-size: 11px; font-weight: 700;
    letter-spacing: 0.12em; text-transform: uppercase; color: %s;
    margin-bottom: 10px; display: block;
  }
  .sidebar .selectize-input {
    background: %s !important; border: 1px solid %s !important;
    color: %s !important; border-radius: 6px !important;
    font-family: 'Open Sans', monospace !important;
    font-size: 15px !important; padding: 9px 12px !important;
    box-shadow: none !important; letter-spacing: 0.06em;
  }
  .sidebar .selectize-input.focus {
    border-color: %s !important;
    box-shadow: 0 0 0 2px rgba(41,128,185,0.15) !important;
  }
  .sidebar .selectize-dropdown {
    background: %s !important; border: 1px solid %s !important;
    color: %s !important; font-size: 14px !important;
    font-family: 'Open Sans', monospace !important;
  }
  .sidebar .selectize-dropdown-content .option { padding: 8px 12px; }
  .sidebar .selectize-dropdown-content .option:hover,
  .sidebar .selectize-dropdown-content .option.active {
    background: #EAF4FB !important; color: %s !important;
  }
  .sidebar .selectize-dropdown-content .option.selected {
    background: rgba(41,128,185,0.1) !important;
  }
  .tissue-select .selectize-input {
    background: %s !important; border: 1px solid %s !important;
    color: %s !important; border-radius: 6px !important;
    font-size: 14px !important; padding: 8px 12px !important;
    box-shadow: none !important;
  }
  .tissue-select .selectize-dropdown {
    background: %s !important; border: 1px solid %s !important;
    color: %s !important; font-size: 14px !important;
  }
  .tissue-select .selectize-dropdown-content .option:hover { background: #EAF4FB !important; }
  .gene-hint  { font-size: 13px; color: %s; margin-top: 6px; line-height: 1.4; }
  .gene-error { font-size: 13px; color: %s; margin-top: 4px; }
  .btn-search {
    width: 100%%;
    background: linear-gradient(135deg, %s 0%%, %s 100%%);
    border: none; border-radius: 6px; color: #fff;
    font-family: 'Open Sans', monospace; font-size: 13px; font-weight: 700;
    letter-spacing: 0.06em; padding: 10px 16px; cursor: pointer;
    transition: opacity .2s; margin-top: 8px;
  }
  .btn-search:hover { opacity: 0.85; }
  .legend-item {
    display: flex; align-items: center; gap: 8px;
    margin-bottom: 7px; font-size: 14px; color: %s;
  }
  .legend-dot { width: 10px; height: 10px; border-radius: 50%%; flex-shrink: 0; }
  .stat-chip {
    background: %s; border: 1px solid %s; border-radius: 20px;
    padding: 4px 12px; font-size: 13px; color: %s;
    display: inline-block; margin: 3px 3px 3px 0;
  }
  .main-content {
    display: flex; flex-direction: column; align-items: stretch;
    overflow-y: auto; background: %s; padding: 24px; gap: 24px;
  }
  .plot-panel {
    background: %s; border: 1px solid %s;
    border-radius: 10px; padding: 20px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05);
    max-width: 1600px; width: 100%%; margin: 0 auto;
  }
  .plot-panel-boxplot { max-height: 600px; overflow: hidden; width: 100%%; max-width: 1600px; margin: 0 auto;}
  .panel-title {
    font-family: 'Open Sans', monospace; font-size: 12px; font-weight: 700;
    letter-spacing: 0.1em; color: %s; margin-bottom: 16px;
  }
  .plots-grid-3 {
    display: grid; grid-template-columns: 1fr 1fr 1fr;
    gap: 20px; max-width: 1600px; width: 100%%; max-height: 600px; margin: 0 auto;
  }
  .plots-grid-3 > .plot-panel { max-width: 100%%; width: 100%%; margin: 0; min-width: 0; }
  .dataTables_wrapper { color: %s !important; font-size: 15px; max-width: 1600px;}
  table.dataTable thead th {
    background: %s !important; color: %s !important;
    border-bottom: 1px solid %s !important;
    font-family: 'Open Sans', monospace !important;
    font-size: 12px !important; letter-spacing: 0.08em !important;
    max-width: 1600px !important;
  }
  table.dataTable tbody tr { background: %s !important; }
  table.dataTable tbody tr:nth-child(even) { background: #F8F9FA !important; }
  table.dataTable tbody tr:hover { background: #EAF4FB !important; }
  table.dataTable tbody td { border-top: 1px solid %s !important; color: %s !important; }
  .dataTables_wrapper .dataTables_scroll,
  .dataTables_wrapper .dataTables_scrollBody,
  .dataTables_wrapper .dataTables_scrollHead {
    max-width: 1600px !important; overflow-x: auto !important;
  }
  .dataTables_filter input, .dataTables_length select {
    background: %s !important; border: 1px solid %s !important;
    color: %s !important; border-radius: 4px; padding: 4px 8px;
  }
  .dataTables_info, .dataTables_paginate { color: %s !important; }
  .paginate_button { color: %s !important; }
  .paginate_button.current { background: rgba(41,128,185,0.12) !important; border-radius: 4px; }
  .empty-state {
    display: flex; flex-direction: column;
    align-items: center; justify-content: center;
    height: 260px; gap: 12px; color: %s;
  }
  .empty-state .es-icon { font-size: 40px; opacity: 0.25; }
  .empty-state p { font-size: 15px; }
  ::-webkit-scrollbar { width: 5px; }
  ::-webkit-scrollbar-track { background: transparent; }
  ::-webkit-scrollbar-thumb { background: %s; border-radius: 10px; }
",
                   PAL_BG, PAL_TEXT,
                   PAL_PANEL, PAL_BORDER,
                   PAL_TEXT, PAL_MUTED,
                   PAL_PANEL, PAL_BORDER, PAL_MUTED, PAL_ACCENT, PAL_ACCENT, PAL_TEXT,
                   PAL_PANEL, PAL_BORDER, PAL_MUTED,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT, PAL_ACCENT,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT, PAL_TEXT,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT,
                   PAL_MUTED, PAL_ALS,
                   PAL_ACCENT, "#1A5276",
                   PAL_TEXT,
                   PAL_BG, PAL_BORDER, PAL_MUTED,
                   PAL_BG,
                   PAL_PANEL, PAL_BORDER,
                   PAL_MUTED,
                   PAL_TEXT,
                   PAL_BG, PAL_TEXT, PAL_BORDER,
                   PAL_PANEL,
                   PAL_BORDER, PAL_TEXT,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT,
                   PAL_MUTED, PAL_TEXT,
                   PAL_MUTED,
                   PAL_BORDER
)

responsive_css <- sprintf("
  .about-container {
    background: %s; min-height: calc(100vh - 113px);
    padding: 48px 64px; max-width: 1200px; margin: 0 auto;
  }
  @media (max-width: 768px) {
    .main-layout { display: flex; flex-direction: column; height: auto; }
    .tab-bar { padding: 0 16px; overflow-x: auto; white-space: nowrap; }
    .sidebar { border-right: none; border-bottom: 1px solid %s; max-height: none; }
    .plots-grid-3 { grid-template-columns: 1fr; max-height: none; }
    .plot-panel-boxplot { max-height: none; }
    .app-header { flex-direction: column; align-items: flex-start; gap: 8px; }
    .about-container { padding: 24px 16px !important; }
  }
", PAL_BG, PAL_BORDER)

app_css <- paste0(app_css, responsive_css, "
  .btn-download {
    background: transparent !important; border: 1px solid #D5D8DC !important;
    border-radius: 6px !important; color: #717D7E !important;
    font-family: 'Open Sans', monospace !important; font-size: 11px !important;
    font-weight: 700 !important; letter-spacing: 0.06em !important;
    padding: 6px 12px !important; cursor: pointer !important;
    transition: background-color .15s, border-color .15s, color .15s !important;
    display: inline-flex !important; align-items: center !important;
    gap: 6px !important; text-decoration: none !important;
    box-shadow: none !important; margin-top: 8px;
  }
  .btn-download:hover {
    background-color: #2980B9 !important; border-color: #2980B9 !important; color: #FFFFFF !important;
  }
  .checkbox { margin-top: 10px; margin-bottom: 10px; }
  .checkbox label {
    font-family: 'Open Sans', monospace !important; font-size: 12px !important;
    font-weight: 700 !important; letter-spacing: 0.08em !important;
    text-transform: uppercase !important; color: #717D7E !important; cursor: pointer !important;
  }
  .checkbox input[type='checkbox'] { margin-right: 8px !important; cursor: pointer !important; }
")

# ── UI ───────────────────────────────────────────────────────
ui <- fluidPage(
  title = "ALS Gene Expression Explorer",
  tags$head(
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Open+Sans:wght@300;400;500;600;700&display=swap",
      rel  = "stylesheet"
    ),
    tags$style(HTML(app_css)),
    tags$style(HTML("
      .cvc-plots {display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1fr);gap:20px;width:100%;max-width:1600px;margin:0 auto;}
      .cvc-plots.single {grid-template-columns:minmax(0,1fr);}
      .cvc-plots > .plot-panel {min-width:0;margin:0;}
      @media(max-width:1000px) {.cvc-plots {grid-template-columns:minmax(0,1fr);}}
    ")),
    tags$script(HTML("
      $(document).on('shiny:connected', function(event) {
        var isMobile = window.innerWidth <= 768;
        Shiny.setInputValue('is_mobile', isMobile);
        $(window).resize(function() {
          Shiny.setInputValue('is_mobile', window.innerWidth <= 768);
        });
      });
    "))
  ),
  tags$head(
    tags$script(async = NA, src = "https://www.googletagmanager.com/gtag/js?id=G-5CB8C7RLQ2"),
    tags$script(HTML("
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-5CB8C7RLQ2');
    "))
  ),
  div(class = "app-header",
      tags$a(
        href = "https://www.nygenome.org/science-technology/collaborative-research-programs/neurodegenerative-disease-research/als-consortium/",
        target = "_blank", rel = "noopener noreferrer",
        style = "display: inline-flex; flex-shrink: 0;",
        tags$img(src = "logo.svg", height = "80px", alt = "NYGC ALS Consortium",
                 style = "object-fit: contain; flex-shrink: 0;")
      ),
      div(
        tags$h1("ALS Gene Expression Explorer"),
        div(class = "subtitle", "Amyotrophic Lateral Sclerosis (ALS) Multi-Tissue DEG Analysis")
      )
  ),
  div(class = "tab-bar",
      uiOutput("tab_about_btn"),
      uiOutput("tab_cvc_btn"),
      uiOutput("tab_co_btn")
  ),
  conditionalPanel(
    condition = "input.active_tab == 'cvc'",
    div(class = "main-layout",
        div(class = "sidebar",
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Gene Search"),
                selectizeInput(
                  "gene_input", label = NULL, choices = NULL, selected = NULL,
                  options = list(placeholder = "e.g. SOD1", maxOptions = 50,
                                 searchField = "value", openOnFocus = FALSE,
                                 highlight = TRUE, closeAfterSelect = TRUE, create = FALSE),
                  width = "100%"
                ),
                div(class = "gene-hint", "Enter gene symbol (e.g. SOD1)"),
                tags$p(class = "gene-hint", "Selections update automatically.")
            ),
            div(class = "sidebar-section",
                selectInput("cvc_tissue", label = tags$span(class = "sidebar-label", "Tissue"),
                            choices = cvc_tissue_choices(tpm_long_data$tissue,
                                                         deg_all_data$tissue,
                                                         deg_cell_adjusted_data$tissue),
                            selected = "__all__", width = "100%")
            ),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Colour Legend"),
                div(class = "legend-item", div(class = "legend-dot", style = paste0("background:", PAL_CTRL)), "Control"),
                div(class = "legend-item", div(class = "legend-dot", style = paste0("background:", PAL_ALS)),  "ALS")
            ),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Significance Thresholds"),
                div(class = "stat-chip", "FDR < 0.05"),
                div(class = "stat-chip", "|log\u2082FC| > 1")
            ),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Status"),
                uiOutput("data_status")
            )
        ),
        div(class = "main-content",
            uiOutput("cvc_plots_ui"),
            div(class = "plot-panel",
                div(class = "panel-title", "Differential Expression Results by Tissue"),
                selectInput("deg_model", "DEG model",
                            choices = c("Standard model" = "existing",
                                        "Cell composition adjusted (MuSiC)" = "music"),
                            selected = "existing", width = "360px"),
                tags$p(class = "gene-hint",
                       "Model selection changes the DEG table and volcano plot. Boxplots show observed bulk TPM."),
                uiOutput("deg_model_notice"),
                uiOutput("deg_table_ui")
            )
        )
    )
  ),
  conditionalPanel(
    condition = "input.active_tab == 'co'",
    div(class = "main-layout",
        div(class = "sidebar",
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Gene"),
                selectizeInput(
                  "co_gene_input", label = NULL, choices = NULL, selected = NULL,
                  options = list(placeholder = "e.g. SOD1", maxOptions = 50,
                                 searchField = "value", openOnFocus = FALSE,
                                 highlight = TRUE, closeAfterSelect = TRUE, create = FALSE),
                  width = "100%"
                ),
                div(class = "gene-hint", "Enter gene symbol (e.g. SOD1)")
            ),
            div(class = "sidebar-section",
                selectInput("co_variable", label = tags$span(class = "sidebar-label", "Analysis"),
                            choices = c("C9orf72 Status" = "c9orf72", "SOD1 Status" = "sod1",
                                        "Age at Death" = "age_at_death", "Disease Duration" = "disease_duration"),
                            selected = "c9orf72", width = "100%"),
                selectInput("co_tissue", label = tags$span(class = "sidebar-label", "Tissue"),
                            choices = cvc_tissue_choices(
                              tpm_als_only$tissue,
                              unlist(lapply(case_only_data, function(x) as.character(x$tissue)), use.names = FALSE)),
                            selected = "__all__", width = "100%")
            ),
            conditionalPanel(
              condition = "input.co_tissue == '__all__' && (input.co_variable == 'age_at_death' || input.co_variable == 'disease_duration')",
              checkboxInput("co_shared_y", "Shared y-axis across tissues", value = FALSE)
            ),
            tags$p(class = "gene-hint", "Selections update automatically."),
            div(class = "sidebar-section", style = "margin-top:26px", uiOutput("co_legend")),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Significance Thresholds"),
                div(class = "stat-chip", "FDR < 0.05"),
                uiOutput("co_effect_threshold")
            )
        ),
        div(class = "main-content",
            uiOutput("co_plots_ui"),
            div(class = "plot-panel",
                div(class = "panel-title", textOutput("co_table_title")),
                DTOutput("co_deg_table")
            )
        )
    )
  ),
  conditionalPanel(
    condition = "input.active_tab == 'about'",
    div(class = "about-container",
        div(style = "margin-bottom: 40px;",
            div(class = "panel-title", "ABOUT THIS APP"),
            div(class = "plot-panel",
                tags$h3(style = paste0("font-family:'Open Sans',monospace; font-size:16px; font-weight:700; color:", PAL_TEXT, "; margin-bottom:12px;"),
                        "ALS Gene Expression Explorer"),
                tags$p(style = paste0("color:", PAL_TEXT, "; font-size:15px; line-height:1.75;"),
                       "This interactive browser enables exploration of differential gene expression
           across multiple tissues in Amyotrophic Lateral Sclerosis (ALS). Users can
           query any gene to visualize expression differences between ALS cases and
           controls, examine associations with clinical variables including C9orf72 and SOD1
           mutation status, age at death, and disease duration, and review full
           differential expression statistics across tissues."),
                tags$br(),
                tags$p(style = paste0("color:", PAL_TEXT, "; font-size:15px; line-height:1.75;"), "The app is organized into two views:"),
                tags$ul(style = paste0("color:", PAL_TEXT, "; font-size:15px; line-height:2; margin-left:20px;"),
                        tags$li(tags$b("ALS vs Control:"), " Boxplots and DEG statistics comparing ALS samples (bulk RNA-seq expression) to healthy controls across all available tissues."),
                        tags$li(tags$b("ALS Only Analysis:"), " Expression scatter plots and DEG results within ALS cases, stratified by C9orf72 or SOD1 mutation status, age at death, and disease duration, where available."))
            )
        ),
        div(style = "margin-bottom: 40px;",
            div(class = "panel-title", "ABSTRACT"),
            div(class = "plot-panel",
                tags$p(style = paste0("color:", PAL_TEXT, "; font-size:15px; line-height:1.75;"),
                       "Amyotrophic lateral sclerosis (ALS) is a devastating neurodegenerative disease with substantial genetic and clinical heterogeneity that impedes
                       therapeutic development. Large-scale multi-tissue genomic resources have transformed the study of neuropsychiatric and neurodegenerative diseases,
                       but no equivalent resource exists for ALS. Here we present the full NYGC ALS Consortium dataset, combining whole-genome sequencing from 4,746 donors
                       and bulk RNA-seq from 2,574 samples across 8 brain and spinal cord regions from 695 donors across the ALS disease spectrum. Our catalogue of small variants,
                       structural variants, and short tandem repeats identified likely pathogenic mutations in 15.6% of ALS cases. Gene expression and mRNA splicing analysis
                       across 5 major tissues reveals shared and region-specific features, highlighting microglial and T-cell dysregulation in the spinal cord. Mapping the
                       genetic regulation of expression and splicing across tissues identified associations with 6 ALS risk loci, whereas allele-specific rare variant analysis
                       detected expression effects for C9orf72 and OPTN. All data are immediately publicly available."),
                tags$br(),
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:14px; line-height:1.75;"), "Pre-print available at:"),
                tags$a(href = "https://www.medrxiv.org/content/10.64898/2026.04.29.26350889v1",
                       target = "_blank",
                       style = paste0("color:", PAL_ACCENT, "; font-size:14px; font-weight:600; text-decoration:none;"),
                       "\U0001f517  medRxiv"),
                tags$br(),
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:14px; line-height:1.75;"), "Code available at:"),
                tags$a(href = "https://github.com/Eastmanmd/AlsGeneBroswer/",
                       target = "_blank",
                       style = paste0("color:", PAL_ACCENT, "; font-size:14px; font-weight:600; text-decoration:none;"),
                       icon("github", class = "fa-lg"), " https://github.com/Eastmanmd/AlsGeneBroswer/"),
                tags$br(),
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:14px; line-height:1.75;"), "Full data available at:"),
                tags$a(href = "https://zenodo.org/records/20090702", target = "_blank",
                       style = paste0("color:", PAL_ACCENT, "; font-size:14px; font-weight:600; text-decoration:none;"), "\U0001f517  Metadata"),
                tags$a(href = "https://zenodo.org/records/18989237", target = "_blank",
                       style = paste0("color:", PAL_ACCENT, "; font-size:14px; font-weight:600; text-decoration:none;"), "\U0001f517  Gene Expression Data"),
                tags$a(href = "https://zenodo.org/records/19223231", target = "_blank",
                       style = paste0("color:", PAL_ACCENT, "; font-size:14px; font-weight:600; text-decoration:none;"), "\U0001f517  Splicing Data")
            )
        ),
        div(style = "margin-bottom: 40px;",
            div(class = "panel-title", "APP DEVELOPMENT"),
            div(class = "plot-panel",
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:12px; letter-spacing:0.08em; text-transform:uppercase; margin-bottom:16px;"), "Written & Designed By"),
                tags$a(href = "https://github.com/Eastmanmd/", target = "_blank", style = "text-decoration:none;",
                       div(style = paste0("display:flex; align-items:flex-start; gap:16px; padding:16px; background:", PAL_BG, "; border-radius:8px; border:1px solid ", PAL_BORDER, ";"),
                           div(style = paste0("width:44px; height:44px; border-radius:50%; flex-shrink:0; background:linear-gradient(135deg,", PAL_ALS, ",", PAL_C9_YES, "); display:flex; align-items:center; justify-content:center; color:#fff; font-weight:700; font-size:16px;"), "AO"),
                           div(
                             tags$p(style = paste0("font-weight:700; color:", PAL_TEXT, "; font-size:15px; margin-bottom:4px;"), "Ali Oku"),
                             tags$p(style = paste0("color:", PAL_MUTED, "; font-size:13px; line-height:1.6;"), "Senior Bioinformatics Analyst, Computational Biology"),
                             tags$p(style = paste0("color:", PAL_MUTED, "; font-size:13px;"), "New York Genome Center (NYGC)")
                           )
                       ))
            )
        ),
        tags$hr(style = paste0("border:none; border-top:1px solid ", PAL_BORDER, "; margin-bottom:16px;")),
        tags$p(style = paste0("color:", PAL_MUTED, "; font-size:12px; text-align:center;"),
               paste0("\u00a9 ", format(Sys.Date(), "%Y"), " New York Genome Center \u00b7 ALS Consortium"))
    )
  ),
  tags$div(style = "display:none",
           textInput("active_tab", label = NULL, value = "cvc"))
)

# ── Server ────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── Tab switching ──────────────────────────────────────────
  output$tab_cvc_btn <- renderUI({
    cls <- if (isTRUE(input$active_tab == "cvc")) "tab-btn active" else "tab-btn"
    actionButton("switch_cvc", "ALS vs Control", class = cls)
  })
  output$tab_co_btn <- renderUI({
    cls <- if (isTRUE(input$active_tab == "co")) "tab-btn active" else "tab-btn"
    actionButton("switch_co", "ALS Only Analysis", class = cls)
  })
  output$tab_about_btn <- renderUI({
    cls <- if (isTRUE(input$active_tab == "about")) "tab-btn active" else "tab-btn"
    actionButton("switch_about", "About", class = cls)
  })
  observeEvent(input$switch_cvc,   updateTextInput(session, "active_tab", value = "cvc"))
  observeEvent(input$switch_co,    updateTextInput(session, "active_tab", value = "co"))
  observeEvent(input$switch_about, updateTextInput(session, "active_tab", value = "about"))
  
  # Populate gene dropdowns server-side (avoids large JSON payload on load)
  updateSelectizeInput(session, "gene_input",    choices = all_genes,
                       selected = if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1], server = TRUE)
  updateSelectizeInput(session, "co_gene_input", choices = all_genes,
                       selected = if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1], server = TRUE)
  
  # ════════════════════════════════════════════════════════
  # TAB 2 — ALS vs Control
  # ════════════════════════════════════════════════════════
  selected_gene <- reactive({
    g <- trimws(input$gene_input)
    req(g != "", g %in% all_genes)
    g
  })
  
  output$data_status <- renderUI({
    items <- list()
    if (!is.null(tpm_long_data)) {
      items <- c(items, list(
        div(class = "stat-chip", paste0(length(unique(tpm_long_data$gene)),   " genes")),
        div(class = "stat-chip", paste0(length(unique(tpm_long_data$sample)), " samples")),
        div(class = "stat-chip", paste0(length(unique(tpm_long_data$tissue)), " tissues"))
      ))
    }
    if (!is.null(deg_all_data))
      items <- c(items, list(div(class = "stat-chip", paste0(length(unique(deg_all_data$tissue)), " DEG tissues"))))
    if (length(items) == 0)
      div(style = "color:#C0392B; font-size:12px;", "\u26a0 Data files not found")
    else tagList(items)
  })
  
  output$plot_title <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "") return("Enter a gene symbol to begin")
    if (is.null(input$cvc_tissue) || input$cvc_tissue == "__all__") paste0(g, " Expression Across Tissues")
    else paste0(g, " Expression \u2014 ", input$cvc_tissue)
  })
  
  # Gene-level TPM subset — recomputed only when gene or tissue changes
  cvc_gene_tpm <- reactive({
    req(selected_gene(), tpm_long_data)
    g <- selected_gene()
    d <- tpm_long_data[tpm_long_data$gene == g & !is.na(tpm_long_data$tissue) & !is.na(tpm_long_data$condition), , drop = FALSE]
    if (!is.null(input$cvc_tissue) && input$cvc_tissue != "__all__") {
      d <- d[clean_identifier(d$tissue) == clean_identifier(input$cvc_tissue), , drop = FALSE]
      validate(need(nrow(d) > 0, paste0("No expression observations for ", g, " in ", input$cvc_tissue, ".")))
    } else {
      validate(need(nrow(d) > 0, paste0("'", g, "' was not found in the expression data.")))
    }
    d
  })
  
  boxplot_plot <- reactive({
    d          <- cvc_gene_tpm()
    single_t   <- !is.null(input$cvc_tissue) && input$cvc_tissue != "__all__"
    if (single_t) return(expression_plot(d, all_tissues = FALSE))
    
    # Mobile: abbreviate tissue labels and flip axes
    if (isTRUE(input$is_mobile)) {
      tiss <- as.character(d$tissue)
      d$tissue <- ifelse(tiss %in% names(TISSUE_LABELS), TISSUE_LABELS[tiss], tiss)
    }
    p <- ggplot(d, aes(x = tissue, y = log2(tpm + 1), fill = condition, color = condition)) +
      geom_boxplot(alpha = 0.25, outlier.shape = NA, linewidth = 0.6,
                   width = 0.55, position = position_dodge(0.7), color = "black") +
      scale_fill_manual(values  = c("ALS" = PAL_ALS,  "Control" = PAL_CTRL)) +
      scale_color_manual(values = c("ALS" = PAL_ALS,  "Control" = PAL_CTRL)) +
      labs(x = NULL, y = "log\u2082(TPM + 1)", fill = NULL, color = NULL) +
      als_theme() +
      theme(axis.text.x = element_text(hjust = 0.5, size = 12))
    if (isTRUE(input$is_mobile)) {
      p <- p + coord_flip()
    } else {
      p <- p + geom_point(aes(group = condition),
                          position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.7),
                          size = 1.0, alpha = 0.75)
    }
    p
  })
  
  output$boxplot_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(tpm_long_data))
      return(div(class = "empty-state", div(class = "es-icon", "\U0001f9ec"), tags$p("Select a gene")))
    tagList(
      plotOutput("boxplot", height = "350px"),
      div(style = "text-align: right;", downloadButton("download_boxplot", "Download Plot", class = "btn-download"))
    )
  })
  output$boxplot <- renderPlot(boxplot_plot(), bg = PAL_PANEL)
  output$download_boxplot <- downloadHandler(
    filename = function() paste0(selected_gene(), "_",
                                 if (is.null(input$cvc_tissue) || input$cvc_tissue == "__all__") "all_tissues" else gsub(" ", "_", input$cvc_tissue),
                                 "_als_vs_control_boxplot.png"),
    content = function(file) ggsave(file, plot = boxplot_plot(), width = 12, height = 6, dpi = 300, device = "png"))
  
  # Use pre-normalized tables — no per-session normalize_deg() call
  selected_deg_norm <- reactive({
    if (identical(input$deg_model, "music")) .deg_music_norm else .deg_norm
  })
  
  output$deg_model_notice <- renderUI({
    if (identical(input$deg_model, "music") && is.null(.deg_music_norm$data))
      div(class = "gene-hint", style = "padding:12px;background:#FFF6E5;border-radius:6px;",
          "Cell composition adjusted results are not available yet. Add data/deg_als_vs_control_music.rds and restart.")
  })
  
  output$deg_table_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(selected_deg_norm()$data))
      return(div(class = "empty-state", div(class = "es-icon", "\U0001f4cb"),
                 tags$p("Select a gene and an available DEG model to see results")))
    DTOutput("deg_table")
  })
  
  cvc_model_results <- reactive({
    entry <- selected_deg_norm()
    validate(need(!is.null(entry$data), entry$message %||% "DEG data unavailable."))
    entry$data
  })
  
  cvc_gene_results <- reactive({
    g  <- selected_gene()
    df <- cvc_model_results()
    df <- df[df$gene == g, , drop = FALSE]
    if (!is.null(input$cvc_tissue) && input$cvc_tissue != "__all__")
      df <- df[df$tissue == input$cvc_tissue, , drop = FALSE]
    validate(need(nrow(df) > 0, "No DEG results for this gene/tissue selection."))
    validate(need(!anyDuplicated(df$tissue),
                  "Multiple distinct DEG rows match this gene symbol in one tissue."))
    df
  })
  
  output$deg_table <- renderDT({
    display <- prepare_cvc_table(cvc_gene_results(), mobile = isTRUE(input$is_mobile))
    datatable(display, rownames = FALSE, escape = 1, width = "100%",
              options = list(dom = "t", ordering = FALSE, scrollX = TRUE,
                             columnDefs = list(list(className = "dt-center", targets = seq_len(ncol(display) - 1)))),
              class = "display compact")
  })
  
  # Volcano — tissue-level, independent of gene filter
  volcano_data <- reactive({
    req(input$cvc_tissue, input$cvc_tissue != "__all__")
    df <- cvc_model_results()
    df <- df[df$tissue == input$cvc_tissue & is.finite(df$effect) &
               is.finite(df$fdr) & df$fdr >= 0 & df$fdr <= 1, , drop = FALSE]
    validate(need(nrow(df) > 0, "No usable DEG effects/FDR values for this tissue and model."))
    positive    <- df$fdr[df$fdr > 0]
    floor_value <- if (length(positive)) max(min(positive) / 10, .Machine$double.xmin) else .Machine$double.xmin
    df$minus_log10_fdr <- -log10(pmax(df$fdr, floor_value))
    df$status <- ifelse(df$fdr < .05 & df$effect >  1, "Higher in ALS",
                        ifelse(df$fdr < .05 & df$effect < -1, "Lower in ALS", "Below thresholds"))
    df
  })
  
  volcano_plot <- reactive({
    df     <- volcano_data()
    chosen <- df[df$gene == selected_gene(), , drop = FALSE]
    p <- ggplot(df, aes(x = effect, y = minus_log10_fdr, color = status)) +
      geom_point(alpha = .65, size = 1.5) +
      geom_vline(xintercept = c(-1, 1),   linetype = "dashed", color = PAL_MUTED, linewidth = .4) +
      geom_hline(yintercept = -log10(.05), linetype = "dashed", color = PAL_MUTED, linewidth = .4) +
      scale_color_manual(values = c("Higher in ALS" = PAL_ALS, "Lower in ALS" = PAL_CTRL, "Below thresholds" = "#B5B9BE")) +
      labs(x = "log2 fold change (ALS vs Control)", y = "-log10(FDR)", color = NULL,
           subtitle = paste(input$cvc_tissue, "\u2014",
                            if (identical(input$deg_model, "music")) "MuSiC adjusted model" else "Standard model")) +
      als_theme() + theme(legend.text = element_text(size = 10), legend.position = "bottom",
                          plot.subtitle = element_text(size = 11))
    if (nrow(chosen)) {
      p <- p +
        geom_point(data = chosen, shape = 21, fill = "#F4D03F", color = "black", size = 3, stroke = 1) +
        geom_text(data = chosen, aes(label = gene), color = PAL_TEXT, vjust = -1, size = 3.5, show.legend = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(.03, .15)))
    }
    p
  })
  
  output$volcano_plot <- renderPlot(volcano_plot(), bg = PAL_PANEL)
  
  output$cvc_plots_ui <- renderUI({
    single <- !is.null(input$cvc_tissue) && input$cvc_tissue != "__all__"
    box_panel <- div(class = "plot-panel",
                     div(class = "panel-title", uiOutput("plot_title")),
                     uiOutput("boxplot_ui"))
    if (!single) return(div(class = "cvc-plots single", box_panel))
    div(class = "cvc-plots", box_panel,
        div(class = "plot-panel",
            div(class = "panel-title", "Volcano Plot \u2014 ALS vs Control"),
            plotOutput("volcano_plot", height = "350px"),
            div(style = "text-align:right;",
                downloadButton("download_volcano", "Download Plot", class = "btn-download"))))
  })
  
  output$download_volcano <- downloadHandler(
    filename = function() paste0(gsub(" ", "_", input$cvc_tissue), "_",
                                 if (identical(input$deg_model, "music")) "music" else "standard", "_volcano.png"),
    content = function(file) ggsave(file, plot = volcano_plot(), width = 8, height = 6, dpi = 300, device = "png"))
  
  # ════════════════════════════════════════════════════════
  # TAB 3 — ALS Only Analysis
  # ════════════════════════════════════════════════════════
  co_selected <- reactive({
    req(input$co_tissue, input$co_variable)
    g <- trimws(input$co_gene_input)
    req(g != "", g %in% all_genes)
    list(gene = g, all_tissues = identical(input$co_tissue, "__all__"),
         tissue = input$co_tissue, variable = input$co_variable)
  })
  
  # Gene-level ALS-only TPM; cached per gene per session
  co_gene_data <- bindCache(reactive({
    g <- trimws(input$co_gene_input)
    req(g != "", g %in% all_genes)
    validate(need(!is.null(tpm_als_only), "ALS expression data is unavailable."))
    tpm_als_only[tpm_als_only$gene == g, , drop = FALSE]
  }), input$co_gene_input, cache = "session")
  
  co_gene_tissue_data <- reactive({
    req(input$co_tissue)
    d <- co_gene_data()
    if (identical(input$co_tissue, "__all__")) return(d)
    d[clean_identifier(d$tissue) == clean_identifier(input$co_tissue), , drop = FALSE]
  })
  
  output$co_legend <- renderUI({
    spec <- ANALYSES[[input$co_variable]]
    req(spec)
    if (spec$type != "categorical") return(NULL)
    tagList(tags$span(class = "sidebar-label", spec$label),
            lapply(spec$groups, function(group)
              div(class = "legend-item",
                  div(class = "legend-dot", style = paste0("background:", spec$colors[[group]])), group)))
  })
  
  co_scatter_layout <- reactive({
    sel  <- co_selected()
    spec <- ANALYSES[[sel$variable]]
    d    <- co_gene_tissue_data()
    if (spec$column %in% names(d)) {
      d <- d[!is.na(d[[spec$column]]) & is.finite(d$tpm) & !is.na(d$tissue), , drop = FALSE]
    }
    width <- session$clientData$output_co_selected_plot_width
    if (is.null(width) || !is.finite(width) || width <= 0)
      width <- if (isTRUE(input$is_mobile)) 500L else 1200L
    scatter_layout(width, length(unique(d$tissue)))
  })
  
  output$co_effect_threshold <- renderUI({
    req(input$co_variable)
    div(class = "stat-chip", paste0("|log\u2082FC| > ", ANALYSES[[input$co_variable]]$effect_cutoff))
  })
  
  output$co_table_title <- renderText({
    sel <- co_selected()
    paste(ANALYSES[[sel$variable]]$label, "\u2014 DEG results \u2014",
          if (isTRUE(sel$all_tissues)) "All tissues" else sel$tissue)
  })
  
  co_selected_plot <- reactive({
    sel <- co_selected()
    d   <- co_gene_tissue_data()
    validate(need(nrow(d) > 0, "No expression data for this gene/tissue selection."))
    if (sel$variable == "sod1" && "sod1" %in% names(d)) {
      value <- tolower(trimws(as.character(d$sod1)))
      value <- gsub("\u2212", "-", value, fixed = TRUE)
      d$sod1 <- ifelse(value %in% c("1","true","yes","+","positive","sod1+","sod1 +"), "SOD1 +",
                       ifelse(value %in% c("0","false","no","-","negative","sod1-","sod1 -"), "SOD1 -", NA_character_))
    }
    expression_plot(d, analysis = sel$variable, all_tissues = sel$all_tissues,
                    shared_y = isTRUE(input$co_shared_y), facet_columns = co_scatter_layout()$columns)
  })
  
  output$co_selected_plot <- bindCache(
    renderPlot(co_selected_plot(), bg = PAL_PANEL),
    input$co_gene_input, input$co_tissue, input$co_variable,
    isTRUE(input$co_shared_y), co_scatter_layout()$columns, cache = "session")
  
  output$co_selected_plot_title <- renderText({
    sel <- co_selected()
    paste(sel$gene, "\u2014", ANALYSES[[sel$variable]]$label, "\u2014",
          if (isTRUE(sel$all_tissues)) "Across All Tissues" else sel$tissue)
  })
  
  co_plot_height <- reactive({
    if (identical(input$co_tissue, "__all__") &&
        input$co_variable %in% c("age_at_death", "disease_duration"))
      paste0(co_scatter_layout()$height, "px") else "380px"
  })
  
  output$co_plots_ui <- bindCache(renderUI({
    all_t <- identical(input$co_tissue, "__all__")
    expression_panel <- div(class = "plot-panel",
                            div(class = "panel-title", textOutput("co_selected_plot_title")),
                            plotOutput("co_selected_plot", height = co_plot_height()),
                            div(style = "text-align:right;",
                                downloadButton("download_co_selected", "Download Plot", class = "btn-download")))
    if (all_t) return(div(class = "cvc-plots single", expression_panel))
    div(class = "cvc-plots", expression_panel,
        div(class = "plot-panel",
            div(class = "panel-title", "Volcano Plot \u2014 ALS Only"),
            plotOutput("co_volcano_plot", height = "380px"),
            div(style = "text-align:right;",
                downloadButton("download_co_volcano", "Download Plot", class = "btn-download"))))
  }), identical(input$co_tissue, "__all__"), co_plot_height(), cache = "session")
  
  output$download_co_selected <- downloadHandler(
    filename = function() {
      sel <- co_selected()
      paste0(sel$gene, "_", if (isTRUE(sel$all_tissues)) "all_tissues" else gsub(" ", "_", sel$tissue),
             "_", sel$variable, ".png")
    },
    content = function(file) {
      sel <- co_selected()
      ggsave(file, plot = co_selected_plot(), width = 12,
             height = if (isTRUE(sel$all_tissues) && ANALYSES[[sel$variable]]$type == "continuous")
               max(5, co_scatter_layout()$height / 90) else 6,
             dpi = 300, device = "png")
    })
  
  # Use pre-normalized case-only tables — no per-session normalize_deg() call
  co_analysis_results <- reactive({
    req(input$co_variable)
    entry <- .case_norm[[input$co_variable]]
    validate(need(!is.null(entry$data),
                  paste(ANALYSES[[input$co_variable]]$label, "\u2014", entry$message %||% "unavailable.")))
    entry$data
  })
  
  co_volcano_data <- reactive({
    req(input$co_tissue, input$co_tissue != "__all__", input$co_variable)
    df <- co_analysis_results()
    df <- df[df$tissue == clean_identifier(input$co_tissue) & is.finite(df$effect) &
               is.finite(df$fdr) & df$fdr >= 0 & df$fdr <= 1, , drop = FALSE]
    validate(need(nrow(df) > 0, "No usable DEG effects/FDR values for this analysis and tissue."))
    positive    <- df$fdr[df$fdr > 0]
    floor_value <- if (length(positive)) max(min(positive) / 10, .Machine$double.xmin) else .Machine$double.xmin
    df$minus_log10_fdr <- -log10(pmax(df$fdr, floor_value))
    cutoff <- ANALYSES[[input$co_variable]]$effect_cutoff
    sig    <- df$fdr < .05 & abs(df$effect) > cutoff
    df$status <- ifelse(sig & df$effect > 0, "Positive effect",
                        ifelse(sig & df$effect < 0, "Negative effect", "Below thresholds"))
    df
  })
  
  co_volcano_plot <- reactive({
    sel    <- co_selected()
    spec   <- ANALYSES[[sel$variable]]
    df     <- co_volcano_data()
    chosen <- df[df$gene == sel$gene, , drop = FALSE]
    p <- ggplot(df, aes(x = effect, y = minus_log10_fdr, color = status)) +
      geom_point(alpha = .65, size = 1.5) +
      geom_hline(yintercept = -log10(.05), linetype = "dashed", color = PAL_MUTED, linewidth = .4) +
      geom_vline(xintercept = c(-spec$effect_cutoff, spec$effect_cutoff), linetype = "dashed", color = PAL_MUTED, linewidth = .4) +
      scale_color_manual(values = c("Positive effect" = PAL_ALS, "Negative effect" = PAL_CTRL, "Below thresholds" = "#B5B9BE")) +
      labs(x = "log2 fold change", y = "-log10(FDR)", color = NULL,
           subtitle = paste(sel$tissue, "\u2014", spec$label)) +
      als_theme() + theme(legend.text = element_text(size = 10), legend.position = "bottom",
                          plot.subtitle = element_text(size = 11))
    if (nrow(chosen)) {
      p <- p +
        geom_point(data = chosen, shape = 21, fill = "#F4D03F", color = "black", size = 3, stroke = 1) +
        geom_text(data = chosen, aes(label = gene), color = PAL_TEXT, vjust = -1, size = 3.5, show.legend = FALSE) +
        scale_y_continuous(expand = expansion(mult = c(.03, .15)))
    }
    p
  })
  
  output$co_volcano_plot <- bindCache(
    renderPlot(co_volcano_plot(), bg = PAL_PANEL),
    input$co_gene_input, input$co_tissue, input$co_variable, cache = "session")
  
  output$download_co_volcano <- downloadHandler(
    filename = function() {
      sel <- co_selected()
      paste0(gsub(" ", "_", sel$tissue), "_", sel$variable, "_volcano.png")
    },
    content = function(file) ggsave(file, plot = co_volcano_plot(), width = 8, height = 6, dpi = 300, device = "png"))
  
  co_table_data <- reactive({
    sel <- co_selected()
    df  <- co_analysis_results()
    df  <- df[df$gene == sel$gene, , drop = FALSE]
    if (!isTRUE(sel$all_tissues))
      df <- df[df$tissue == clean_identifier(sel$tissue), , drop = FALSE]
    validate(need(nrow(df) > 0, "No results for this gene/tissue in the selected analysis."))
    validate(need(!anyDuplicated(df$tissue), "Multiple distinct DEG rows match this symbol in one tissue."))
    display <- prepare_cvc_table(df, mobile = isTRUE(input$is_mobile), effect_metric = "log2FC",
                                 effect_cutoff = ANALYSES[[sel$variable]]$effect_cutoff)
    if (!isTRUE(sel$all_tissues) && !isTRUE(input$is_mobile)) names(display)[2] <- "Value"
    display
  })
  
  output$co_deg_table <- renderDT({
    display         <- co_table_data()
    sel             <- co_selected()
    display_headers <- names(display)
    if (!isTRUE(sel$all_tissues) && !isTRUE(input$is_mobile)) display_headers[2] <- sel$tissue
    text_cols <- which(names(display) %in% c("Analysis", "Effect type", "Availability", "Metric", "Tissue"))
    datatable(display, colnames = display_headers, rownames = FALSE, escape = text_cols, width = "100%",
              options = list(dom = "t", ordering = FALSE, scrollX = TRUE, stateSave = FALSE,
                             columnDefs = list(list(className = "dt-center", targets = seq_len(ncol(display) - 1)))),
              class = "display compact")
  }, server = FALSE)
}

# Null-coalescing operator (base R does not have one)
`%||%` <- function(a, b) if (!is.null(a)) a else b

shinyApp(ui = ui, server = server)
