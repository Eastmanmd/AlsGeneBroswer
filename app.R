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

# load_deg_results <- function(deg_dir = "data/case_vs_control") {
#   files <- list.files(deg_dir, pattern = "_als_vs_control\\.csv$", full.names = TRUE)
#   if (length(files) == 0) { warning("No DEG files found in: ", deg_dir); return(NULL) }
#   purrr::map_dfr(files, function(f) {
#     tissue <- sub("_als_vs_control\\.csv$", "", basename(f))
#     df <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
#     if (!"symbol" %in% colnames(df)) df <- dplyr::rename(df, symbol = 1)
#     df$tissue <- tissue
#     df
#   })
# }

# load_case_only_results <- function(deg_dir = "data/case_only/") {
#   analyses <- list(
#     c9orf72          = "_c9orf72\\.csv$",
#     age_at_death     = "_age_at_death\\.csv$",
#     disease_duration = "_disease_duration\\.csv$"
#   )
#   results <- list()
#   for (nm in names(analyses)) {
#     files <- list.files(deg_dir, pattern = analyses[[nm]], full.names = TRUE)
#     if (length(files) == 0) next
#     results[[nm]] <- purrr::map_dfr(files, function(f) {
#       tissue <- sub(analyses[[nm]], "", basename(f))
#       df <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
#       df <- df[, !colnames(df) %in% c("baseMean", "lfcSE", "stat", "AveExpr", "z.std", "...1", "t")]
#       if (!"symbol" %in% colnames(df)) df <- dplyr::rename(df, symbol = 1)
#       df$tissue <- tissue
#       df
#     })
#   }
#   results
# }

# tidy_tpm <- function(tpm_path = "data/tpm.csv", meta_path = "data/metadata.csv") {
#   tpm  <- suppressMessages(readr::read_csv(tpm_path,  show_col_types = FALSE)) %>%
#     column_to_rownames(var = "...1")
#   meta <- suppressMessages(readr::read_csv(meta_path, show_col_types = FALSE)) %>%
#     column_to_rownames(var = "...1")
#   meta$sample <- rownames(meta)
#   tpm_long <- tpm %>%
#     rownames_to_column(var = "gene") %>%
#     tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "tpm")
#   keep_cols <- intersect(
#     c("sample", "Sample.Source", "Subject.Group",
#       "c9orf72", "Age.at.Death", "Disease.Duration.in.Months"),
#     colnames(meta)
#   )
#   meta_slim <- meta[, keep_cols]
#   colnames(meta_slim)[colnames(meta_slim) == "Sample.Source"]            <- "tissue"
#   colnames(meta_slim)[colnames(meta_slim) == "Subject.Group"]            <- "condition"
#   if ("C9orf72.Status" %in% colnames(meta_slim))
#     colnames(meta_slim)[colnames(meta_slim) == "C9orf72.Status"]         <- "c9orf72"
#   if ("Age.at.Death" %in% colnames(meta_slim))
#     colnames(meta_slim)[colnames(meta_slim) == "Age.at.Death"]           <- "age_at_death"
#   if ("Disease.Duration.in.Months" %in% colnames(meta_slim))
#     colnames(meta_slim)[colnames(meta_slim) == "Disease.Duration.in.Months"] <- "disease_duration"
#   dplyr::left_join(tpm_long, meta_slim, by = "sample")
# }

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
                 effect = "log2FC", direction = "C9+ vs C9-; confirm against the fitted model"),
  sod1 = list(label = "SOD1 status", column = "sod1", type = "categorical",
              groups = c("SOD1 -", "SOD1 +"), colors = c("SOD1 -" = "#D68910", "SOD1 +" = "#117864"),
              effect = "log2FC", direction = "SOD1+ vs SOD1-; confirm against the fitted model"),
  age_at_death = list(label = "Age at death", column = "age_at_death", type = "continuous",
                      axis = "Age at death (years)", effect = "Coefficient", direction = "Effect unit: verify the scaling used in the fitted model"),
  disease_duration = list(label = "Disease duration", column = "disease_duration", type = "continuous",
                          axis = "Disease duration (months)", effect = "Coefficient", direction = "Effect unit: verify the scaling used in the fitted model")
)
MODEL_LABELS <- c(existing = "Standard model", music = "MuSiC proportion adjusted")
TISSUE_LABELS <- c("Cerebellum" = "CBL", "Lumbar Spinal Cord" = "LSC",
                   "Cervical Spinal Cord" = "CSC", "Frontal Cortex" = "FCX", "Motor Cortex" = "MCX")


# Optional files can be added later without changing app code.
load_optional_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  tryCatch(readRDS(path), error = function(e) {
    warning("Could not load ", path, ": ", conditionMessage(e)); NULL
  })
}

expression_plot <- function(d, analysis = NULL, all_tissues = TRUE) {
  if (is.null(analysis)) {
    column <- 'condition'; groups <- c('Control', 'ALS')
    colors <- c(Control = PAL_CTRL, ALS = PAL_ALS)
    label <- 'Condition'; categorical <- TRUE
  } else {
    spec <- ANALYSES[[analysis]]
    column <- spec$column; groups <- spec$groups; colors <- spec$colors
    label <- spec$label; categorical <- spec$type == 'categorical'
  }
  validate(need(column %in% names(d), paste(label, 'metadata is not available yet.')))
  d <- d[!is.na(d[[column]]) & !is.na(d$tpm) & !is.na(d$tissue), , drop = FALSE]
  if (categorical) d <- d[d[[column]] %in% groups, , drop = FALSE]
  validate(need(nrow(d) > 0, paste('No usable expression observations for', label)))
  d$value <- d[[column]]
  if (categorical) {
    d$value <- factor(d$value, levels = groups)
    p <- ggplot(d, aes(y = log2(tpm + 1), fill = value, color = value))
    if (all_tissues) {
      p <- p + aes(x = tissue) +
        geom_boxplot(alpha = .25, outlier.shape = NA, width = .55,
                     position = position_dodge(.7), color = 'black') +
        geom_point(position = position_jitterdodge(jitter.width = .2, dodge.width = .7, seed = 42),
                   size = 1.2, alpha = .65)
    } else {
      p <- p + aes(x = value) +
        geom_boxplot(alpha = .25, outlier.shape = NA, width = .45, color = 'black') +
        geom_point(position = position_jitter(width = .18, height = 0, seed = 42), size = 1.5, alpha = .7)
    }
    p <- p + scale_fill_manual(values = colors, drop = FALSE) +
      scale_color_manual(values = colors, drop = FALSE) +
      labs(x = if (all_tissues) NULL else label, fill = NULL, color = NULL)
    if (all_tissues) p <- p + scale_x_discrete(labels = function(x) {
      ifelse(x %in% names(TISSUE_LABELS), TISSUE_LABELS[x], x)
    })
  } else {
    validate(need(is.numeric(d$value), paste(label, 'must be numeric.')))
    d <- d[is.finite(d$value), , drop = FALSE]
    validate(need(nrow(d) > 0, 'No finite clinical measurements available.'))
    p <- ggplot(d, aes(x = value, y = log2(tpm + 1))) +
      geom_point(color = PAL_ALS, size = 1.5, alpha = .7)
    # Fit only where a line is supported by the observed values.
    eligible <- d %>% group_by(tissue) %>%
      filter(n() >= 3, n_distinct(value) >= 2) %>% ungroup()
    if (nrow(eligible)) p <- p + geom_smooth(data = eligible, method = 'lm', formula = y ~ x,
                                             color = PAL_ACCENT, fill = PAL_ACCENT, alpha = .15, linewidth = .8)
    if (all_tissues) p <- p + facet_wrap(~tissue, ncol = 2)
    p <- p + labs(x = spec$axis)
  }
  p + labs(y = 'log\u2082(TPM + 1)') + als_theme()
}


# Convert identifiers before concatenating factors with character choices.
# Otherwise R can turn tissue factor labels into their integer level codes.
clean_identifier <- function(x) {
  value <- trimws(gsub("\u00a0", " ", as.character(x), fixed = TRUE))
  value[!is.na(value) & !nzchar(value)] <- NA_character_
  value
}
cvc_tissue_choices <- function(...) {
  tissues <- sort(unique(unlist(lapply(list(...), clean_identifier), use.names = FALSE)))
  tissues <- tissues[!is.na(tissues)]
  c("All tissues" = "__all__", setNames(tissues, tissues))
}

normalize_deg <- function(x) {
  x <- as.data.frame(x)
  aliases <- list(gene = c('symbol', 'gene'), tissue = 'tissue',
                  effect = c('effect', 'log2FoldChange', 'logFC', 'LFC', 'coefficient'),
                  pvalue = c('pvalue', 'PValue', 'pval', 'P.Value'),
                  fdr = c('fdr', 'padj', 'FDR', 'adj.P.Val', 'p.adjust'))
  out <- lapply(names(aliases), function(nm) {
    found <- intersect(aliases[[nm]], names(x))
    if (!length(found)) stop('Missing DEG column: ', nm)
    x[[found[1]]]
  })
  names(out) <- names(aliases)
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  out$gene <- clean_identifier(out$gene)
  out$tissue <- clean_identifier(out$tissue)
  # Unannotated DEG rows should not prevent annotated rows from being displayed.
  missing_id <- is.na(out$gene) | is.na(out$tissue)
  if (any(missing_id)) warning(sum(missing_id),
                               " DEG row(s) excluded because the gene or tissue identifier is missing.", call. = FALSE)
  out <- out[!missing_id, , drop = FALSE]
  for (nm in c('effect', 'pvalue', 'fdr')) {
    old <- out[[nm]]
    out[[nm]] <- suppressWarnings(as.numeric(as.character(old)))
    if (any(!is.na(old) & is.na(out[[nm]]))) stop('Non-numeric DEG values: ', nm)
  }
  # Identical rows can be repeated by input joins. Distinct rows with the same
  # symbol are retained for the volcano; ambiguous table matches are checked later.
  out <- unique(out)
  for (nm in c('pvalue', 'fdr')) {
    if (any(out[[nm]] < 0 | out[[nm]] > 1, na.rm = TRUE)) stop(nm, ' outside [0, 1]')
  }
  out
}



# Display only: keep unrounded numbers for threshold decisions and plotting.
format_deg_number <- function(x) {
  vapply(x, function(value) {
    if (is.na(value) || !is.finite(value)) return("NA")
    if (value == 0) return("0")
    format(signif(value, 3), scientific = abs(value) < 0.01, trim = TRUE)
  }, character(1), USE.NAMES = FALSE)
}
style_deg_value <- function(value, metric) {
  text <- format_deg_number(value)
  if (metric == "log2FC") {
    color <- ifelse(is.na(value), PAL_MUTED,
                    ifelse(value < -1, PAL_CTRL, ifelse(value > 1, PAL_ALS, PAL_MUTED)))
    paste0('<span style="color:', color, '">', text, '</span>')
  } else {
    ifelse(!is.na(value) & value < .05, paste0('<span style="font-weight:600">', text, '</span>'), text)
  }
}
prepare_cvc_table <- function(df, mobile = FALSE) {
  if (mobile) {
    result <- data.frame(Tissue = df$tissue, check.names = FALSE)
    result$Tissue <- ifelse(result$Tissue %in% names(TISSUE_LABELS), TISSUE_LABELS[result$Tissue], result$Tissue)
    result$log2FC <- style_deg_value(df$effect, "log2FC")
    result[["p-value"]] <- style_deg_value(df$pvalue, "p-value")
    result$FDR <- style_deg_value(df$fdr, "FDR")
    return(result)
  }
  result <- data.frame(Metric = c("log2FC", "p-value", "FDR"), check.names = FALSE)
  for (i in seq_len(nrow(df))) {
    result[[as.character(df$tissue[i])]] <- c(
      style_deg_value(df$effect[i], "log2FC"),
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

# Optional MuSiC-proportion-adjusted ALS vs Control results.
# Same columns as deg_all_data_subset_final.rds. No fallback to existing results.
deg_cell_adjusted_data <- load_optional_rds("./data/deg_als_vs_control_music.rds")
# Optional SOD1 DEG table; alternatively add a sod1 entry to the existing case-only list.
sod1_results <- load_optional_rds("./data/deg_case_only_sod1.rds")
if (!is.null(sod1_results)) case_only_data$sod1 <- sod1_results
# Optional sample-level SOD1 metadata; unknown status remains missing.
sod1_metadata <- load_optional_rds("./data/sod1_metadata.rds")
if (!is.null(tpm_long_data) && !is.null(sod1_metadata)) {
  if (all(c("sample", "sod1") %in% names(sod1_metadata)) &&
      !anyNA(sod1_metadata$sample) && !anyDuplicated(sod1_metadata$sample)) {
    added <- sod1_metadata$sod1[match(tpm_long_data$sample, sod1_metadata$sample)]
    if (!"sod1" %in% names(tpm_long_data)) tpm_long_data$sod1 <- added else {
      missing <- is.na(tpm_long_data$sod1)
      tpm_long_data$sod1 <- as.character(tpm_long_data$sod1)
      tpm_long_data$sod1[missing] <- as.character(added[missing])
    }
  } else warning("sod1_metadata.rds requires unique sample IDs and sample/sod1 columns")
}

# ALS-only subset for case-only plots
tpm_als_only <- if (!is.null(tpm_long_data) && "condition" %in% colnames(tpm_long_data))
  dplyr::filter(tpm_long_data, condition == "ALS") else NULL

# Pre-compute gene list once for autocomplete
all_genes <- if (!is.null(tpm_long_data) && "gene" %in% colnames(tpm_long_data))
  sort(unique(tpm_long_data$gene)) else character(0)

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

  /* ── Tab bar using Shiny navbarPage override ── */
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

  /* ── Layout ── */
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

  /* ── selectize gene search ── */
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

  /* tissue dropdown (case-only) */
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
    max-width: 1600px;
    width: 100%%;
    margin: 0 auto;
  }
  .plot-panel-boxplot { max-height: 600px; overflow: hidden; width: 100%%; max-width: 1600px; margin: 0 auto;}
  .panel-title {
    font-family: 'Open Sans', monospace; font-size: 12px; font-weight: 700;
    letter-spacing: 0.1em; color: %s; margin-bottom: 16px;
  }
  .plots-grid-3 {
    display: grid; grid-template-columns: 1fr 1fr 1fr; 
    gap: 20px; 
    max-width: 1600px; 
    width: 100%%;
    max-height: 600px;
    margin: 0 auto;
  }
  
  /* prevent individual panels inside the grid from overriding grid sizing */
  .plots-grid-3 > .plot-panel {
    max-width: 100%%;
    width: 100%%;
    margin: 0;
    min-width: 0;    /* critical — prevents grid blowout */
  }

  /* DataTables */
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
  
  /* The scroll body div DT injects when scrollX is enabled */
  .dataTables_wrapper .dataTables_scroll,
  .dataTables_wrapper .dataTables_scrollBody,
  .dataTables_wrapper .dataTables_scrollHead {
    max-width: 1600px !important;
    overflow-x: auto !important;
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
                   # tab bar
                   PAL_PANEL, PAL_BORDER, PAL_MUTED, PAL_ACCENT, PAL_ACCENT, PAL_TEXT,
                   # sidebar
                   PAL_PANEL, PAL_BORDER, PAL_MUTED,
                   # selectize gene
                   PAL_PANEL, PAL_BORDER, PAL_TEXT, PAL_ACCENT,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT, PAL_TEXT,
                   # tissue select
                   PAL_PANEL, PAL_BORDER, PAL_TEXT,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT,
                   # hints
                   PAL_MUTED, PAL_ALS,
                   # btn
                   PAL_ACCENT, "#1A5276",
                   PAL_TEXT,
                   # chips
                   PAL_BG, PAL_BORDER, PAL_MUTED,
                   # content
                   PAL_BG,
                   PAL_PANEL, PAL_BORDER,
                   PAL_MUTED,
                   # tables
                   PAL_TEXT,
                   PAL_BG, PAL_TEXT, PAL_BORDER,
                   PAL_PANEL,
                   PAL_BORDER, PAL_TEXT,
                   PAL_PANEL, PAL_BORDER, PAL_TEXT,
                   PAL_MUTED, PAL_TEXT,
                   # empty state / scrollbar
                   PAL_MUTED,
                   PAL_BORDER
)

responsive_css <- sprintf("
  .about-container {
    background: %s;
    min-height: calc(100vh - 113px);
    padding: 48px 64px;
    max-width: 1200px;
    margin: 0 auto;
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
    background: transparent !important;
    border: 1px solid #D5D8DC !important;
    border-radius: 6px !important;
    color: #717D7E !important;
    font-family: 'Open Sans', monospace !important;
    font-size: 11px !important;
    font-weight: 700 !important;
    letter-spacing: 0.06em !important;
    padding: 6px 12px !important;
    cursor: pointer !important;
    transition: background-color .15s, border-color .15s, color .15s !important;
    display: inline-flex !important;
    align-items: center !important;
    gap: 6px !important;
    text-decoration: none !important;
    box-shadow: none !important;
    margin-top: 8px;
  }
  .btn-download:hover {
    background-color: #2980B9 !important;
    border-color: #2980B9 !important;
    color: #FFFFFF !important;
  }
  .checkbox {
    margin-top: 10px;
    margin-bottom: 10px;
  }
  .checkbox label {
    font-family: 'Open Sans', monospace !important;
    font-size: 12px !important;
    font-weight: 700 !important;
    letter-spacing: 0.08em !important;
    text-transform: uppercase !important;
    color: #717D7E !important;
    cursor: pointer !important;
  }
  .checkbox input[type='checkbox'] {
    margin-right: 8px !important;
    cursor: pointer !important;
  }
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
  
  # ── Google Analytics  ────────────────────────────────────────────────
  tags$head(
    tags$script(async = NA, src = "https://www.googletagmanager.com/gtag/js?id=G-5CB8C7RLQ2"),
    tags$script(HTML("
      window.dataLayer = window.dataLayer || [];
      function gtag(){dataLayer.push(arguments);}
      gtag('js', new Date());
      gtag('config', 'G-5CB8C7RLQ2');
    "))
  ),
  
  # ── Header ────────────────────────────────────────────────
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
  
  # ── Pure-Shiny tab bar ────────────────────────────────────
  div(class = "tab-bar",
      uiOutput("tab_about_btn"), # About tab
      uiOutput("tab_cvc_btn"), # ALS vs Control
      uiOutput("tab_co_btn") # ALS only analysis
  ),
  
  # ── Active tab content (conditionalPanel — always rendered) ──
  # Both panels are kept in the DOM; only visibility toggled via
  # CSS display so plots don't need to re-render on tab switch.
  
  # TAB 1 — ALS vs Control
  conditionalPanel(
    condition = "input.active_tab == 'cvc'",
    div(class = "main-layout",
        div(class = "sidebar",
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Gene Search"),
                # selectizeInput with autocomplete (search-as-you-type)
                selectizeInput(
                  "gene_input",
                  label   = NULL,
                  choices = NULL, #all_genes,
                  selected = NULL, #if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1],
                  options = list(
                    placeholder      = "e.g. SOD1",
                    maxOptions       = 50,
                    searchField      = "value",
                    openOnFocus      = FALSE,
                    highlight        = TRUE,
                    closeAfterSelect = TRUE,
                    create           = FALSE
                    #plugins          = list("clear_button")
                  ),
                  width = "100%"
                ),
                div(class = "gene-hint",
                    "Enter gene symbol (e.g. SOD1 )"),
                tags$p(class = "gene-hint", "Selections update automatically.")
            ),
            div(class = "sidebar-section",
                selectInput("cvc_tissue", label = tags$span(class = "sidebar-label", "Tissue"),
                            choices = cvc_tissue_choices(tpm_long_data$tissue,
                                                         deg_all_data$tissue, deg_cell_adjusted_data$tissue),
                            selected = "__all__", width = "100%")
            ),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Colour Legend"),
                div(class = "legend-item",
                    div(class = "legend-dot", style = paste0("background:", PAL_CTRL)), "Control"),
                div(class = "legend-item",
                    div(class = "legend-dot", style = paste0("background:", PAL_ALS)), "ALS")
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
                            choices = c("Standard model" = "existing", "Cell composition adjusted (MuSiC)" = "music"),
                            selected = "existing", width = "360px"),
                tags$p(class = "gene-hint", "Model selection changes the DEG table and volcano plot. Boxplots show observed bulk TPM."),
                uiOutput("deg_model_notice"),
                uiOutput("deg_table_ui")
            )
        )
    )
  ),
  
  # TAB 2 — ALS Only Analysis
  conditionalPanel(
    condition = "input.active_tab == 'co'",
    div(class = "main-layout",
        div(class = "sidebar",
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Gene"),
                selectizeInput(
                  "co_gene_input",
                  label    = NULL,
                  choices  = NULL, #all_genes,
                  selected = NULL, #if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1],
                  options  = list(
                    placeholder      = "e.g. SOD1",
                    maxOptions       = 50,
                    searchField      = "value",
                    openOnFocus      = FALSE,
                    highlight        = TRUE,
                    closeAfterSelect = TRUE,
                    create           = FALSE
                    #plugins          = list("clear_button")
                  ),
                  width = "100%"
                ),
                div(class = "gene-hint",
                    "Enter gene symbol (e.g. SOD1 )")
            ),
            div(class = "sidebar-section",
                checkboxInput("co_all_tissues", "All Tissues", value = FALSE)
            ),
            div(class = "sidebar-section",
                selectInput("co_variable", label = tags$span(class = "sidebar-label", "Analysis"),
                            choices = c("C9orf72 Status" = "c9orf72", "SOD1 Status" = "sod1",
                                        "Age at Death" = "age_at_death", "Disease Duration" = "disease_duration"),
                            selected = "c9orf72", width = "100%"),
                conditionalPanel(condition = "!input.co_all_tissues",
                                 selectInput("co_tissue", label = tags$span(class = "sidebar-label", "Tissue"),
                                             choices = if (!is.null(tpm_als_only)) sort(unique(tpm_als_only$tissue)) else character(),
                                             selected = if (!is.null(tpm_als_only) && "Cerebellum" %in% tpm_als_only$tissue) "Cerebellum" else NULL,
                                             width = "100%"))
            ),
            tags$p(class = "gene-hint", "Selections update automatically."),
            div(class = "sidebar-section", style = "margin-top:26px", uiOutput("co_legend")),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Significance Thresholds"),
                div(class = "stat-chip", "FDR < 0.05"),
                div(class = "stat-chip", "|log\u2082FC| > 1")
            )
        ),
        div(class = "main-content",
            uiOutput("co_plots_ui"),
            div(class = "plot-panel",
                div(class = "panel-title", "Differential Expression Results Across Select Clinical Variables"),
                uiOutput("co_deg_table_ui")
            )
        )
    )
  ),
  
  # TAB 1 — About App Tab
  conditionalPanel(
    condition = "input.active_tab == 'about'",
    div(class = "about-container",
        
        # ── App Overview ─────────────────────────────────────
        div(style = "margin-bottom: 40px;",
            div(class = "panel-title", "ABOUT THIS APP"),
            div(class = "plot-panel",
                tags$h3(style = paste0("font-family:'Open Sans',monospace; font-size:16px;
                                font-weight:700; color:", PAL_TEXT, "; margin-bottom:12px;"),
                        "ALS Gene Expression Explorer"),
                tags$p(style = paste0("color:", PAL_TEXT, "; font-size:15px; line-height:1.75;"),
                       "This interactive browser enables exploration of differential gene expression
           across multiple tissues in Amyotrophic Lateral Sclerosis (ALS). Users can
           query any gene to visualize expression differences between ALS cases and
           controls, examine associations with clinical variables including C9orf72 and SOD1
           mutation status, age at death, and disease duration, and review full
           differential expression statistics across tissues."
                ),
                tags$br(),
                tags$p(style = paste0("color:", PAL_TEXT, "; font-size:15px; line-height:1.75;"),
                       "The app is organized into two  views:"
                ),
                tags$ul(style = paste0("color:", PAL_TEXT, "; font-size:15px;
                                line-height:2; margin-left:20px;"),
                        tags$li(tags$b("ALS vs Control:"),
                                " Boxplots and DEG statistics comparing ALS samples (bulk RNA-seq expression) to healthy controls
              across all available tissues."),
                        tags$li(tags$b("ALS Only Analysis:"),
                                " Expression scatter plots and DEG results within ALS cases, stratified
              by C9orf72 or SOD1 mutation status, age at death, and disease duration, where available.")
                )
            )
        ),
        
        # ── Study Goal ───────────────────────────────────────
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
                       detected expression effects for C9orf72 and OPTN. All data are immediately publicly available."
                ),
                tags$br(),
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:14px; line-height:1.75;"),
                       "Pre-print available at:"
                ),
                tags$a(
                  href   = "https://www.medrxiv.org/content/10.64898/2026.04.29.26350889v1",   
                  target = "_blank",
                  style  = paste0("color:", PAL_ACCENT, "; font-size:14px;
                           font-weight:600; text-decoration:none;"),
                  "\U0001f517  medRxiv"
                ),
                tags$br(),
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:14px; line-height:1.75;"),
                       "Code available at:"
                ),
                tags$a(
                  href   = "https://github.com/Eastmanmd/AlsGeneBroswer/",   
                  target = "_blank",
                  style  = paste0("color:", PAL_ACCENT, "; font-size:14px;
                           font-weight:600; text-decoration:none;"),
                  icon("github", class = "fa-lg"),
                  " https://github.com/Eastmanmd/AlsGeneBroswer/"
                ),
                tags$br(),
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:14px; line-height:1.75;"),
                       "Full data available at:"
                ),
                tags$a(
                  href   = "https://zenodo.org/records/20090702", 
                  target = "_blank",
                  style  = paste0("color:", PAL_ACCENT, "; font-size:14px;
                           font-weight:600; text-decoration:none;"),
                  "\U0001f517  Metadata"
                ),
                tags$a(
                  href   = "https://zenodo.org/records/18989237",   
                  target = "_blank",
                  style  = paste0("color:", PAL_ACCENT, "; font-size:14px;
                           font-weight:600; text-decoration:none;"),
                  "\U0001f517  Gene Expression Data"
                ),
                tags$a(
                  href   = "https://zenodo.org/records/19223231",   
                  target = "_blank",
                  style  = paste0("color:", PAL_ACCENT, "; font-size:14px;
                           font-weight:600; text-decoration:none;"),
                  "\U0001f517  Splicing Data"
                )
            )
        ),
        
        # ── Team ─────────────────────────────────────────────
        div(style = "margin-bottom: 40px;",
            div(class = "panel-title", "APP DEVELOPMENT"),
            div(class = "plot-panel",
                tags$p(style = paste0("color:", PAL_MUTED, "; font-size:12px;
                               letter-spacing:0.08em; text-transform:uppercase;
                               margin-bottom:16px;"),
                       "Written & Designed By"),
                # Developer card
                tags$a(
                  href = "https://github.com/Eastmanmd/", 
                  target = "_blank",  # opens in new tab
                  style = "text-decoration:none;",
                  
                  div(style = paste0("display:flex; align-items:flex-start; gap:16px;
                            padding:16px; background:", PAL_BG, ";
                            border-radius:8px; border:1px solid ", PAL_BORDER, ";"),
                      div(style = paste0("width:44px; height:44px; border-radius:50%; flex-shrink:0;
                              background:linear-gradient(135deg,", PAL_ALS, ",", PAL_C9_YES, ");
                              display:flex; align-items:center; justify-content:center;
                              color:#fff; font-weight:700; font-size:16px;"),
                          "AO"),  
                      div(
                        tags$p(style = paste0("font-weight:700; color:", PAL_TEXT, ";
                                   font-size:15px; margin-bottom:4px;"),
                               "Ali Oku"),
                        tags$p(style = paste0("color:", PAL_MUTED, "; font-size:13px;
                                   line-height:1.6;"),
                               "Bioinformatics Analyst, Computational Biology"),
                        tags$p(style = paste0("color:", PAL_MUTED, "; font-size:13px;"),
                               "New York Genome Center (NYGC)")
                      )
                  )
                ))
        ),
        
        # ── Footer ───────────────────────────────────────────
        tags$hr(style = paste0("border:none; border-top:1px solid ", PAL_BORDER, ";
                            margin-bottom:16px;")),
        tags$p(style = paste0("color:", PAL_MUTED, "; font-size:12px; text-align:center;"),
               paste0("\u00a9 ", format(Sys.Date(), "%Y"),
                      " New York Genome Center \u00b7 ALS Consortium"))
    )
  ),
  
  # Hidden input that holds the active tab — driven purely by Shiny
  # actionButtons above write to this via updateTextInput in server
  tags$div(style = "display:none",
           textInput("active_tab", label = NULL, value = "cvc"))
)

# ── Server ────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── Tab switching (no JS) ─────────────────────────────────
  # Two invisible action buttons; pressing either updates the
  # hidden text input that conditionalPanel reads.
  output$tab_cvc_btn <- renderUI({
    active <- isTRUE(input$active_tab == "cvc")
    cls <- if (active) "tab-btn active" else "tab-btn"
    actionButton("switch_cvc", "ALS vs Control", class = cls)
  })
  output$tab_co_btn <- renderUI({
    active <- isTRUE(input$active_tab == "co")
    cls <- if (active) "tab-btn active" else "tab-btn"
    actionButton("switch_co", "ALS Only Analysis", class = cls)
  })
  
  observeEvent(input$switch_cvc, {
    updateTextInput(session, "active_tab", value = "cvc")
  })
  observeEvent(input$switch_co, {
    updateTextInput(session, "active_tab", value = "co")
  })
  
  updateSelectizeInput(
    session, "gene_input",
    choices  = all_genes,
    selected = if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1],
    server   = TRUE
    #options=list(plugins=list("remove_button","clear_button"))
  )
  
  updateSelectizeInput(
    session, "co_gene_input",
    choices  = all_genes,
    selected = if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1],
    server   = TRUE
  )
  
  # ════════════════════════════════════════════════════════
  # TAB 1 — About the app
  # ════════════════════════════════════════════════════════
  output$tab_about_btn <- renderUI({
    active <- isTRUE(input$active_tab == "about")
    cls <- if (active) "tab-btn active" else "tab-btn"
    actionButton("switch_about", "About", class = cls)
  })
  
  observeEvent(input$switch_about, {
    updateTextInput(session, "active_tab", value = "about")
  })
  
  # ════════════════════════════════════════════════════════
  # TAB 2 — ALS vs Control
  # ════════════════════════════════════════════════════════
  
  # selected_gene <- eventReactive(input$go, {
  #   trimws(toupper(input$gene_input))
  # }, ignoreNULL = FALSE)
  
  selected_gene <- reactive({
    #g <- trimws(toupper(input$gene_input))
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
      items <- c(items, list(
        div(class = "stat-chip", paste0(length(unique(deg_all_data$tissue)), " DEG tissues"))
      ))
    if (length(items) == 0)
      div(style = "color:#C0392B; font-size:12px;", "\u26a0 Data files not found")
    else tagList(items)
  })
  
  output$plot_title <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "") return("Enter a gene symbol to begin")
    if (is.null(input$cvc_tissue) || input$cvc_tissue == "__all__") paste0(g, " Expression Across Tissues")
    else paste0(g, " Expression — ", input$cvc_tissue)
  })
  
  boxplot_plot <- reactive({
    req(selected_gene(), tpm_long_data)
    gene_data <- tpm_long_data %>%
      dplyr::filter(gene == selected_gene(), !is.na(tissue), !is.na(condition))
    if (!is.null(input$cvc_tissue) && input$cvc_tissue != "__all__") {
      gene_data <- dplyr::filter(gene_data, clean_identifier(tissue) == clean_identifier(input$cvc_tissue))
      validate(need(nrow(gene_data) > 0,
                    paste0("No expression observations for ", selected_gene(), " in ", input$cvc_tissue,
                           ". This tissue may have DEG results but no expression observations for the selected gene.")))
      return(expression_plot(gene_data, all_tissues = FALSE))
    }
    validate(need(nrow(gene_data) > 0,
                  paste0("'", selected_gene(), "' was not found in the expression data.")))
    
    if (isTRUE(input$is_mobile)) {
      t_map <- c("Cerebellum" = "CBL", 
                 "Lumbar Spinal Cord" = "LSC", 
                 "Cervical Spinal Cord" = "CSC", 
                 "Frontal Cortex" = "FCX", 
                 "Motor Cortex" = "MCX")
      if (is.factor(gene_data$tissue)) {
        lvl <- levels(gene_data$tissue)
        lvl <- ifelse(lvl %in% names(t_map), t_map[lvl], lvl)
        levels(gene_data$tissue) <- lvl
      } else {
        gene_data$tissue <- ifelse(gene_data$tissue %in% names(t_map), 
                                   t_map[gene_data$tissue], 
                                   gene_data$tissue)
      }
    }
    
    p <- ggplot(gene_data,
                aes(x = tissue, y = log2(tpm + 1), fill = condition, color = condition)) +
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
      p <- p + geom_jitter(aes(group = condition),
                           position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.7),
                           size = 1.0, alpha = 0.75)
    }
    p
  })
  
  output$boxplot_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(tpm_long_data))
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f9ec"),
                 tags$p("Select a gene")))
    tagList(
      plotOutput("boxplot", height = "350px"),
      div(style = "text-align: right;",
          downloadButton("download_boxplot", "Download Plot", class = "btn-download")
      )
    )
  })
  
  output$boxplot <- renderPlot({
    boxplot_plot()
  }, bg = PAL_PANEL)
  
  output$download_boxplot <- downloadHandler(
    filename = function() {
      paste0(selected_gene(), "_", if (is.null(input$cvc_tissue) || input$cvc_tissue == "__all__") "all_tissues" else gsub(" ", "_", input$cvc_tissue), "_als_vs_control_boxplot.png")
    },
    content = function(file) {
      ggsave(file, plot = boxplot_plot(), width = 12, height = 6, dpi = 300, device = "png")
    }
  )
  
  selected_deg_data <- reactive({
    if (identical(input$deg_model, "music")) deg_cell_adjusted_data else deg_all_data
  })
  output$deg_model_notice <- renderUI({
    if (identical(input$deg_model, "music") && is.null(deg_cell_adjusted_data))
      div(class = "gene-hint", style = "padding:12px;background:#FFF6E5;border-radius:6px;",
          "Cell composition adjusted results are not available yet. Add data/deg_als_vs_control_music.rds and restart. The TPM boxplots remain available.")
  })
  output$deg_table_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(selected_deg_data()))
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f4cb"),
                 tags$p("Select a gene and an available DEG model to see results")))
    DTOutput("deg_table")
  })
  
  cvc_model_results <- reactive({
    df <- selected_deg_data()
    validate(need(!is.null(df), "Selected model results are not available yet."))
    result <- tryCatch(normalize_deg(df), error = function(e) e)
    validate(need(!inherits(result, "error"),
                  if (inherits(result, "error")) conditionMessage(result) else ""))
    result
  })
  cvc_gene_results <- reactive({
    g <- selected_gene()
    df <- cvc_model_results()
    df <- df[df$gene == g, , drop = FALSE]
    if (!is.null(input$cvc_tissue) && input$cvc_tissue != "__all__")
      df <- df[df$tissue == input$cvc_tissue, , drop = FALSE]
    validate(need(nrow(df) > 0, "No DEG results for this gene/tissue selection."))
    validate(need(!anyDuplicated(df$tissue),
                  "Multiple distinct DEG rows match this gene symbol in one tissue. Resolve these input rows to a unique gene identifier before displaying the table."))
    df
  })
  output$deg_table <- renderDT({
    display <- prepare_cvc_table(cvc_gene_results(), mobile = isTRUE(input$is_mobile))
    datatable(display, rownames = FALSE, escape = 1, width = "100%",
              options = list(dom = "t", ordering = FALSE, scrollX = TRUE,
                             columnDefs = list(list(className = "dt-center", targets = seq_len(ncol(display) - 1)))),
              class = "display compact")
  })
  
  # Volcano uses every supplied DEG row for this tissue in the selected model.
  # It is deliberately independent of the selected-gene filter used by the table.
  volcano_data <- reactive({
    req(input$cvc_tissue, input$cvc_tissue != "__all__")
    df <- cvc_model_results()
    df <- df[df$tissue == input$cvc_tissue & is.finite(df$effect) &
               is.finite(df$fdr) & df$fdr >= 0 & df$fdr <= 1, , drop = FALSE]
    validate(need(nrow(df) > 0, "No usable DEG effects/FDR values for this tissue and model."))
    positive <- df$fdr[df$fdr > 0]
    floor_value <- if (length(positive)) max(min(positive) / 10, .Machine$double.xmin) else .Machine$double.xmin
    df$minus_log10_fdr <- -log10(pmax(df$fdr, floor_value))
    df$status <- ifelse(df$fdr < .05 & df$effect > 1, "Higher in ALS",
                        ifelse(df$fdr < .05 & df$effect < -1, "Lower in ALS", "Below thresholds"))
    df
  })
  volcano_plot <- reactive({
    df <- volcano_data()
    chosen <- df[df$gene == selected_gene(), , drop = FALSE]
    p <- ggplot(df, aes(x = effect, y = minus_log10_fdr, color = status)) +
      geom_point(alpha = .65, size = 1.5) +
      geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = PAL_MUTED, linewidth = .4) +
      geom_hline(yintercept = -log10(.05), linetype = "dashed", color = PAL_MUTED, linewidth = .4) +
      scale_color_manual(values = c("Higher in ALS" = PAL_ALS, "Lower in ALS" = PAL_CTRL, "Below thresholds" = "#B5B9BE")) +
      labs(x = "log2 fold change (ALS vs Control)", y = "-log10(FDR)", color = NULL,
           subtitle = paste(input$cvc_tissue, "—", if (identical(input$deg_model, "music")) "MuSiC adjusted model" else "Standard model")) +
      als_theme() + theme(legend.text = element_text(size = 10), legend.position = "bottom",
                          plot.subtitle = element_text(size = 11))
    if (nrow(chosen)) p <- p +
      geom_point(data = chosen, shape = 21, fill = "#F4D03F", color = "black", size = 3, stroke = 1) +
      geom_text(data = chosen, aes(label = gene), color = PAL_TEXT, vjust = -1, size = 3.5, show.legend = FALSE) +
      scale_y_continuous(expand = expansion(mult = c(.03, .15)))
    p
  })
  output$volcano_plot <- renderPlot(volcano_plot(), bg = PAL_PANEL)
  output$cvc_plots_ui <- renderUI({
    single <- !is.null(input$cvc_tissue) && input$cvc_tissue != "__all__"
    box_panel <- div(class = "plot-panel",
                     div(class = "panel-title", uiOutput("plot_title")), uiOutput("boxplot_ui"))
    if (!single) return(div(class = "cvc-plots single", box_panel))
    div(class = "cvc-plots", box_panel,
        div(class = "plot-panel",
            div(class = "panel-title", "Volcano Plot — ALS vs Control"),
            plotOutput("volcano_plot", height = "350px"),
            div(style = "text-align:right;", downloadButton("download_volcano", "Download Plot", class = "btn-download"))))
  })
  output$download_volcano <- downloadHandler(
    filename = function() paste0(gsub(" ", "_", input$cvc_tissue), "_", if (identical(input$deg_model, "music")) "music" else "standard", "_volcano.png"),
    content = function(file) ggsave(file, plot = volcano_plot(), width = 8, height = 6, dpi = 300, device = "png"))
  
  # ════════════════════════════════════════════════════════
  # TAB 3 — ALS Only Analysis
  # ════════════════════════════════════════════════════════
  
  # co_selected <- eventReactive(list(input$co_go, input$co_tissue), {
  #   req(input$co_tissue)
  #   list(gene = trimws(toupper(input$co_gene_input)), tissue = input$co_tissue)
  # }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  co_selected <- reactive({
    if (isTRUE(input$co_all_tissues)) {
      g <- trimws(input$co_gene_input)
      req(g != "", g %in% all_genes)
      req(input$co_variable)
      list(gene = g, all_tissues = TRUE, variable = input$co_variable)
    } else {
      req(input$co_tissue, input$co_variable)
      g <- trimws(input$co_gene_input)
      req(g != "", g %in% all_genes)
      list(gene = g, all_tissues = FALSE, tissue = input$co_tissue, variable = input$co_variable)
    }
  })
  
  # ALS data filtered to selected gene + tissue
  co_gene_tissue_data <- reactive({
    sel <- co_selected()
    req(sel$gene != "")
    validate(need(!is.null(tpm_als_only), "ALS expression data is unavailable."))
    if (isTRUE(sel$all_tissues)) {
      dplyr::filter(tpm_als_only, gene == sel$gene)
    } else {
      dplyr::filter(tpm_als_only, gene == sel$gene, tissue == sel$tissue)
    }
  })
  
  output$co_legend <- renderUI({
    spec <- ANALYSES[[input$co_variable]]
    req(spec)
    if (spec$type != "categorical") return(NULL)
    tagList(tags$span(class = "sidebar-label", spec$label),
            lapply(spec$groups, function(group) div(class = "legend-item",
                                                    div(class = "legend-dot", style = paste0("background:", spec$colors[[group]])), group)))
  })
  
  co_selected_plot <- reactive({
    sel <- co_selected()
    d <- co_gene_tissue_data()
    validate(need(nrow(d) > 0, "No expression data for this gene/tissue selection."))
    if (sel$variable == "sod1" && "sod1" %in% names(d)) {
      value <- tolower(trimws(as.character(d$sod1)))
      value <- gsub("\u2212", "-", value, fixed = TRUE)
      d$sod1 <- ifelse(value %in% c("1", "true", "yes", "+", "positive", "sod1+", "sod1 +"), "SOD1 +",
                       ifelse(value %in% c("0", "false", "no", "-", "negative", "sod1-", "sod1 -"), "SOD1 -", NA_character_))
    }
    expression_plot(d, analysis = sel$variable, all_tissues = sel$all_tissues)
  })
  output$co_selected_plot <- renderPlot(co_selected_plot(), bg = PAL_PANEL)
  output$co_selected_plot_title <- renderText({
    sel <- co_selected()
    paste(sel$gene, "\u2014", ANALYSES[[sel$variable]]$label, "\u2014",
          if (isTRUE(sel$all_tissues)) "Across All Tissues" else sel$tissue)
  })
  output$co_plots_ui <- renderUI({
    sel <- co_selected()
    div(class = "plot-panel",
        div(class = "panel-title", textOutput("co_selected_plot_title")),
        plotOutput("co_selected_plot", height = if (isTRUE(sel$all_tissues) && ANALYSES[[sel$variable]]$type == "continuous") "600px" else "380px"),
        tags$p(class = "gene-hint", "Observed bulk expression. Scatterplot lines are descriptive fits; DEG statistics come from the fitted analysis model."),
        div(style = "text-align:right;", downloadButton("download_co_selected", "Download Plot", class = "btn-download")))
  })
  output$download_co_selected <- downloadHandler(
    filename = function() {
      sel <- co_selected()
      paste0(sel$gene, "_", if (isTRUE(sel$all_tissues)) "all_tissues" else gsub(" ", "_", sel$tissue), "_", sel$variable, ".png")
    },
    content = function(file) {
      sel <- co_selected()
      ggsave(file, plot = co_selected_plot(), width = 12,
             height = if (isTRUE(sel$all_tissues) && ANALYSES[[sel$variable]]$type == "continuous") 8 else 6,
             dpi = 300, device = "png")
    })
  
  # ── ALS-only DEG table ──────────────────────────────────
  output$co_deg_table_ui <- renderUI({
    g <- trimws(input$co_gene_input)
    if (g == "" || length(case_only_data) == 0)
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f4cb"),
                 tags$p("Select a gene to see DEG results")))
    DTOutput("co_deg_table")
  })
  
  output$co_deg_table <- renderDT({
    sel <- co_selected()
    req(sel$gene != "", length(case_only_data) > 0)
    g <- sel$gene
    
    if (isTRUE(sel$all_tissues)) {
      # ── CASE 1: All Tissues ──
      df <- case_only_data[[sel$variable]]
      validate(need(!is.null(df), "No data for this variable."))
      df <- data.frame(df)
      
      tbl_base <- df %>%
        dplyr::filter(.data$symbol == .env$g) %>%
        dplyr::select(
          Tissue    = tissue,
          `log2FC`  = dplyr::any_of(c("log2FoldChange", "logFC", "LFC")),
          `p-value` = dplyr::any_of(c("pvalue", "PValue", "pval", "P.Value")),
          `FDR`     = dplyr::any_of(c("padj", "FDR", "adj.P.Val", "p.adjust"))
        ) %>%
        dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 3)))
      
      validate(need(nrow(tbl_base) > 0, paste0("'", g, "' not found in DEG results for ", sel$variable, ".")))
      
      if (isTRUE(input$is_mobile)) {
        t_map <- c("Cerebellum" = "CBL", 
                   "Lumbar Spinal Cord" = "LSC", 
                   "Cervical Spinal Cord" = "CSC", 
                   "Frontal Cortex" = "FCX", 
                   "Motor Cortex" = "MCX")
        
        tbl_display <- tbl_base %>%
          dplyr::mutate(Tissue = ifelse(as.character(Tissue) %in% names(t_map), 
                                        t_map[as.character(Tissue)], 
                                        as.character(Tissue))) %>%
          dplyr::mutate(
            log2FC = {
              val <- suppressWarnings(as.numeric(log2FC))
              dplyr::case_when(
                !is.na(val) & val < -1 ~ paste0('<span style="color:', PAL_CTRL, '">', val, '</span>'),
                !is.na(val) & val >  1 ~ paste0('<span style="color:', PAL_ALS,  '">', val, '</span>'),
                !is.na(val)            ~ paste0('<span style="color:', PAL_MUTED, '">', val, '</span>'),
                TRUE ~ as.character(val)
              )
            },
            `p-value` = {
              val <- suppressWarnings(as.numeric(`p-value`))
              dplyr::case_when(
                !is.na(val) & val < 0.05 ~ paste0('<span style="font-weight:600">', val, '</span>'),
                TRUE ~ as.character(val)
              )
            },
            FDR = {
              val <- suppressWarnings(as.numeric(FDR))
              dplyr::case_when(
                !is.na(val) & val < 0.05 ~ paste0('<span style="font-weight:600">', val, '</span>'),
                TRUE ~ as.character(val)
              )
            }
          ) %>%
          dplyr::rename(` ` = Tissue)
        
        col_targets <- seq_len(ncol(tbl_display) - 1)
        
      } else {
        tbl <- tbl_base %>%
          tidyr::pivot_longer(cols = -Tissue, names_to = "Metric", values_to = "Value") %>%
          tidyr::pivot_wider(names_from = Tissue, values_from = Value) %>%
          dplyr::rename(` ` = Metric)
        
        tissue_cols <- setdiff(colnames(tbl), " ")
        
        tbl_display <- tbl %>%
          dplyr::mutate(dplyr::across(dplyr::all_of(tissue_cols), ~ {
            val   <- suppressWarnings(as.numeric(.x))
            label <- tbl[[" "]][row_number()]
            dplyr::case_when(
              label == "log2FC" & !is.na(val) & val < -1 ~ paste0('<span style="color:', PAL_CTRL, '">', val, '</span>'),
              label == "log2FC" & !is.na(val) & val >  1 ~ paste0('<span style="color:', PAL_ALS,  '">', val, '</span>'),
              label == "log2FC" & !is.na(val) ~ paste0('<span style="color:', PAL_MUTED, '">', val, '</span>'),
              (label == "FDR" | label == "p-value") & !is.na(val) & val < 0.05 ~ paste0('<span style="font-weight:600">', val, '</span>'),
              TRUE ~ as.character(val)
            )
          }))
        
        col_targets <- seq_len(ncol(tbl_display) - 1)
      }
      
      datatable(tbl_display, rownames = FALSE, escape = FALSE, width = "100%",
                options = list(dom = "t", 
                               ordering = FALSE,
                               scrollX   = FALSE,
                               columnDefs = list(list(className = "dt-center", targets = col_targets))),
                class = "display compact")
      
    } else {
      # ── CASE 2: Single Tissue ──
      t <- sel$tissue
      
      pull_tbl <- function(nm, label) {
        df <- case_only_data[[nm]]
        if (is.null(df)) return(data.frame(Analysis = paste0(label, " (not available yet)"),
                                           log2FC = NA_real_, `p-value` = NA_real_, FDR = NA_real_, check.names = FALSE))
        df <- data.frame(df)
        
        df %>%
          dplyr::filter(.data$symbol == .env$g, .data$tissue == .env$t) %>%
          dplyr::select(
            `log2FC`  = dplyr::any_of(c("log2FoldChange", "logFC", "LFC")),
            `p-value` = dplyr::any_of(c("pvalue", "PValue", "pval", "P.Value")),
            `FDR`     = dplyr::any_of(c("padj", "FDR", "adj.P.Val", "p.adjust"))
          ) %>%
          dplyr::mutate(Analysis = label, .before = 1) %>%
          dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4)))
      }
      
      tbl <- dplyr::bind_rows(
        pull_tbl("c9orf72",          "C9orf72 Status"),
        pull_tbl("sod1",             "SOD1 Status"),
        pull_tbl("age_at_death",     "Age at Death"),
        pull_tbl("disease_duration", "Disease Duration")
      )
      
      validate(need(nrow(tbl) > 0,
                    paste0("'", g, "' not found in ALS-only DEG results for ", t, ".")))
      
      datatable(tbl, rownames = FALSE, width = "100%",
                options = list(dom = "t", ordering = FALSE,
                               scrollX   = FALSE,
                               columnDefs = list(list(className = "dt-center", targets = 1:3))),
                class = "display compact") %>%
        formatStyle("FDR",
                    color      = styleInterval(c(0.05), c("#1A7A4A", PAL_MUTED)),
                    fontWeight = styleInterval(c(0.05), c("600", "normal"))) %>%
        formatStyle("log2FC",
                    color = styleInterval(c(-1, 1), c(PAL_CTRL, PAL_MUTED, PAL_ALS)))
    }
  })
}

shinyApp(ui = ui, server = server)

