# ============================================================
#  ALS Differential Gene Expression Explorer — Shiny App
#  Two tabs: ALS vs Control | ALS Only Analysis
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

# ── Load data ONCE ──────────────────────────────────────
tpm_long_data <- tryCatch(
  suppressMessages(readRDS("/gpfs/commons/projects/ALS_Consortium_analysis/compbio/als_browser/als_shiny_browser/data/tpm_long_data_factor.rds")),
  error = function(e) { message("TPM load error: ", e$message); NULL }
)
deg_all_data <- tryCatch(
  suppressMessages(readRDS("/gpfs/commons/projects/ALS_Consortium_analysis/compbio/als_browser/als_shiny_browser/data/deg_all_data_subset.rds")),
  error = function(e) { message("DEG load error: ", e$message); NULL }
)
case_only_data <- tryCatch(
  suppressMessages(readRDS("/gpfs/commons/projects/ALS_Consortium_analysis/compbio/als_browser/als_shiny_browser/data/case_only_data_subset.rds")),
  error = function(e) { message("Case-only load error: ", e$message); list() }
)

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
    display: flex; flex-direction: column;
    overflow-y: auto; background: %s; padding: 24px; gap: 24px;
  }
  .plot-panel {
    background: %s; border: 1px solid %s;
    border-radius: 10px; padding: 20px;
    box-shadow: 0 1px 4px rgba(0,0,0,0.05);
  }
  .plot-panel-boxplot { max-height: 600px; overflow: hidden; }
  .panel-title {
    font-family: 'Open Sans', monospace; font-size: 12px; font-weight: 700;
    letter-spacing: 0.1em; color: %s; margin-bottom: 16px;
  }
  .plots-grid-3 {
    display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 20px; max-width: 1600px; max-height: 600px;
  }

  /* DataTables */
  .dataTables_wrapper { color: %s !important; font-size: 15px; }
  table.dataTable thead th {
    background: %s !important; color: %s !important;
    border-bottom: 1px solid %s !important;
    font-family: 'Open Sans', monospace !important;
    font-size: 12px !important; letter-spacing: 0.08em !important;
  }
  table.dataTable tbody tr { background: %s !important; }
  table.dataTable tbody tr:nth-child(even) { background: #F8F9FA !important; }
  table.dataTable tbody tr:hover { background: #EAF4FB !important; }
  table.dataTable tbody td { border-top: 1px solid %s !important; color: %s !important; }
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

# ── UI ───────────────────────────────────────────────────────
ui <- fluidPage(
  title = "ALS Gene Expression Explorer",
  tags$head(
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Open+Sans:wght@300;400;500;600;700&display=swap",
      rel  = "stylesheet"
    ),
    tags$style(HTML(app_css))
  ),
  
  # ── Header ────────────────────────────────────────────────
  div(class = "app-header",
      tags$img(src = "logo.svg", height = "80px",
               style = "object-fit: contain; flex-shrink: 0;"),
      div(
        tags$h1("ALS Gene Expression Explorer"),
        div(class = "subtitle", "Amyotrophic lateral sclerosis (ALS) Multi-tissue DEG Analysis")
      )
  ),
  
  # ── Pure-Shiny tab bar ────────────────────────────────────
  div(class = "tab-bar",
      uiOutput("tab_cvc_btn"),
      uiOutput("tab_co_btn")
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
                  choices = all_genes,
                  selected = if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1],
                  options = list(
                    placeholder      = "e.g. SOD1",
                    maxOptions       = 50,
                    searchField      = "value",
                    openOnFocus      = FALSE,
                    highlight        = TRUE,
                    closeAfterSelect = TRUE,
                    plugins          = list("clear_button")
                  ),
                  width = "100%"
                ),
                div(class = "gene-hint",
                    "Type to search gene symbols (e.g. SOD1, TARDBP, FUS)"),
                actionButton("go", "\u25b6  EXPLORE GENE",
                             class = "btn-search")
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
            div(class = "plot-panel plot-panel-boxplot",
                div(class = "panel-title", uiOutput("plot_title")),
                uiOutput("boxplot_ui")
            ),
            div(class = "plot-panel",
                div(class = "panel-title", "Differential Expression Results by Tissue"),
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
                  choices  = all_genes,
                  selected = if ("CHIT1" %in% all_genes) "CHIT1" else all_genes[1],
                  options  = list(
                    placeholder      = "e.g. SOD1",
                    maxOptions       = 50,
                    searchField      = "value",
                    openOnFocus      = FALSE,
                    highlight        = TRUE,
                    closeAfterSelect = TRUE,
                    plugins          = list("clear_button")
                  ),
                  width = "100%"
                ),
                div(class = "gene-hint",
                    "Type to search gene symbols (e.g. SOD1, TARDBP, FUS)")
            ),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Tissue"),
                div(class = "tissue-select",
                    uiOutput("co_tissue_select_ui"))
            ),
            actionButton("co_go", "\u25b6  EXPLORE GENE",
                         class = "btn-search"),
            div(class = "sidebar-section", style = "margin-top:26px",
                tags$span(class = "sidebar-label", "C9orf72 Legend"),
                div(class = "legend-item",
                    div(class = "legend-dot", style = paste0("background:", PAL_C9_YES)), "C9+"),
                div(class = "legend-item",
                    div(class = "legend-dot", style = paste0("background:", PAL_C9_NO)),  "C9\u2212")
            ),
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Significance Thresholds"),
                div(class = "stat-chip", "FDR < 0.05"),
                div(class = "stat-chip", "|log\u2082FC| > 1")
            )
        ),
        div(class = "main-content",
            div(class = "plots-grid-3",
                div(class = "plot-panel",
                    div(class = "panel-title", uiOutput("co_boxplot_title")),
                    uiOutput("co_c9_boxplot_ui")
                ),
                div(class = "plot-panel",
                    div(class = "panel-title", uiOutput("co_age_title")),
                    uiOutput("co_age_scatter_ui")
                ),
                div(class = "plot-panel",
                    div(class = "panel-title", uiOutput("co_dur_title")),
                    uiOutput("co_dur_scatter_ui")
                )
            ),
            div(class = "plot-panel",
                div(class = "panel-title", "Differential Expression Results Across Select Clinical Variables"),
                uiOutput("co_deg_table_ui")
            )
        )
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
  
  # ── Tissue dropdown (case-only) — default to Cerebellum ──
  output$co_tissue_select_ui <- renderUI({
    tissues <- if (!is.null(tpm_als_only)) sort(unique(tpm_als_only$tissue)) else character(0)
    default_tissue <- if ("Cerebellum" %in% tissues) "Cerebellum" else tissues[1]
    selectInput("co_tissue", label = NULL, choices = tissues,
                selected = default_tissue, width = "100%")
  })
  
  # ════════════════════════════════════════════════════════
  # TAB 1 — ALS vs Control
  # ════════════════════════════════════════════════════════
  
  selected_gene <- eventReactive(input$go, {
    trimws(toupper(input$gene_input))
  }, ignoreNULL = FALSE)
  
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
    paste0(g, " Expression Across Tissues")
  })
  
  output$boxplot_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(tpm_long_data))
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f9ec"),
                 tags$p("Select a gene and press Explore")))
    plotOutput("boxplot", height = "350px")
  })
  
  output$boxplot <- renderPlot({
    req(selected_gene(), tpm_long_data)
    gene_data <- tpm_long_data %>%
      dplyr::filter(gene == selected_gene(), !is.na(tissue), !is.na(condition))
    validate(need(nrow(gene_data) > 0,
                  paste0("'", selected_gene(), "' was not found in the expression data.")))
    ggplot(gene_data,
           aes(x = tissue, y = log2(tpm + 1), fill = condition, color = condition)) +
      geom_boxplot(alpha = 0.25, outlier.shape = NA, linewidth = 0.6,
                   width = 0.55, position = position_dodge(0.7), color = "black") +
      geom_jitter(aes(group = condition),
                  position = position_jitterdodge(jitter.width = 0.35, dodge.width = 0.7),
                  size = 1.5, alpha = 0.75) +
      scale_fill_manual(values  = c("ALS" = PAL_ALS,  "Control" = PAL_CTRL)) +
      scale_color_manual(values = c("ALS" = PAL_ALS,  "Control" = PAL_CTRL)) +
      labs(x = NULL, y = "log\u2082(TPM + 1)", fill = NULL, color = NULL) +
      als_theme() +
      theme(axis.text.x = element_text(hjust = 0.5, size = 15))
  }, bg = PAL_PANEL)
  
  output$deg_table_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(deg_all_data))
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f4cb"),
                 tags$p("Select a gene and press Explore to see DEG results")))
    DTOutput("deg_table")
  })
  
  output$deg_table <- renderDT({
    req(selected_gene(), deg_all_data)
    g <- selected_gene()
    tbl <- deg_all_data %>%
      dplyr::filter(symbol == g) %>%
      dplyr::select(
        Tissue   = tissue,
        `log2FC` = dplyr::any_of(c("log2FoldChange", "logFC", "LFC")),
        `p-value`= dplyr::any_of(c("pvalue", "PValue", "pval", "P.Value")),
        `FDR`    = dplyr::any_of(c("padj", "FDR", "adj.P.Val", "p.adjust"))
      ) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))) %>%
      tidyr::pivot_longer(cols = -Tissue, names_to = "Metric", values_to = "Value") %>%
      tidyr::pivot_wider(names_from = Tissue, values_from = Value) %>%
      dplyr::rename(` ` = Metric)
    validate(need(nrow(tbl) > 0, paste0("'", g, "' not found in DEG results.")))
    
    # datatable(tbl, rownames = FALSE,
    #           options = list(dom = "t", ordering = FALSE,
    #                          columnDefs = list(list(className = "dt-center",
    #                                                 targets = seq_len(ncol(tbl) - 1)))),
    #           class = "display compact") %>%
    #   formatStyle(" ", target = "row",
    #               backgroundColor = styleEqual("FDR", "rgba(26,122,74,0.06)"))
    
    ## ----------- Test new formating option to color LFC to match theme
    
    tissue_cols <- setdiff(colnames(tbl), " ")
    
    tbl_display <- tbl %>%
      mutate(across(all_of(tissue_cols), ~ {
        val   <- suppressWarnings(as.numeric(.x))
        label <- tbl[[" "]][row_number()]
        dplyr::case_when(
          label == "log2FC" & !is.na(val) & val < -1 ~
            paste0('<span style="color:', PAL_CTRL, '">', sprintf("%.2f", val), '</span>'),
          label == "log2FC" & !is.na(val) & val >  1 ~
            paste0('<span style="color:', PAL_ALS,  '">', sprintf("%.2f", val), '</span>'),
          label == "log2FC" & !is.na(val) ~
            paste0('<span style="color:', PAL_MUTED, '">', sprintf("%.2f", val), '</span>'),
          (label == "FDR" | label == "p-value") & !is.na(val) & val < 0.05 ~
            paste0('<span style="font-weight:600">', sprintf("%.2f", val), '</span>'),
          TRUE ~ .x
        )
      }))
    
    datatable(tbl_display, rownames = FALSE, escape = FALSE,
              options = list(dom = "t", ordering = FALSE,
                             columnDefs = list(list(className = "dt-center",
                                                    targets = seq_len(ncol(tbl_display) - 1)))),
              class = "display compact")
    
  ## ---------------------------------------------------------------------------
    
  })
  
  # ════════════════════════════════════════════════════════
  # TAB 2 — ALS Only Analysis
  # ════════════════════════════════════════════════════════
  
  co_selected <- eventReactive(list(input$co_go, input$co_tissue), {
    req(input$co_tissue)
    list(gene = trimws(toupper(input$co_gene_input)), tissue = input$co_tissue)
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # ALS data filtered to selected gene + tissue
  co_gene_tissue_data <- reactive({
    sel <- co_selected()
    req(sel$gene != "")
    dplyr::filter(tpm_als_only, gene == sel$gene, tissue == sel$tissue)
  })
  
  # ── panel titles ─────────────────────────────────────────
  output$co_boxplot_title <- renderUI({
    sel <- co_selected()
    if (is.null(sel) || is.null(sel$gene) || sel$gene == "") return("Select a gene and tissue")
    paste0(sel$gene, "  \u2014  C9orf72 Status  \u00b7  ", sel$tissue)
  })
  output$co_age_title <- renderUI({
    sel <- co_selected()
    if (is.null(sel) || is.null(sel$gene) || sel$gene == "") return("Age at Death")
    paste0(sel$gene, "  \u2014  Age at Death  \u00b7  ", sel$tissue)
  })
  output$co_dur_title <- renderUI({
    sel <- co_selected()
    if (is.null(sel) || is.null(sel$gene) || sel$gene == "") return("Disease Duration")
    paste0(sel$gene, "  \u2014  Disease Duration  \u00b7  ", sel$tissue)
  })
  
  # ── C9orf72 boxplot ──────────────────────────────────────
  output$co_c9_boxplot_ui <- renderUI({
    sel <- co_selected()
    if (is.null(sel) || is.null(sel$gene) || sel$gene == "")
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f9ec"),
                 tags$p("Select a gene and tissue to begin")))
    plotOutput("co_c9_boxplot", height = "350px")
  })
  
  output$co_c9_boxplot <- renderPlot({
    d <- co_gene_tissue_data()
    validate(
      need(nrow(d) > 0,                    "No data for this gene/tissue combination."),
      need("c9orf72" %in% colnames(d),     "C9orf72 status not available in metadata.")
    )
    d <- dplyr::filter(d, !is.na(c9orf72))
    d$c9orf72 <- factor(d$c9orf72)
    ggplot(d, aes(x = c9orf72, y = log2(tpm + 1), fill = c9orf72, color = c9orf72)) +
      geom_boxplot(alpha = 0.25, outlier.shape = NA, linewidth = 0.6, width = 0.45, color = "black") +
      geom_jitter(width = 0.25, size = 1.5, alpha = 0.75) +
      scale_fill_manual(values  = c("C9 +" = PAL_C9_YES, "C9 -" = PAL_C9_NO)) +
      scale_color_manual(values = c("C9 +" = PAL_C9_YES, "C9 -" = PAL_C9_NO)) +
      labs(x = "C9orf72 Status", y = "log\u2082(TPM + 1)", fill = NULL, color = NULL) +
      als_theme()
  }, bg = PAL_PANEL)
  
  # ── Age at death scatter ─────────────────────────────────
  output$co_age_scatter_ui <- renderUI({
    sel <- co_selected()
    if (is.null(sel) || is.null(sel$gene) || sel$gene == "")
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f4c5"),
                 tags$p("Select a gene and tissue")))
    plotOutput("co_age_scatter", height = "350px")
  })
  
  output$co_age_scatter <- renderPlot({
    d <- co_gene_tissue_data()
    validate(
      need(nrow(d) > 0,                        "No data for this gene/tissue combination."),
      need("age_at_death" %in% colnames(d),    "Age at death not available in metadata.")
    )
    d <- dplyr::filter(d, !is.na(age_at_death))
    ggplot(d, aes(x = age_at_death, y = log2(tpm + 1))) +
      geom_point(color = "#C0392B", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "lm", color = PAL_ACCENT, fill = PAL_ACCENT,
                  alpha = 0.15, linewidth = 0.9) +
      labs(x = "Age at Death (years)", y = "log\u2082(TPM + 1)") +
      als_theme()
  }, bg = PAL_PANEL)
  
  # ── Disease duration scatter ─────────────────────────────
  output$co_dur_scatter_ui <- renderUI({
    sel <- co_selected()
    if (is.null(sel) || is.null(sel$gene) || sel$gene == "")
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f4c5"),
                 tags$p("Select a gene and tissue")))
    plotOutput("co_dur_scatter", height = "350px")
  })
  
  output$co_dur_scatter <- renderPlot({
    d <- co_gene_tissue_data()
    validate(
      need(nrow(d) > 0,                          "No data for this gene/tissue combination."),
      need("disease_duration" %in% colnames(d),  "Disease duration not available in metadata.")
    )
    d <- dplyr::filter(d, !is.na(disease_duration))
    ggplot(d, aes(x = disease_duration, y = log2(tpm + 1))) +
      geom_point(color = "#C0392B", alpha = 0.7, size = 1.5) +
      geom_smooth(method = "lm", color = PAL_MUTED, fill = PAL_MUTED,
                  alpha = 0.15, linewidth = 0.9) +
      labs(x = "Disease Duration (months)", y = "log\u2082(TPM + 1)") +
      als_theme()
  }, bg = PAL_PANEL)
  
  # ── ALS-only DEG table ──────────────────────────────────
  output$co_deg_table_ui <- renderUI({
    sel <- co_selected()
    if (is.null(sel$gene) || sel$gene == "" || length(case_only_data) == 0)
      return(div(class = "empty-state",
                 div(class = "es-icon", "\U0001f4cb"),
                 tags$p("Select a gene and tissue to see DEG results")))
    DTOutput("co_deg_table")
  })
  
  output$co_deg_table <- renderDT({
    sel <- co_selected()
    req(sel$gene != "", length(case_only_data) > 0)
    g <- sel$gene
    t <- sel$tissue
    
    pull_tbl <- function(nm, label) {
      df <- case_only_data[[nm]]
      if (is.null(df)) return(NULL)
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
      pull_tbl("age_at_death",     "Age at Death"),
      pull_tbl("disease_duration", "Disease Duration")
    )
    
    validate(need(nrow(tbl) > 0,
                  paste0("'", g, "' not found in ALS-only DEG results for ", t, ".")))
    
    datatable(tbl, rownames = FALSE,
              options = list(dom = "t", ordering = FALSE,
                             columnDefs = list(list(className = "dt-center", targets = 1:3))),
              class = "display compact") %>%
      formatStyle("FDR",
                  color      = styleInterval(c(0.05), c("#1A7A4A", PAL_MUTED)),
                  fontWeight = styleInterval(c(0.05), c("600", "normal"))) %>%
      formatStyle("log2FC",
                  color = styleInterval(c(-1, 1), c(PAL_CTRL, PAL_MUTED, PAL_ALS)))
  })
}

shinyApp(ui = ui, server = server)
