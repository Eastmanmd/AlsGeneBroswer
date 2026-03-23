# ============================================================
#  ALS Differential Gene Expression Explorer — Shiny App
#  Layout: dark clinical aesthetic with accent colors
# ============================================================

library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(DT)

# ── colour palette ──────────────────────────────────────────
PAL_ALS     <- "#E05C5C"
PAL_CTRL    <- "#5C9BE0"
PAL_BG      <- "#0F1117" 
PAL_PANEL   <-  "#1A1D27"
PAL_BORDER  <- "#2A2D3A"
PAL_TEXT    <- "#E8EAF0"
PAL_MUTED   <- "#7A7E99"
PAL_ACCENT  <- "#7EB8F7"

# ── helpers ─────────────────────────────────────────────────

load_deg_results <- function(deg_dir = "data/case_vs_control") {
  files <- list.files(deg_dir, pattern = "_als_vs_control\\.csv$", full.names = TRUE)
  if (length(files) == 0) { warning("No DEG files found in: ", deg_dir); return(NULL) }
  
  purrr::map_dfr(files, function(f) {
    tissue <- sub("_als_vs_control\\.csv$", "", basename(f))
    df <- suppressMessages(readr::read_csv(f, show_col_types = FALSE))
    if (!"gene" %in% colnames(df)) df <- dplyr::rename(df, gene = 1)
    df$tissue <- tissue
    df
  })
}

tidy_tpm <- function(tpm_path = "data/tpm.csv", meta_path = "data/metadata.csv") {
  tpm  <- suppressMessages(readr::read_csv(tpm_path,  show_col_types = FALSE)) %>%
    column_to_rownames(var = "...1")
  meta <- suppressMessages(readr::read_csv(meta_path, show_col_types = FALSE)) %>%
    column_to_rownames(var = "...1")
  meta$sample <- rownames(meta)
  
  tpm_long <- tpm %>%
    rownames_to_column(var = "gene") %>%
    tidyr::pivot_longer(cols = -gene, names_to = "sample", values_to = "tpm")
  
  meta_slim <- meta[, c("sample", "Sample.Source", "Subject.Group")]
  colnames(meta_slim) <- c("sample", "tissue", "condition")
  
  dplyr::left_join(tpm_long, meta_slim, by = "sample")
}

# ── UI ──────────────────────────────────────────────────────
ui <- fluidPage(
  
  title = "ALS Gene Expression Explorer",
  
  tags$head(
    tags$link(
      href = "https://fonts.googleapis.com/css2?family=Space+Mono:wght@400;700&family=DM+Sans:wght@300;400;600&display=swap",
      rel  = "stylesheet"
    ),
    tags$style(HTML(sprintf("
      * { box-sizing: border-box; margin: 0; padding: 0; }
      body {
        background: %s; color: %s;
        font-family: 'DM Sans', sans-serif; font-size: 14px; min-height: 100vh;
      }
      .app-header {
        background: %s; border-bottom: 1px solid %s;
        padding: 18px 32px; display: flex; align-items: center; gap: 16px;
      }
      .app-header .logo-mark {
        width: 36px; height: 36px;
        background: linear-gradient(135deg, %s 0%%, %s 100%%);
        border-radius: 8px; display: flex; align-items: center; justify-content: center;
        font-family: 'Space Mono', monospace; font-size: 13px; font-weight: 700;
        color: #fff; flex-shrink: 0;
      }
      .app-header h1 {
        font-family: 'Space Mono', monospace; font-size: 16px; font-weight: 700;
        letter-spacing: 0.04em; color: %s;
      }
      .app-header .subtitle { font-size: 12px; color: %s; margin-top: 2px; }
      .main-layout {
        display: grid; grid-template-columns: 280px 1fr;
        gap: 0; height: calc(100vh - 73px);
      }
      .sidebar {
        background: %s; border-right: 1px solid %s;
        padding: 24px 20px; overflow-y: auto;
      }
      .sidebar-section { margin-bottom: 28px; }
      .sidebar-label {
        font-family: 'Space Mono', monospace; font-size: 10px; font-weight: 700;
        letter-spacing: 0.12em; text-transform: uppercase; color: %s;
        margin-bottom: 10px; display: block;
      }

      /* ── gene text input ── */
      #gene_input {
        width: 100%%;
        background: %s !important; border: 1px solid %s !important;
        color: %s !important; border-radius: 6px !important;
        font-family: 'Space Mono', monospace !important;
        font-size: 13px !important; padding: 9px 12px !important;
        outline: none; transition: border-color .15s, box-shadow .15s;
        letter-spacing: 0.06em;
      }
      #gene_input:focus {
        border-color: %s !important;
        box-shadow: 0 0 0 2px rgba(126,184,247,0.15) !important;
      }
      .gene-hint { font-size: 11px; color: %s; margin-top: 6px; line-height: 1.4; }
      .gene-error { font-size: 11px; color: %s; margin-top: 4px; display: none; }

      .btn-search {
        width: 100%%;
        background: linear-gradient(135deg, %s 0%%, %s 100%%);
        border: none; border-radius: 6px; color: #fff;
        font-family: 'Space Mono', monospace; font-size: 12px; font-weight: 700;
        letter-spacing: 0.06em; padding: 10px 16px; cursor: pointer;
        transition: opacity .2s; margin-top: 8px;
      }
      .btn-search:hover { opacity: 0.85; }
      .legend-item {
        display: flex; align-items: center; gap: 8px;
        margin-bottom: 7px; font-size: 12px; color: %s;
      }
      .legend-dot { width: 10px; height: 10px; border-radius: 50%%; flex-shrink: 0; }
      .stat-chip {
        background: %s; border: 1px solid %s; border-radius: 20px;
        padding: 4px 12px; font-size: 11px; color: %s;
        display: inline-block; margin: 3px 3px 3px 0;
      }
      .main-content { display: flex; flex-direction: column; overflow: hidden; background: %s; }
      .tabs-bar {
        display: flex; border-bottom: 1px solid %s;
        padding: 0 24px; background: %s; flex-shrink: 0;
      }
      .tab-btn {
        background: none; border: none; padding: 14px 18px;
        font-family: 'Space Mono', monospace; font-size: 11px; font-weight: 700;
        letter-spacing: 0.08em; text-transform: uppercase; color: %s;
        cursor: pointer; border-bottom: 2px solid transparent;
        transition: color .15s, border-color .15s; position: relative; top: 1px;
      }
      .tab-btn.active { color: %s; border-bottom-color: %s; }
      .tab-btn:hover:not(.active) { color: %s; }
      .tab-content { flex: 1; overflow: auto; padding: 24px; }
      .plot-panel {
        background: %s; border: 1px solid %s;
        border-radius: 10px; padding: 20px; height: 100%%;
      }
      .panel-title {
        font-family: 'Space Mono', monospace; font-size: 11px; font-weight: 700;
        letter-spacing: 0.1em; text-transform: uppercase; color: %s; margin-bottom: 16px;
      }
      .dataTables_wrapper { color: %s !important; font-size: 13px; }
      table.dataTable thead th {
        background: #1e2130 !important; color: %s !important;
        border-bottom: 1px solid %s !important;
        font-family: 'Space Mono', monospace !important;
        font-size: 10px !important; letter-spacing: 0.08em !important;
        text-transform: uppercase !important;
      }
      table.dataTable tbody tr { background: %s !important; }
      table.dataTable tbody tr:nth-child(even) { background: #1c1f2e !important; }
      table.dataTable tbody tr:hover { background: rgba(126,184,247,0.06) !important; }
      table.dataTable tbody td { border-top: 1px solid %s !important; color: %s !important; }
      .dataTables_filter input, .dataTables_length select {
        background: %s !important; border: 1px solid %s !important;
        color: %s !important; border-radius: 4px; padding: 4px 8px;
      }
      .dataTables_info, .dataTables_paginate { color: %s !important; }
      .paginate_button { color: %s !important; }
      .paginate_button.current { background: rgba(126,184,247,0.15) !important; border-radius: 4px; }
      .empty-state {
        display: flex; flex-direction: column;
        align-items: center; justify-content: center;
        height: 300px; gap: 12px; color: %s;
      }
      .empty-state .es-icon { font-size: 40px; opacity: 0.3; }
      .empty-state p { font-size: 13px; }
      ::-webkit-scrollbar { width: 5px; }
      ::-webkit-scrollbar-track { background: transparent; }
      ::-webkit-scrollbar-thumb { background: %s; border-radius: 10px; }
    ",
                            PAL_BG, PAL_TEXT,
                            PAL_PANEL, PAL_BORDER,
                            PAL_ALS, "#9B5CE0",
                            PAL_TEXT, PAL_MUTED,
                            PAL_PANEL, PAL_BORDER,
                            PAL_MUTED,
                            PAL_BG, PAL_BORDER, PAL_TEXT,    # gene input
                            PAL_ACCENT,                        # gene input focus
                            PAL_MUTED,                         # hint
                            PAL_ALS,                           # error
                            PAL_ALS, "#C1555F",               # btn gradient
                            PAL_TEXT,
                            PAL_BG, PAL_BORDER, PAL_MUTED,
                            PAL_TEXT,
                            PAL_BG, PAL_BORDER, PAL_PANEL,
                            PAL_MUTED, PAL_ACCENT, PAL_ACCENT, PAL_TEXT,
                            PAL_PANEL, PAL_BORDER,
                            PAL_MUTED,
                            PAL_TEXT,
                            PAL_ACCENT, PAL_BORDER,
                            PAL_TEXT, PAL_BORDER, PAL_TEXT,
                            PAL_BG, PAL_BORDER, PAL_TEXT,
                            PAL_MUTED, PAL_TEXT,
                            PAL_MUTED,
                            PAL_BORDER
    ))
    ),
    
    # ── Header ────────────────────────────────────────────────
    div(class = "app-header",
        div(class = "logo-mark", "ALS"),
        div(
          tags$h1("Gene Expression Explorer"),
          div(class = "subtitle", "ALS vs Control · Multi-tissue DEG Analysis")
        )
    ),
    
    # ── Main layout ───────────────────────────────────────────
    div(class = "main-layout",
        
        # ── Sidebar ─────────────────────────────────────────────
        div(class = "sidebar",
            
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Gene Search"),
                
                # Manual text entry — no dropdown
                textInput("gene_input", label = NULL, value = "FGR", width = "100%"),
                div(class = "gene-hint", "Enter an official gene symbol (e.g. SOD1, TARDBP, FUS)"),
                div(class = "gene-error", id = "gene_error_msg", "⚠ Gene not found in dataset"),
                
                tags$button(
                  class   = "btn-search",
                  id      = "go_btn",
                  onclick = "Shiny.setInputValue('go', Math.random())",
                  "▶  EXPLORE GENE"
                )
            ),
            
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Colour Legend"),
                div(class = "legend-item",
                    div(class = "legend-dot", style = paste0("background:", PAL_ALS)), "ALS"),
                div(class = "legend-item",
                    div(class = "legend-dot", style = paste0("background:", PAL_CTRL)), "Control")
            ),
            
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Significance Thresholds"),
                div(class = "stat-chip", "padj < 0.05"),
                div(class = "stat-chip", "|log₂FC| > 1")
            ),
            
            div(class = "sidebar-section",
                tags$span(class = "sidebar-label", "Status"),
                uiOutput("data_status")
            )
        ),
        
        # ── Main content ─────────────────────────────────────────
        div(class = "main-content",
            
            div(class = "tabs-bar",
                tags$button(class = "tab-btn active", id = "tab_plot",
                            onclick = "switchTab('plot')", "Expression Plot"),
                tags$button(class = "tab-btn", id = "tab_table",
                            onclick = "switchTab('table')", "DEG Results Table")
            ),
            
            div(class = "tab-content", id = "content_plot",
                div(class = "plot-panel",
                    div(class = "panel-title", uiOutput("plot_title")),
                    uiOutput("boxplot_ui")
                )
            ),
            
            div(class = "tab-content", id = "content_table", style = "display:none",
                div(class = "plot-panel",
                    div(class = "panel-title", "Differential Expression Results by Tissue"),
                    uiOutput("deg_table_ui")
                )
            )
        )
    ),
    
    tags$script(HTML("
    // Tab switching
    function switchTab(name) {
      ['plot','table'].forEach(function(t) {
        document.getElementById('content_' + t).style.display = (t === name) ? 'block' : 'none';
        document.getElementById('tab_' + t).className = 'tab-btn' + (t === name ? ' active' : '');
      });
    }
    // Enter key triggers search
    document.addEventListener('keydown', function(e) {
      if (e.key === 'Enter' && document.activeElement.id === 'gene_input') {
        Shiny.setInputValue('go', Math.random());
      }
    });
    // Handler to show/hide the 'gene not found' message
    Shiny.addCustomMessageHandler('toggleGeneError', function(msg) {
      var el = document.getElementById('gene_error_msg');
      if (el) el.style.display = msg.show ? 'block' : 'none';
    });
  "))
  )
)
  
# ── Server ────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  # ── Load data ONCE at startup ──────────────────────────────
  tpm_long_data <- tryCatch(
    suppressMessages(tidy_tpm("data/tpm.csv", "data/metadata.csv")),
    error = function(e) { message("TPM load error: ", e$message); NULL }
  )
  deg_all_data <- tryCatch(
    suppressMessages(load_deg_results("data/case_vs_control")),
    error = function(e) { message("DEG load error: ", e$message); NULL }
  )
  
  tpm_long <- reactive({ tpm_long_data })
  deg_all  <- reactive({ deg_all_data  })
  
  # ── Selected gene — fires on button click OR Enter key ─────
  selected_gene <- eventReactive(list(input$go, input$gene_input), {
    trimws(toupper(input$gene_input))   # uppercase + strip whitespace
  }, ignoreNULL = FALSE)
  
  # ── Validate gene and show inline error if not found ───────
  observe({
    g   <- selected_gene()
    tpm <- tpm_long()
    if (is.null(tpm) || is.null(g) || g == "") return()
    session$sendCustomMessage("toggleGeneError", list(show = !(g %in% tpm$gene)))
  })
  
  # ── Status chips ──────────────────────────────────────────
  output$data_status <- renderUI({
    tpm <- tpm_long(); deg <- deg_all()
    items <- list()
    if (!is.null(tpm)) {
      items <- c(items, list(
        div(class = "stat-chip", paste0(length(unique(tpm$gene)),   " genes")),
        div(class = "stat-chip", paste0(length(unique(tpm$sample)), " samples")),
        div(class = "stat-chip", paste0(length(unique(tpm$tissue)), " tissues"))
      ))
    }
    if (!is.null(deg))
      items <- c(items, list(
        div(class = "stat-chip", paste0(length(unique(deg$tissue)), " DEG tissues"))
      ))
    if (length(items) == 0)
      div(style = "color:#E05C5C; font-size:12px;", "⚠ Data files not found")
    else
      tagList(items)
  })
  
  # ── Plot title ────────────────────────────────────────────
  output$plot_title <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "") return("Enter a gene symbol to begin")
    paste0("Expression of  ", g, "  across tissues")
  })
  
  # ── Boxplot UI ────────────────────────────────────────────
  output$boxplot_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(tpm_long()))
      return(div(class = "empty-state",
                 div(class = "es-icon", "🧬"),
                 tags$p("Enter a gene symbol and press Explore or hit Enter")))
    plotOutput("boxplot", height = "460px")
  })
  
  output$boxplot <- renderPlot({
    req(selected_gene(), tpm_long())
    gene_data <- tpm_long() %>%
      dplyr::filter(gene == selected_gene(), !is.na(tissue), !is.na(condition))
    
    validate(need(nrow(gene_data) > 0,
                  paste0("'", selected_gene(), "' was not found in the TPM data.\nCheck spelling and try again.")))
    
    ggplot(gene_data,
           aes(x = tissue, y = log2(tpm + 1), fill = condition, color = condition)) +
      geom_boxplot(alpha = 0.35, outlier.shape = NA, linewidth = 0.6,
                   width = 0.55, position = position_dodge(0.7)) +
      geom_jitter(aes(group = condition),
                  position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.7),
                  size = 1.8, alpha = 0.7) +
      scale_fill_manual(values  = c("ALS" = PAL_ALS, "Control" = PAL_CTRL)) +
      scale_color_manual(values = c("ALS" = PAL_ALS, "Control" = PAL_CTRL)) +
      labs(x = NULL, y = "log\u2082(TPM + 1)", fill = NULL, color = NULL) +
      theme_minimal(base_family = "sans") +
      theme(
        plot.background   = element_rect(fill = PAL_PANEL, color = NA),
        panel.background  = element_rect(fill = PAL_PANEL, color = NA),
        panel.grid.major  = element_line(color = PAL_BORDER, linewidth = 0.4),
        panel.grid.minor  = element_blank(),
        axis.text         = element_text(color = PAL_TEXT, size = 11),
        axis.text.x       = element_text(angle = 35, hjust = 1, size = 11),
        axis.title.y      = element_text(color = PAL_MUTED, size = 12, margin = margin(r = 10)),
        legend.position   = "top",
        legend.text       = element_text(color = PAL_TEXT, size = 12),
        legend.background = element_rect(fill = PAL_PANEL, color = NA),
        legend.key        = element_rect(fill = "transparent"),
        plot.margin       = margin(12, 12, 12, 12)
      )
  }, bg = PAL_PANEL)
  
  # ── DEG table UI ──────────────────────────────────────────
  output$deg_table_ui <- renderUI({
    g <- selected_gene()
    if (is.null(g) || g == "" || is.null(deg_all()))
      return(div(class = "empty-state",
                 div(class = "es-icon", "📋"),
                 tags$p("Enter a gene symbol and press Explore to see DEG results")))
    DTOutput("deg_table")
  })
  
  output$deg_table <- renderDT({
    req(selected_gene(), deg_all())
    g <- selected_gene()
    
    tbl <- deg_all() %>%
      dplyr::filter(gene == g) %>%
      dplyr::select(
        Tissue        = tissue,
        `log2FC` = dplyr::any_of(c("log2FoldChange", "logFC", "LFC")),
        `p-value`     = dplyr::any_of(c("pvalue", "PValue", "pval", "P.Value")),
        `adj p-value` = dplyr::any_of(c("padj", "FDR", "adj.P.Val", "p.adjust"))
      ) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4)))
    
    validate(need(nrow(tbl) > 0,
                  paste0("'", g, "' not found in DEG results.")))
    
    datatable(
      tbl, rownames = FALSE,
      options = list(
        pageLength = 15, dom = "ftip",
        order      = list(list(3, "asc")),
        columnDefs = list(list(className = "dt-center", targets = 1:3))
      ),
      class = "display compact"
    ) %>%
      formatStyle("adj p-value",
                  color      = styleInterval(c(0.05), c("#6EE7B7", PAL_MUTED)),
                  fontWeight = styleInterval(c(0.05), c("600", "normal"))
      ) %>%
      formatStyle("log\u2082FC",
                  color = styleInterval(c(-1, 1), c(PAL_CTRL, PAL_MUTED, PAL_ALS))
      )
  })
}

shinyApp(ui = ui, server = server)