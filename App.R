#############################################################
#### ----------------------------------------------------------
#### 🧠 Enhancer Network Explorer 
#### ----------------------------------------------------------
#### Version 2.9 (Stable)
#### Author: Juan Jair Santillan-Cigales
#### Date: 23.10.2025
####
#### Description:
#### This Shiny app enables interactive exploration of 
#### enhancer–TF–gene regulatory networks derived from 
#### integrative CUT&RUN and expression analyses.
#############################################################


# All dependencies are loaded at the start to ensure reproducibility
library(shiny)          # Core Shiny framework
library(shinyjqui)      # Resizable and draggable UI elements
library(shinyWidgets)   # Enhanced UI widgets (switchInput, etc.)
library(shinyalert)     # Custom alert messages within Shiny
library(visNetwork)     # Interactive network visualization
library(dplyr)          # Data manipulation and filtering
library(scales)         # Rescaling and color scaling utilities
library(tibble)         # Modern data frame handling
library(igraph)         # Graph layouts and algorithms
library(htmlwidgets)    # Widget saving and export utilities
library(DT)             # Interactive tables
library(bslib)          # Bootstrap theming support


# ============================================================
# Load datasets
# ============================================================
# The app can load default datasets stored in the working directory.
# Additional datasets can be loaded manually by the user before launching
# the app (they will be detected automatically by the next function).
df_filtered_d8 <- readRDS("df_filtered_d8_MERGED_ENCODE_JASPAR.rds")
df_filtered_d16 <- readRDS("df_filtered_d16_MERGED_ENCODE_JASPAR.rds")

# ===============================================
# Dynamic detection of datasets in the global environment
# ===============================================
# This function scans the current R environment for valid dataframes 
# containing the essential columns: enhancer_id, tf_name, connected_gene.
# It returns a vector with the names of all valid datasets available.
get_available_datasets <- function() {
  objs <- ls(envir = .GlobalEnv)
  valid_dfs <- objs[sapply(objs, function(x) {
    obj <- get(x, envir = .GlobalEnv)
    is.data.frame(obj) &&
      all(c("enhancer_id", "tf_name", "connected_gene") %in% colnames(obj))
  })]
  return(valid_dfs)
}

# Example output (when run interactively):
# > get_available_datasets()
# [1] "df_filtered_d8" "df_filtered_d16"

# ===============================================
# Definition of names for dataset selection
# ===============================================
# These labels are shown in the Shiny UI selector instead of 
# raw object names. Additional datasets can be added here manually
# just specify the name of each object preloaded.
dataset_labels <- c(
  "df_filtered_d8"  = "MBO day 8 (BHB vs Ctrl)",
  "df_filtered_d16" = "MBO day 16 (BHB vs Ctrl)"
)

# ============================================================
# Universal Enhancer–TF–Gene Network Generator
# ============================================================
# Function: generate_enhancer_network()
# ------------------------------------------------------------
# Purpose:
#   Builds an interactive visNetwork visualization integrating
#   enhancer–TF–gene relationships from a filtered dataset.
#
# Arguments:
#   df                Data frame containing enhancer–TF–gene relationships.
#   enhancers_of_interest  Character vector of enhancer IDs to include.
#   layout_method     Layout algorithm: "kk" (static) or "physics" (dynamic).
#   cluster_colors    Named color palette for enhancer clusters.
#   fc_limit          Absolute limit for log2FC color scaling.
#   node_font_size    Font size for node labels.
#   node_size         Default node size (if not scaled automatically).
#   node_threshold    Placeholder (legacy parameter; currently not in use).
#   layout_scale      Scaling factor for static layout (Kamada–Kawai).
#   filter_edges      List of numeric thresholds for edge filtering:
#                     list(tf_enh = X, enh_gene = Y, tf_gene = Z)
#
# Returns:
#   An interactive visNetwork object representing the regulatory network.
# ------------------------------------------------------------

generate_enhancer_network <- function(
    df,
    enhancers_of_interest,
    layout_method = c("kk", "physics"),
    cluster_colors = c("1" = "#e6b400", "2" = "#6b6b6b", "3" = "#7aa874"),
    fc_limit = 1.5,
    node_font_size = 20,
    node_size = 25,
    node_threshold = 40,
    layout_scale = NULL,
    filter_edges = list(tf_enh = NULL, enh_gene = NULL, tf_gene = NULL)
) {
  # Initial setup
  layout_method <- match.arg(layout_method)
  
  # Subset only the enhancers of interest
  df_sub <- df %>% filter(enhancer_id %in% enhancers_of_interest)
  
  # Normalize key columns and handle missing values
  # Detect and standardize the log2FC column (for genes)
  log2fc_col <- grep("^log2FoldChange$|^log2FC$", names(df_sub), value = TRUE)
  
  if (length(log2fc_col) == 0) {
    df_sub$log2FC <- NA_real_
  } else {
    chosen_col <- log2fc_col[1]
    df_sub$log2FC <- df_sub[[chosen_col]]
  }
  
  # Ensure TF fold change column exists
  if (!"tf_l2fc" %in% names(df_sub)) df_sub$tf_l2fc <- NA_real_
  
  # Define color palette for log2FC
  fc_palette <- scales::col_numeric(
    palette = c("#377EB8", "white", "#E41A1C"),
    domain = c(-fc_limit, fc_limit)
  )
  
  # Automatically scale node size depending on network size
  if (is.null(node_size)) {
    n_total <- length(unique(c(df_sub$enhancer_id, df_sub$tf_name, df_sub$connected_gene)))
    node_size <- case_when(
      n_total < 100  ~ 40,
      n_total < 300  ~ 30,
      n_total < 800  ~ 20,
      TRUE            ~ 10
    )
  }
  
  # ====================== NODES ======================
  # Nodes construction
  nodes <- tibble(
    id = unique(c(df_sub$enhancer_id, df_sub$tf_name, df_sub$connected_gene))
  ) %>%
    mutate(
      # Define node type
      type = case_when(
        id %in% df_sub$enhancer_id ~ "Enhancer",
        id %in% df_sub$tf_name     ~ "TF",
        TRUE                       ~ "Gene"
      ),
      # Optional motif source annotation for TFs
      motif_source_node = if ("motif_source" %in% names(df_sub)) {
        df_sub$motif_source[match(id, df_sub$tf_name)]
      } else {
        NA_character_
      },
      label = id,
      # Map fold-change and cluster information
      log2FC = ifelse(type == "Gene", df_sub$log2FC[match(id, df_sub$connected_gene)], NA_real_),
      tf_l2fc_node = ifelse(type == "TF", df_sub$tf_l2fc[match(id, df_sub$tf_name)], NA_real_),
      cluster_node = ifelse(type == "Enhancer", df_sub$cluster[match(id, df_sub$enhancer_id)], NA_real_),
      # Truncate FC values to avoid over-saturation of colors
      log2FC_trunc  = pmax(pmin(log2FC, fc_limit), -fc_limit),
      tf_l2fc_trunc = pmax(pmin(tf_l2fc_node, fc_limit), -fc_limit),
      # Assign node colors
      color = case_when(
        type == "Enhancer" ~ cluster_colors[as.character(cluster_node)],
        type == "Gene" & !is.na(log2FC)      ~ fc_palette(log2FC_trunc),
        type == "TF"   & !is.na(tf_l2fc_node)~ fc_palette(tf_l2fc_trunc),
        TRUE ~ "#bdbdbd"
      ),
      
      # Tooltip (hover text)
      title = case_when(
        type == "Gene" ~ paste0("<b>", label, "</b><br>Type: Gene<br>log2FC: ", round(log2FC, 3)),
        type == "TF"   ~ paste0("<b>", label, "</b><br>Type: TF",
                                ifelse(!is.na(motif_source_node),
                                       paste0("<br>Motif source: ", motif_source_node), ""),
                                "<br>log2FC: ", round(tf_l2fc_node, 3)),
        type == "Enhancer" ~ paste0("<b>", label, "</b><br>Type: Enhancer<br>Cluster: ", cluster_node)
      ),
      font.size = node_font_size,
      size = node_size
    )
  
  # ====================== EDGES ======================
  # Edges construction
  # Define three edge categories:
  #   TF → Enhancer, Enhancer → Gene, TF → Gene (indirect)
  edges_tf_enh <- df_sub %>%
    select(from = tf_name, to = enhancer_id, weight = motif_score) %>%
    group_by(from, to) %>%
    summarise(weight = max(weight, na.rm = TRUE), .groups = "drop") %>%
    mutate(edge_type = "tf_enh")
  
  edges_enh_gene <- df_sub %>%
    select(from = enhancer_id, to = connected_gene, weight = score) %>%
    group_by(from, to) %>%
    summarise(weight = max(weight, na.rm = TRUE), .groups = "drop") %>%
    mutate(edge_type = "enh_gene")
  
  edges_tf_gene <- df_sub %>%
    select(from = tf_name, to = connected_gene, weight = motif_score) %>%
    group_by(from, to) %>%
    summarise(weight = max(weight, na.rm = TRUE), .groups = "drop") %>%
    mutate(edge_type = "tf_gene")
  
  edges <- bind_rows(edges_tf_enh, edges_enh_gene, edges_tf_gene)
  
  # ====================== FILTERS ======================
  # Apply edge filters
  if (!is.null(filter_edges$tf_enh)) {
    edges <- edges %>% filter(!(edge_type == "tf_enh" & weight < filter_edges$tf_enh))
  }
  if (!is.null(filter_edges$enh_gene)) {
    edges <- edges %>% filter(!(edge_type == "enh_gene" & weight < filter_edges$enh_gene))
  }
  if (!is.null(filter_edges$tf_gene)) {
    edges <- edges %>% filter(!(edge_type == "tf_gene" & weight < filter_edges$tf_gene))
  }
  
  # Add visualization parameters to edges
  edges <- edges %>%
    mutate(
      title = case_when(
        edge_type == "tf_enh"   ~ paste0("TF→Enhancer | motif_score: ", round(weight, 2)),
        edge_type == "enh_gene" ~ paste0("Enhancer→Gene | score: ", round(weight, 2)),
        edge_type == "tf_gene"  ~ paste0("TF→Gene (indirect) | motif_score: ", round(weight, 2))
      ),
      color = case_when(
        edge_type == "tf_enh"   ~ "#8D69B8",
        edge_type == "enh_gene" ~ "#FFA500",
        edge_type == "tf_gene"  ~ "transparent"
      ),
      width = case_when(
        edge_type == "tf_gene" ~ 0.1,
        TRUE ~ scales::rescale(weight, to = c(1, 8), from = range(weight, na.rm = TRUE))
      ),
      id = paste(edge_type, from, to, sep = "|")
    )
  
  # ====================== VISUALIZATION ======================
  # Dynamic force-based layout (interactive)
  if (layout_method == "physics") {
    vis <- visNetwork(nodes, edges) %>%
      visEdges(arrows = "to") %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visInteraction(navigationButtons = TRUE) %>%
      visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(
          gravitationalConstant = -150,
          centralGravity = 0.003,
          springLength = 200,
          springConstant = 0.02,
          avoidOverlap = 1
        ),
        stabilization = list(enabled = TRUE, iterations = 2500)
      ) %>%
      visEvents(stabilized = "function () { this.setOptions({ physics: {enabled: false} }); }")
  
  # Static Kamada–Kawai layout  
  } else if (layout_method == "kk") {
    edges_for_ig <- edges %>% select(from, to) %>% distinct()
    g <- graph_from_data_frame(edges_for_ig, vertices = nodes, directed = TRUE)
    
    lay <- layout_with_kk(g, weights = NA)
    if (is.null(layout_scale)) {
      layout_scale <- ifelse(nrow(nodes) <= 50, 300,
                             ifelse(nrow(nodes) <= 150, 600, 1000))
    }
    
    nodes$x <- lay[, 1] * layout_scale
    nodes$y <- lay[, 2] * layout_scale
    
    vis <- visNetwork(nodes, edges) %>%
      visNodes(physics = FALSE, fixed = TRUE) %>%
      visEdges(arrows = "to", smooth = FALSE) %>%
      visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
      visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE) %>%
      visPhysics(enabled = FALSE)
  }
  
  return(vis)
}


### Enhancer Network Explorer v2.9 ----

# ===========================================================
# Enhancer Network Explorer — UI (multi-step flow)
# ===========================================================
# This user interface defines a clear, guided workflow divided into:
#   (1) Dataset selection
#   (2) Guided search
#   (3) Filtering parameters
#   (4) Visualization
#
# Each step is visually and functionally separated for better usability.
# The design uses a 'flatly' bootswatch theme and custom CSS 
# for styling and consistency.
# =========================================================

ui <- fluidPage(
  # Theme and global setup
  theme = bs_theme(bootswatch = "flatly"),
  titlePanel("🧠 Enhancer Network Explorer"),
  
  tags$head(
    tags$script(src = "https://html2canvas.hertzen.com/dist/html2canvas.min.js")
  ),
  
  # Custom CSS (styling adjustments)
  tags$style(HTML("
  /* --- Graph container style --- */
  #enh_network {
    background-color: white !important;
    border-radius: 8px;
    box-shadow: 0 0 10px rgba(0,0,0,0.1);
  }

  /* --- Switch aesthetic customization --- */
  .bootstrap-switch {
    border: 1.5px solid #D5DBDB !important;
    border-radius: 6px !important;
    box-shadow: inset 0 1px 2px rgba(0,0,0,0.05);
  }

  .bootstrap-switch .bootstrap-switch-handle-on {
    background: #1ABC9C !important;   /* teal tone when ON */
    color: white !important;
    font-weight: 500;
  }

  .bootstrap-switch .bootstrap-switch-handle-off {
    background: #E5E8E8 !important;   /* light gray when OFF */
    color: #34495E !important;
    font-weight: 500;
  }

  .bootstrap-switch .bootstrap-switch-label {
    background: #F8F9F9 !important;
    color: #5D6D7E !important;
  }
  
  /* === Global style tweaks === */
  h4 {
    font-size: 16px !important;
    font-weight: 600 !important;
    margin-top: 15px;
    margin-bottom: 10px;
  }

  /* Slightly smaller for subheadings inside boxes */
  .box h4 {
    font-size: 15px !important;
    font-weight: 500 !important;
  }
")),
  
  # Layout: Sidebar + Main panel
  sidebarLayout(
    
    # ======================================================
    # Sidebar Panel
    # ======================================================
    sidebarPanel(
      width = 3,
      
      # --- Instructions ---
      h4("📘 Instructions"),
      tags$details(
        tags$summary("Click to expand / collapse instructions"),
        p("This interactive app allows you to explore transcriptional regulatory networks 
           between transcription factors (TFs), enhancers, and their connected genes."),
        tags$ul(
          tags$li(strong("1. Select dataset:"), " Choose the dataset to analyze."),
          tags$li(strong("2. Guided search:"), " Select TF, Gene, or Enhancer to start exploration."),
          tags$li(strong("3. Add or combine networks:"), " You can add multiple contexts before generating the network."),
          tags$li(strong("4. Adjust filters and visualization:"), " Refine thresholds or enable physics."),
          tags$li(strong("5. Generate network:"), " Visualize and inspect interactions in detail."),
          tags$li(strong("⚠️ Tip:"), " Large networks (>200 nodes) may take longer to render.")
        )
      ),
      hr(),
      
      # --- Dataset selection (dynamic) ---
      h4("1️⃣ Dataset selection"),
      uiOutput("dataset_selector"),
      hr(),
      
      # --- Guided search (multi-step) ---
      h4("2️⃣ Guided search (multi-step)"),
      uiOutput("primary_type_ui"),
      uiOutput("primary_value_ui"),
      uiOutput("secondary_type_ui"),
      uiOutput("secondary_value_ui"),
      uiOutput("enhancer_selection_ui"),
      
      # --- Add / Clear networks side-by-side ---
      br(),
      fluidRow(
        column(
          6,
          actionButton(
            "add_network",
            "➕ Add network",
            class = "btn-secondary",
            style = "width:100%; font-weight:500; background-color:#5D6D7E; border:none; color:white;"
          )
        ),
        column(
          6,
          actionButton(
            "clear_all",
            "Clear all",
            class = "btn-warning",
            style = "width:100%; font-weight:500; background-color:#F39C12; border:none; color:white;"
          )
        )
      ),
      
      # --- Summary of added network contexts ---
      br(),
      div(
        style = "margin-top: 10px; padding: 10px; border: 1px solid #dcdcdc; border-radius: 6px; background-color: #fcfcfc;",
        h5("🧩 Added network contexts"),
        uiOutput("added_networks_summary")
      ),
      br(),
      
      # ======================================================
      # Filtering parameters
      # ======================================================
      h4("3️⃣ Filtering parameters"),
      tags$details(
        tags$summary(HTML("<b style='color:#34495E;'>Click to expand / collapse filters</b>")),
        br(),
        sliderInput("motif_cutoff", "Minimum motif score (TF → Enhancer):",
                    min = 0, max = 50, value = 0, step = 0.5),
        sliderInput("enhgene_cutoff", "Minimum enhancer score (Enhancer → Gene):",
                    min = 0, max = 700, value = 0, step = 1),
        sliderInput("tf_fc_range", HTML("<b>TF log2FC range:</b>"),
                    min = -5, max = 5, value = c(-5, 5), step = 0.1),
        sliderInput("gene_fc_range", HTML("<b>Target gene log2FC range:</b>"),
                    min = -5, max = 5, value = c(-5, 5), step = 0.1),
        sliderInput("font_size", HTML("<b>Label size:</b>"),
                    min = 10, max = 60, value = 35, step = 2),
        sliderInput("node_size", HTML("<b>Node size:</b>"), 
                    min = 5, max = 80, value = 25, step = 1)
      ),
      hr(),
      
      # ======================================================
      # Visualization controls
      # ======================================================
      
      h4("4️⃣ Visualization"),
      
      fluidRow(
        column(
          7,
          actionButton(
            "generate",
            "🔍 Generate network",
            class = "btn-primary",
            style = "width:100%; height:50px; font-weight:500; background-color:#34495E; border:none; color:white; font-size:15px;"
          )
        ),
        column(
          5,
          switchInput(
            inputId = "use_physics",
            label = HTML("<b>Switch mode:</b>"),
            value = FALSE,
            onLabel = "Physics",
            offLabel = "Static",
            size = "small",
            handleWidth = 80,
            labelWidth = 90
          )
        )
      ),
      hr(),
      
      # --- Footer ---
      tags$footer(
        align = "center",
        style = "font-size: 12px; color: grey;",
        HTML("Developed by <b>Juan Jair Santillán-Cigales</b><br>
             PhD Candidate, Instituto de Fisiología Celular (IFC-UNAM)")
      )
    ),
    
    # ======================================================
    # Main Panel
    # ======================================================
    mainPanel(
      width = 9,
      h4("🌐 Network visualization"),
      textOutput("current_dataset_label", inline = TRUE),
      br(),
      
      # Snapshot button (HTML2Canvas)
      fluidRow(
        column(
          12,
          div(
            style = "display: flex; justify-content: flex-end; gap: 10px; margin-bottom: 8px;",
            actionButton("snapshot_btn", "📸 Snapshot", class = "btn-info")
          )
        )
      ),
      
      # Custom JS to capture the visNetwork snapshot
      tags$script(HTML("
        Shiny.addCustomMessageHandler('captureVisNetwork', function(message) {
          var networkDiv = document.getElementById('enh_network');
          if (networkDiv) {
            html2canvas(networkDiv).then(function(canvas) {
              var link = document.createElement('a');
              link.download = 'enhancer_network_snapshot.png';
              link.href = canvas.toDataURL('image/png');
              link.click();
            });
          }
        });
      ")),
      
      # Network container
      textOutput("network_summary"),
      jqui_resizable(
        div(
          id = "network_container",
          style = "
            height: 70vh;
            overflow: hidden;
            resize: vertical;
            width: 100%;
            margin: 0 auto;
            position: relative;
            border-radius: 8px;
            box-shadow: 0 0 8px rgba(0,0,0,0.15);
          ",
          visNetworkOutput("enh_network", height = "100%"),
          div(
            style = "
              position: absolute;
              bottom: -6px;
              left: 50%;
              transform: translateX(-50%);
              width: 40px;
              height: 6px;
              background-color: #bbb;
              border-radius: 3px;
              cursor: ns-resize;
            "
          )
        ),
        options = list(handles = "s")
      ),
      
      # Interaction table
      br(),
      h4("📊 TF–Enhancer–Gene interaction table"),
      DTOutput("interaction_table")
    )
  )
)

# =========================================================
# Enhancer Network Explorer — SERVER
# =========================================================

server <- function(input, output, session) {
  
  # ===============================================
  # Dynamic dataset selector
  # ===============================================
  # This block renders the dataset selection dropdown dynamically
  # based on datasets available in the R environment.
  # It automatically detects valid data frames (via get_available_datasets()).
  # The user can load additional datasets externally, and they will appear here.
  output$dataset_selector <- renderUI({
    available_datasets <- get_available_datasets()
    
    # # Combine detected datasets with pre-defined names
    choices_named <- setNames(
      available_datasets,
      ifelse(
        available_datasets %in% names(dataset_labels),
        dataset_labels[available_datasets],
        available_datasets  # fallback: use raw object name
      )
    )
    
    # If no valid dataset is found, display a warning message
    if (length(choices_named) == 0) {
      return(tags$div(
        style = "color:#A93226; font-weight:500;",
        "⚠️ No valid datasets found. Please load at least one dataframe containing enhancer_id, tf_name, and connected_gene."
      ))
    }
    
    # Render dropdown menu
    selectInput(
      inputId = "dataset_choice",
      label = "Select dataset:",
      choices = choices_named,
      selected = ifelse(
        "df_filtered_d8" %in% available_datasets, "df_filtered_d8",
        available_datasets[1]
      )
    )
  })
  
  
  # =========================================================
  # Dataset loading and normalization
  # =========================================================
  # This reactive expression retrieves the selected dataset,
  # ensures it contains essential columns, and standardizes
  # common variable names for compatibility with downstream steps.
  df_base <- reactive({
    req(input$dataset_choice)  # ensure user selection
    df <- get(input$dataset_choice, envir = .GlobalEnv)
    
    # --- Check required columns ---
    required_cols <- c("tf_name", "connected_gene", "enhancer_id")
    missing_cols <- setdiff(required_cols, names(df))
    if (length(missing_cols) > 0)
      stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
    
    # --- Normalize log2FC for genes ---
    log2fc_col <- grep("^log2FoldChange$|^log2FC$", names(df), value = TRUE)
    if (length(log2fc_col) == 0) {
      df$log2FC <- NA_real_
    } else {
      df$log2FC <- df[[log2fc_col[1]]]
    }
    
    # --- Normalize TF log2FC ---
    if (!"tf_l2fc" %in% names(df)) df$tf_l2fc <- NA_real_
    
    # --- Ensure motif_source column exists ---
    if (!"motif_source" %in% names(df)) df$motif_source <- NA_character_
    
    df
  })
  
  # Backward compatibility alias (legacy support)
  current_df <- df_base
  
  # =========================================================
  # Dynamic multi-step guided selection
  # =========================================================
  # This section implements a structured exploration workflow:
  # Step 1: choose exploration type (TF / Gene / Enhancer)
  # Step 2: select specific element(s)
  # Step 3: choose secondary type to refine relationships
  # Step 4: select corresponding secondary element(s)
  # Step 5: optionally choose connecting enhancers / TFs / genes
  
  # ---------------------------------------------------------
  # Step 1: Primary exploration choice
  output$primary_type_ui <- renderUI({
    selectInput(
      "primary_type", "Step 1️⃣ What would you like to explore?",
      choices = c("Transcription factor (TF)" = "TF",
                  "Target gene" = "Gene",
                  "Enhancer" = "Enhancer"),
      selected = "TF"
    )
  })
  
  # Step 2: Primary selection input (dependent on type)
  output$primary_value_ui <- renderUI({
    req(input$primary_type)
    df <- current_df()
    
    if (input$primary_type == "TF") {
      selectizeInput("primary_value", "Select transcription factor(s):",
                     choices = sort(unique(df$tf_name)), multiple = TRUE)
    } else if (input$primary_type == "Gene") {
      selectizeInput("primary_value", "Select target gene(s):",
                     choices = sort(unique(df$connected_gene)), multiple = TRUE)
    } else {
      selectizeInput("primary_value", "Select enhancer(s):",
                     choices = sort(unique(df$enhancer_id)), multiple = TRUE)
    }
  })
  
  # Step 3: Secondary type (what to connect to)
  output$secondary_type_ui <- renderUI({
    req(input$primary_value)
    secondary_choices <- c("TF", "Gene", "Enhancer")
    secondary_choices <- setdiff(secondary_choices, input$primary_type)
    
    display_names <- c("TF" = "Transcription factor",
                       "Gene" = "Target gene",
                       "Enhancer" = "Enhancer")
    
    selectInput(
      "secondary_type",
      "Step 2️⃣ How would you like to continue?",
      choices = setNames(secondary_choices, display_names[secondary_choices]),
      selected = secondary_choices[1]
    )
  })
  
  # Step 4: Secondary selection input (depends on combination)
  output$secondary_value_ui <- renderUI({
    req(input$primary_value, input$secondary_type)
    df <- current_df()
    
    if (input$primary_type == "TF" && input$secondary_type == "Gene") {
      sub <- df %>% filter(tf_name %in% input$primary_value)
      options <- sort(unique(sub$connected_gene))
      label <- "Select connected target gene(s):"
    } else if (input$primary_type == "TF" && input$secondary_type == "Enhancer") {
      sub <- df %>% filter(tf_name %in% input$primary_value)
      options <- sort(unique(sub$enhancer_id))
      label <- "Select connected enhancer(s):"
    } else if (input$primary_type == "Gene" && input$secondary_type == "TF") {
      sub <- df %>% filter(connected_gene %in% input$primary_value)
      options <- sort(unique(sub$tf_name))
      label <- "Select regulating transcription factor(s):"
    } else if (input$primary_type == "Gene" && input$secondary_type == "Enhancer") {
      sub <- df %>% filter(connected_gene %in% input$primary_value)
      options <- sort(unique(sub$enhancer_id))
      label <- "Select associated enhancer(s):"
    } else if (input$primary_type == "Enhancer" && input$secondary_type == "TF") {
      sub <- df %>% filter(enhancer_id %in% input$primary_value)
      options <- sort(unique(sub$tf_name))
      label <- "Select transcription factor(s) binding this enhancer:"
    } else {
      sub <- df %>% filter(enhancer_id %in% input$primary_value)
      options <- sort(unique(sub$connected_gene))
      label <- "Select connected target gene(s):"
    }
    
    selectizeInput("secondary_value", label, choices = options, multiple = TRUE)
  })
  
  # Step 5: Tertiary selection (depends on combination)
  output$enhancer_selection_ui <- renderUI({
    req(input$primary_type, input$secondary_type,
        input$primary_value, input$secondary_value)
    df <- current_df()
    
    # --- Subsets based on primary and secondary inputs ---
    if (input$primary_type == "TF") {
      base_1 <- df %>% filter(tf_name %in% input$primary_value)
    } else if (input$primary_type == "Gene") {
      base_1 <- df %>% filter(connected_gene %in% input$primary_value)
    } else {
      base_1 <- df %>% filter(enhancer_id %in% input$primary_value)
    }
    
    if (input$secondary_type == "TF") {
      base_2 <- df %>% filter(tf_name %in% input$secondary_value)
    } else if (input$secondary_type == "Gene") {
      base_2 <- df %>% filter(connected_gene %in% input$secondary_value)
    } else {
      base_2 <- df %>% filter(enhancer_id %in% input$secondary_value)
    }
    
    # --- Conditional logic: determine connecting elements ---
    if ((input$primary_type == "TF" && input$secondary_type == "Gene") ||
        (input$primary_type == "Gene" && input$secondary_type == "TF")) {
      shared <- intersect(base_1$enhancer_id, base_2$enhancer_id)
      label <- "Step 3️⃣ Select enhancer(s) connecting both conditions:"
      
    } else if ((input$primary_type == "Gene" && input$secondary_type == "Enhancer") ||
               (input$primary_type == "Enhancer" && input$secondary_type == "Gene")) {
      
      # Caso Gene ↔ Enhancer → mostrar TFs que conectan ambos
      shared <- intersect(base_1$tf_name, base_2$tf_name)
      label <- "Step 3️⃣ Select transcription factor(s) connecting both conditions:"
      
    } else if ((input$primary_type == "TF" && input$secondary_type == "Enhancer") ||
               (input$primary_type == "Enhancer" && input$secondary_type == "TF")) {
      
      # Caso TF ↔ Enhancer → mostrar genes regulados por ambos
      shared <- intersect(base_1$connected_gene, base_2$connected_gene)
      label <- "Step 3️⃣ Select target gene(s) connecting both conditions:"
      
    } else {
      return(NULL)
    }
    
    if (length(shared) == 0) return(NULL)
    
    selectizeInput(
      "selected_enhancers",
      label,
      choices = sort(unique(shared)),
      multiple = TRUE,
      selected = shared,
      options = list(plugins = list('remove_button'))
    )
  })
  
  
  # =========================================================
  # Add network button behavior (refactor estable)
  # =========================================================
  # This section controls how each “network context” is created, stored,
  # edited, and removed before the final visualization is generated.
  # Each context represents a specific TF–Enhancer–Gene subset defined
  # by user selections and active filter thresholds.
  # =========================================================
  
  # Reactive storage for user-added network contexts
  added_networks <- reactiveVal(list()) # stores network sets added by user
  bound_ids <- reactiveVal(character(0)) # keeps track of which edit/remove buttons are bound
  
  # Summary display of added networks
  output$added_networks_summary <- renderUI({
    sets <- added_networks()
    
    # No networks added yet
    if (length(sets) == 0) {
      return(div(
        style = "font-style: italic; color: #7f8c8d; margin-top: 5px;",
        "No networks added yet."
      ))
    }
    
    # Render summary cards for each stored network
    summaries <- lapply(seq_along(sets), function(i) {
      s <- sets[[i]]
      
      primary_str   <- if (!is.null(s$primary_value)) paste(s$primary_value, collapse = ", ") else "-"
      secondary_str <- if (!is.null(s$secondary_value)) paste(s$secondary_value, collapse = ", ") else "-"
      enh_str       <- if (!is.null(s$enhancers)) paste(unique(s$enhancers), collapse = ", ") else "-"
      tf_str        <- if (!is.null(s$tfs)) paste(unique(s$tfs), collapse = ", ") else "-"
      
      div(
        style = "margin-bottom: 8px; background: #f8f9fa; border-left: 4px solid #2c3e50; padding: 6px 10px; border-radius: 6px;",
        fluidRow(
          column(
            9,
            HTML(paste0(
              "<b>Network ", i, "</b><br>",
              "• <b>Primary (", s$primary_type, "):</b> ", primary_str, "<br>",
              "• <b>Secondary (", s$secondary_type, "):</b> ", secondary_str, "<br>",
              if (tf_str != "-") paste0("• <b>TF(s):</b> ", tf_str, "<br>") else "",
              if (enh_str != "-") paste0("• <b>Enhancer(s):</b> ", enh_str) else ""
            ))
          ),
          column(
            3,
            align = "right",
            div(
              style = "display:flex; gap:6px; justify-content:flex-end; margin-top:10px;",
              actionButton(
                inputId = paste0("edit_net_", s$id),
                label = NULL,
                icon = icon("pen"),
                class = "btn btn-info btn-sm"
              ),
              actionButton(
                inputId = paste0("remove_net_", s$id),
                label = NULL,
                icon = icon("trash"),
                class = "btn btn-danger btn-sm"
              )
            )
          )
        )
      )
    })
    
    tagList(summaries)
  })
  
  # Dynamic observers for edit/remove buttons
  observe({
    sets <- added_networks()
    ids  <- vapply(sets, function(s) s$id, FUN.VALUE = character(1), USE.NAMES = FALSE)
    new_ids <- setdiff(ids, bound_ids())
    if (length(new_ids) == 0 || is.null(new_ids)) return(invisible())
    
    # Bind new observers for each new network context
    for (id in new_ids) {
      local({
        this_id <- id
        
        # Remove network
        observeEvent(input[[paste0("remove_net_", this_id)]], {
          current <- added_networks()
          keep <- vapply(current, function(s) s$id != this_id, FUN.VALUE = logical(1))
          added_networks(current[keep])
          showNotification(paste0("🗑️ Network removed: ", this_id), type = "warning")
        }, ignoreInit = TRUE)
        
        # Edit network
        observeEvent(input[[paste0("edit_net_", this_id)]], {
          sets_now <- added_networks()
          s <- Filter(function(x) x$id == this_id, sets_now)[[1]]
          
          # Restore previous selections in UI
          updateSelectInput(session, "primary_type", selected = s$primary_type)
          session$onFlushed(function() {
            updateSelectizeInput(session, "primary_value", selected = s$primary_value)
            updateSelectInput(session, "secondary_type", selected = s$secondary_type)
            session$onFlushed(function() {
              updateSelectizeInput(session, "secondary_value", selected = s$secondary_value)
              if (!is.null(s$enhancers)) {
                updateSelectizeInput(session, "selected_enhancers", selected = s$enhancers)
              }
            }, once = TRUE)
          }, once = TRUE)
          
          showNotification(paste0("✏️ Editing network ", this_id, " — values loaded."), type = "message")
        }, ignoreInit = TRUE)
      })
    }
    
    # Update the list of bound IDs to avoid rebinding
    bound_ids(c(bound_ids(), new_ids))
  })
  
  # Add new network context
  observeEvent(input$add_network, {
    df <- df_base()
    req(df)
    
    # --- Minimal validation ---
    req(input$primary_type, input$secondary_type)
    req(input$primary_value, length(input$primary_value) > 0)
    req(input$secondary_value, length(input$secondary_value) > 0)
    
    # ===============================
    # Apply active filter thresholds
    # ===============================
    df_filt <- df %>%
      dplyr::filter(
        motif_score >= input$motif_cutoff,
        score       >= input$enhgene_cutoff,
        between(tf_l2fc, input$tf_fc_range[1], input$tf_fc_range[2]),
        between(log2FC,  input$gene_fc_range[1], input$gene_fc_range[2])
      )
    
    # Safe replacement for missing values (avoid NA issues)
    df_filt <- df_filt %>%
      mutate(
        tf_l2fc = ifelse(is.na(tf_l2fc), 0, tf_l2fc),
        log2FC  = ifelse(is.na(log2FC), 0, log2FC)
      )
    
    # ===============================
    # Compute enhancers for current set
    # ===============================
    enh_list <- character(0)
    pt <- input$primary_type
    st <- input$secondary_type
    
    # Handle all possible pair combinations
    if (pt == "TF" && st == "Gene") {
      base_tf   <- df_filt %>% dplyr::filter(tf_name %in% input$primary_value)
      base_gene <- df_filt %>% dplyr::filter(connected_gene %in% input$secondary_value)
      
      enh_common <- intersect(base_tf$enhancer_id, base_gene$enhancer_id)
      if (!is.null(input$selected_enhancers) && length(input$selected_enhancers) > 0) {
        enh_common <- intersect(enh_common, input$selected_enhancers)
      }
      enh_list <- unique(enh_common)
      
    } else if (pt == "TF" && st == "Enhancer") {
      sub <- df_filt %>%
        dplyr::filter(tf_name %in% input$primary_value,
                      enhancer_id %in% input$secondary_value)
      enh_list <- unique(sub$enhancer_id)
      
    } else if (pt == "Gene" && st == "TF") {
      sub <- df_filt %>%
        dplyr::filter(connected_gene %in% input$primary_value,
                      tf_name %in% input$secondary_value)
      enh_list <- unique(sub$enhancer_id)
      
    } else if (pt == "Gene" && st == "Enhancer") {
      sub <- df_filt %>%
        dplyr::filter(connected_gene %in% input$primary_value,
                      enhancer_id %in% input$secondary_value)
      enh_list <- unique(sub$enhancer_id)
      
    } else if (pt == "Enhancer" && st == "TF") {
      sub <- df_filt %>%
        dplyr::filter(enhancer_id %in% input$primary_value,
                      tf_name %in% input$secondary_value)
      enh_list <- unique(sub$enhancer_id)
      
    } else { # pt == "Enhancer" && st == "Gene"
      sub <- df_filt %>%
        dplyr::filter(enhancer_id %in% input$primary_value,
                      connected_gene %in% input$secondary_value)
      enh_list <- unique(sub$enhancer_id)
    }
    
    # --- No matching enhancers ---
    if (length(enh_list) == 0) {
      shinyalert::shinyalert(
        title = "⚠️ No data found",
        text  = "No enhancers matched this combination under current thresholds.",
        type  = "warning"
      )
      return()
    }
    
    # -------------------------------
    # Store new network in reactive memory
    # -------------------------------
    new_id <- paste0("net_", as.integer(as.numeric(Sys.time())*1000), "_", sample(1e6, 1))
    
    new_set <- list(
      id              = new_id,
      primary_type    = input$primary_type,
      primary_value   = input$primary_value,
      secondary_type  = input$secondary_type,
      secondary_value = input$secondary_value,
      enhancers       = enh_list,
      df_filtered     = df_filt
    )
    
    all_sets <- added_networks()
    added_networks(append(all_sets, list(new_set)))
    
    showNotification(
      paste0("✅ Network context added (", length(added_networks()), " total). ",
             "Enhancers in this set: ", length(enh_list)),
      type = "message"
    )
    
    # -------------------------------
    # Reset UI inputs for next exploration
    # -------------------------------
    updateSelectInput(session, "primary_type", selected = "TF")
    updateSelectizeInput(session, "primary_value",  selected = character(0))
    updateSelectInput(session, "secondary_type",    selected = NULL)
    updateSelectizeInput(session, "secondary_value",selected = character(0))
    if (!is.null(input$selected_enhancers)) {
      updateSelectizeInput(session, "selected_enhancers", selected = character(0))
    }
  })
  
  # =========================================================
  # generate_network_plot() — Core rendering helper
  # =========================================================
  # This helper function takes the final filtered subset of
  # interactions (`sub_all`) and:
  #   1. Extracts enhancer IDs.
  #   2. Calls generate_enhancer_network() to build the visNetwork.
  #   3. Renders both network and interaction table.
  #   4. Displays node/edge summary statistics.
  # =========================================================
  generate_network_plot <- function(sub_all) {
    # --- Prepare base dataset and enhancer list ---
    df <- df_base()
    enh_list <- unique(sub_all$enhancer_id)
    
    # --- Layout mode selection based on user toggle ---
    layout_mode <- if (isTRUE(input$use_physics)) "physics" else "kk"
    
    # --- Build visNetwork object using unified generator ---
    vis <- generate_enhancer_network(
      df = sub_all,
      enhancers_of_interest = enh_list,
      node_font_size = input$font_size,
      node_size = input$node_size,
      fc_limit = 1.5,
      layout_method = layout_mode,
      filter_edges = list(
        tf_enh = input$motif_cutoff,
        enh_gene = input$enhgene_cutoff
      )
    )
    
    # Render network
    output$enh_network <- renderVisNetwork({ vis })
    
    # Render TF–Enhancer–Gene interaction table
    output$interaction_table <- renderDT({
      sub_all <- sub_all %>%
        mutate(
          tf_l2fc = ifelse(is.na(tf_l2fc), 0, tf_l2fc),
          log2FC  = ifelse(is.na(log2FC), 0, log2FC)
        )
      
      # Define columns to display (dynamic inclusion)
      cols_to_show <- c("tf_name", "enhancer_id", "connected_gene", "motif_score", "score")
      if ("cluster" %in% colnames(sub_all)) cols_to_show <- c(cols_to_show, "cluster")
      if ("log2FC" %in% colnames(sub_all)) cols_to_show <- c(cols_to_show, "log2FC")
      if ("tf_l2fc" %in% colnames(sub_all)) cols_to_show <- c(cols_to_show, "tf_l2fc")
      
      tbl <- sub_all %>%
        select(any_of(cols_to_show)) %>%
        rename(
          "Transcription factor (TF)"        = tf_name,
          "Enhancer ID"                      = enhancer_id,
          "Target gene"                      = connected_gene,
          "Motif score (TF → Enhancer)"      = motif_score,
          "Enhancer score (Enhancer → Gene)" = score,
          "log2FC (Target gene)"             = log2FC,
          "log2FC (TF)"                      = tf_l2fc
        )
      
      datatable(
        tbl,
        rownames = FALSE,
        filter   = "top",              # allows on-the-fly filtering
        options  = list(pageLength = 10, scrollX = TRUE),
        caption  = htmltools::tags$caption(
          style = 'caption-side: top; text-align: left; font-weight:bold;',
          "TF–Enhancer–Gene interactions for current network"
        )
      )
    })
    
    # Display network summary statistics
    output$network_summary <- renderText({
      nodes <- length(unique(c(sub_all$tf_name, sub_all$connected_gene, sub_all$enhancer_id)))
      edges <- nrow(sub_all)
      paste0("🔹 Estimated nodes: ", nodes, " | 🔸 Estimated edges: ", edges)
    })
  }
  
  # =========================================================
  # Generate network based on all added contexts
  # =========================================================
  # This observer listens for the "Generate network" button.
  # It merges all user-added contexts (added_networks()),
  # applies filters, checks network size, and passes the
  # consolidated subset to generate_network_plot().
  # =========================================================
  observeEvent(input$generate, {
    df <- df_base()              # current dataset
    sets <- added_networks()     # list of added contexts
    
    # Fallback mode: use current selection if no sets added
    if (length(sets) == 0) {
      sets <- list(list(
        primary_type = input$primary_type,
        primary_value = input$primary_value,
        secondary_type = input$secondary_type,
        secondary_value = input$secondary_value
      ))
    }
    
    # Merge filtered data from all contexts
    if (all(vapply(sets, function(s) "df_filtered" %in% names(s), logical(1)))) {
      df <- do.call(dplyr::bind_rows, lapply(sets, `[[`, "df_filtered"))
    } else {
      df <- df %>%
        dplyr::filter(
          motif_score >= input$motif_cutoff,
          score       >= input$enhgene_cutoff,
          between(tf_l2fc, input$tf_fc_range[1], input$tf_fc_range[2]),
          between(log2FC,  input$gene_fc_range[1], input$gene_fc_range[2])
        )
    }
    
    # Iterate through sets to build subset for visualization
    sub_all <- tibble()
    
    for (s in sets) {
      if (is.null(s$primary_type) || is.null(s$secondary_type)) next
      
      # Step 1 — Filter by primary type
      if (s$primary_type == "TF") {
        base_tf <- df %>% dplyr::filter(tf_name %in% s$primary_value)
      } else if (s$primary_type == "Gene") {
        base_tf <- df %>% dplyr::filter(connected_gene %in% s$primary_value)
      } else {
        base_tf <- df %>% dplyr::filter(enhancer_id %in% s$primary_value)
      }
      
      # Step 2 — Filter by secondary type
      if (s$secondary_type == "TF") {
        base_gene <- df %>% dplyr::filter(tf_name %in% s$secondary_value)
      } else if (s$secondary_type == "Gene") {
        base_gene <- df %>% dplyr::filter(connected_gene %in% s$secondary_value)
      } else {
        base_gene <- df %>% dplyr::filter(enhancer_id %in% s$secondary_value)
      }
      
      # Step 3 — Identify enhancers shared between both sets
      enh_common <- intersect(base_tf$enhancer_id, base_gene$enhancer_id)
      if (length(enh_common) == 0) next
      
      sub <- df %>% dplyr::filter(enhancer_id %in% enh_common)
      sub_all <- dplyr::bind_rows(sub_all, sub)
    }
    
    sub_all <- sub_all %>%
      distinct(tf_name, enhancer_id, connected_gene, .keep_all = TRUE)
    
    # Step 4: Safety check for large networks
    if (nrow(sub_all) > 1000) {
      shinyalert::shinyalert(
        title = "⚠️ Large network detected",
        text = paste0(
          "This selection will generate a large network with ",
          nrow(sub_all), " interactions.\n\n",
          "Rendering may take several seconds and affect responsiveness.\n\n",
          "Would you like to continue?"
        ),
        type = "warning",
        showCancelButton = TRUE,
        confirmButtonText = "Yes, generate network",
        cancelButtonText = "Cancel",
        callbackR = function(x) {
          if (!isTRUE(x)) {
            showNotification("🚫 Network generation cancelled.", type = "message")
            return()
          }
          # Continue only if confirmed
          isolate({
            generate_network_plot(sub_all)
          })
        }
      )
      return()  # prevent automatic continuation
    }
    
    # -------------------------------------------------------
    # Stop if empty
    # -------------------------------------------------------
    # If no valid interactions were found after all filters
    # and selections, stop execution and notify the user.
    # This prevents empty network or table renders.
    if (nrow(sub_all) == 0) {
      output$interaction_table <- renderDT(NULL)
      output$network_summary <- renderText("Estimated nodes: 0 | Estimated edges: 0")
      shinyalert::shinyalert(
        title = "⚠️ No data found",
        text = "No interactions found for this combination. Try adjusting filters or selecting another element.",
        type = "warning"
      )
      return() # early exit
    }
    
    # -------------------------------------------------------
    # Generate network visualization
    # -------------------------------------------------------
    # Extract unique enhancer IDs to feed the network generator
    enh_list <- unique(sub_all$enhancer_id)
    
    # -------------------------------------------------------
    # Layout mode selection (auto-switch for large networks)
    # -------------------------------------------------------
    # If the user has physics mode enabled but the network is very large,
    # automatically switch to the static KK layout for performance stability.
    if (isTRUE(input$use_physics) && nrow(sub_all) > 1000) {
      showNotification("⚠️ Large network — switching to static layout for stability.", type = "warning")
      layout_mode <- "kk"
    } else {
      layout_mode <- if (isTRUE(input$use_physics)) "physics" else "kk"
    }
    
    # -------------------------------------------------------
    # Generate visNetwork object
    # -------------------------------------------------------
    # Build the interactive TF–Enhancer–Gene network visualization
    vis <- generate_enhancer_network(
      df = sub_all,
      enhancers_of_interest = enh_list,
      node_font_size = input$font_size,
      node_size = input$node_size,
      fc_limit = 1.5,
      layout_method = layout_mode,
      filter_edges = list(
        tf_enh = input$motif_cutoff,
        enh_gene = input$enhgene_cutoff
      )
    )
    
    # Render the visNetwork in the main UI panel
    output$enh_network <- renderVisNetwork({ vis })
    
    # -------------------------------------------------------
    # Render interaction table (safe re-render)
    # -------------------------------------------------------
    # Build the detailed table of all TF–Enhancer–Gene relationships
    # associated with the current network. Includes optional columns
    # like cluster, log2FC for genes, and log2FC for TFs.
    output$interaction_table <- renderDT({
      sub_all <- sub_all %>%
        mutate(
          tf_l2fc = ifelse(is.na(tf_l2fc), 0, tf_l2fc),
          log2FC  = ifelse(is.na(log2FC), 0, log2FC)
        )
      
      cols_to_show <- c("tf_name", "enhancer_id", "connected_gene", "motif_score", "score")
      if ("cluster" %in% colnames(sub_all)) cols_to_show <- c(cols_to_show, "cluster")
      if ("log2FC" %in% colnames(sub_all)) cols_to_show <- c(cols_to_show, "log2FC")
      if ("tf_l2fc" %in% colnames(sub_all)) cols_to_show <- c(cols_to_show, "tf_l2fc")
      
      tbl <- sub_all %>%
        select(any_of(cols_to_show)) %>%
        rename(
          "Transcription factor (TF)"        = tf_name,
          "Enhancer ID"                      = enhancer_id,
          "Target gene"                      = connected_gene,
          "Motif score (TF → Enhancer)"      = motif_score,
          "Enhancer score (Enhancer → Gene)" = score,
          "log2FC (Target gene)"             = log2FC,
          "log2FC (TF)"                      = tf_l2fc
        )
      
      datatable(
        tbl,
        rownames = FALSE,
        filter   = "top",
        options  = list(pageLength = 10, scrollX = TRUE),
        caption  = htmltools::tags$caption(
          style = 'caption-side: top; text-align: left; font-weight:bold;',
          "TF–Enhancer–Gene interactions for current network"
        )
      )
    })
    
    # -------------------------------------------------------
    # Network summary
    # -------------------------------------------------------
    # Display the estimated number of nodes and edges currently
    # visualized, shown above the network in the main panel.
    output$network_summary <- renderText({
      nodes <- length(unique(c(sub_all$tf_name, sub_all$connected_gene, sub_all$enhancer_id)))
      edges <- nrow(sub_all)
      paste0("🔹 Estimated nodes: ", nodes, " | 🔸 Estimated edges: ", edges)
    })
  })
  
  
  # =========================================================
  # Clear All (full reset)
  # =========================================================
  # This observer performs a full reset of the app state when
  # the user clicks the "Clear all" button.
  #
  # Workflow:
  #   1. Reset all reactive values and lists.
  #   2. Clear visual outputs (network, table, summary).
  #   3. Reset all filters, sliders, and select inputs to defaults.
  #   4. Handle optional UI elements dynamically (if present).
  #   5. Restore the empty "No networks added yet" message.
  #   6. Notify the user of successful cleanup.
  # =========================================================
  observeEvent(input$clear_all, {
    # 1. Reset main reactive storage
    added_networks(list()) # clears all added networks in memory
    
    # 2. Clear outputs: network visualization, table, and summary text
    output$enh_network <- renderVisNetwork(NULL)
    output$interaction_table <- renderDT(NULL)
    output$network_summary <- renderText("Estimated nodes: 0 | Estimated edges: 0")
    
    # 3. Reset core selection inputs to default state
    updateSelectInput(session, "primary_type", selected = "TF")
    updateSelectizeInput(session, "primary_value", selected = character(0))
    updateSelectInput(session, "secondary_type", selected = NULL)
    updateSelectizeInput(session, "secondary_value", selected = character(0))
    
    # 4. Conditional reset of optional UI elements
    # These inputs are only present in specific app modes,
    # so check their existence before attempting to reset.
    if ("selected_enhancers" %in% names(input)) {
      updateSelectizeInput(session, "selected_enhancers", selected = character(0))
    }
    if ("add_mode" %in% names(input)) {
      updateRadioButtons(session, "add_mode", selected = character(0))
    }
    if ("add_elements" %in% names(input)) {
      updateSelectizeInput(session, "add_elements", selected = character(0))
    }
    
    # 5. Reset numeric sliders to base values
    updateSliderInput(session, "motif_cutoff", value = 0)
    updateSliderInput(session, "enhgene_cutoff", value = 0)
    updateSliderInput(session, "font_size", value = 35)
    
    # 6. Restore sidebar summary to default message
    if ("added_networks_summary" %in% names(output)) {
      output$added_networks_summary <- renderUI({
        div(style = "font-style: italic; color: #7f8c8d; margin-top: 5px;",
            "No networks added yet.")
      })
    }
    
    # 7. User feedback notification
    showNotification("♻️ All selections, networks, and filters cleared.", type = "message")
  })
  # =========================================================
  # Snapshot button — Capture current visNetwork view
  # =========================================================
  # Sends a custom message to the client (JavaScript side)
  # to trigger an image capture of the current visNetwork graph.
  # The JS handler is defined in the UI or www/scripts.js file.
  # =========================================================
  observeEvent(input$snapshot_btn, {
    session$sendCustomMessage("captureVisNetwork", list())
  })
  
  # =========================================================
  # Current dataset label
  # =========================================================
  # Displays the name of the active dataset in use.
  # Uses a human-friendly label if available, otherwise
  # falls back to the raw dataset identifier.
  # =========================================================
  output$current_dataset_label <- renderText({
    req(input$dataset_choice)
    # Resolve user-friendly label if defined in `dataset_labels`
    label <- if (input$dataset_choice %in% names(dataset_labels)) {
      dataset_labels[[input$dataset_choice]]
    } else {
      input$dataset_choice
    }
    paste0("Active dataset: ", label)
  })
  
} # <-- End of server function

# =========================================================
# Run Shiny App
# =========================================================
# This final command launches the Enhancer Network Explorer.
# It binds the user interface (ui) and the server logic (server)
# defined above into a single reactive Shiny application.
#
# The app dynamically:
#   • Detects available datasets in the working directory.
#   • Allows guided multi-step filtering and visualization.
#   • Renders TF–Enhancer–Gene interaction networks interactively.
#
# To run locally:
#   - Ensure all dependencies are installed (see DESCRIPTION file).
#   - Confirm that .RData or .rds datasets are in the same directory.
#   - Execute this script directly from RStudio or the terminal.
#
# Once launched, open the URL provided in the R console
# (typically http://127.0.0.1:XXXX) to interact with the app on 
# internet explorer of your preference (Tested in Chrome/Firefox/Safari).
# =========================================================

shinyApp(ui, server)

