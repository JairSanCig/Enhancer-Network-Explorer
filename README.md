# Enhancer-Network-Explorer

Version 2.9 — Stable Release
Author: Juan Jair Santillán-Cigales
Institution: Instituto de Fisiología Celular (IFC-UNAM)
Date: 23 October 2025

# Overview

Enhancer Network Explorer is an interactive Shiny application designed to visualize and explore transcriptional regulatory networks between Transcription Factors (TFs), Enhancers, and Target Genes.

This tool enables researchers to dynamically select datasets, filter interactions by threshold values, and generate publication-ready network visualizations.

# Features

Automatic dataset detection: identifies valid .rds or .RData files containing enhancer–TF–gene relationships.

Guided multi-step search: interactively select TFs, genes, or enhancers to explore regulatory connections.

Filter controls: adjust thresholds for motif scores, enhancer–gene scores, and log₂FC ranges.

Interactive visualization: networks rendered with visNetwork and enhanced by igraph.

Multi-context analysis: combine or compare multiple network contexts before visualization.

Dynamic tables: interactive filtering and sorting of TF–Enhancer–Gene interactions using DT.

Snapshot export: export network views as PNG files.

# Structure
Component Description
App.R Main Shiny app (UI + Server combined).
get_available_datasets()  Detects and validates datasets loaded into the environment.
generate_enhancer_network() Core network-building function using visNetwork.
UI Section  Defines layout, controls, and visualization panels.
Server Section  Handles reactivity, filtering, and data logic.
/data/  Folder containing example .rds datasets.
/www/ (Optional) Folder for icons, CSS, or external JS resources.

# Installation
1. Requirements

Tested under:

R 4.4.1

Shiny 1.9.1

2. Dependencies

Install missing libraries with:

install.packages(c(
  "shiny", "shinyWidgets", "shinyalert", "shinyjqui",
  "visNetwork", "dplyr", "tibble", "scales",
  "igraph", "htmlwidgets", "DT", "bslib"
))

3. Run the App

Place your datasets (e.g. df_filtered_d8_MERGED_ENCODE_JASPAR.rds, df_filtered_d16_MERGED_ENCODE_JASPAR.rds) inside the data/ directory or in the same folder as App.R.

Open RStudio or your terminal in that directory.

Run:

shiny::runApp("path/to/Enhancer_Network_Explorer")


Once launched, open the local URL shown in your R console (typically http://127.0.0.1:XXXX) to access the app in your browser.

# Input Data Requirements

Each dataset must include at least the following columns:

Column  Description
enhancer_id Unique enhancer identifier
tf_name Transcription factor name
connected_gene  Associated target gene
motif_score TF–Enhancer binding score
score Enhancer–Gene linkage score
log2FC  Target gene log₂ fold change (optional)
tf_l2fc TF log₂ fold change (optional)
cluster Cluster or category label (optional)

# Example Datasets

df_filtered_d8_MERGED_ENCODE_JASPAR.rds — MBO Day 8 (BHB vs Ctrl)

df_filtered_d16_MERGED_ENCODE_JASPAR.rds — MBO Day 16 (BHB vs Ctrl)

Additional datasets can be loaded dynamically;
their names will appear automatically in the dataset selection menu.

# License

This project is licensed under the
Creative Commons Attribution–NonCommercial 4.0 International (CC BY-NC 4.0) license.

You are free to share and adapt this material for academic and non-commercial purposes, provided that proper credit is given to:
Juan Jair Santillán-Cigales, Instituto de Fisiología Celular (IFC-UNAM)

Full license text:
https://creativecommons.org/licenses/by-nc/4.0/

# Citation

If you use this tool in your research, please cite:

Santillán-Cigales JJ (2025). Enhancer Network Explorer: An interactive Shiny framework for TF–Enhancer–Gene network visualization and analysis. Version 2.9. Instituto de Fisiología Celular (IFC-UNAM).
This release corresponds to the version of the app described in the manuscript currently under review. DOI:

# Contact

For questions, feedback, or collaborations:
Juan Jair Santillán-Cigales:
juan.jair.santillan@ifc.unam.mx | j.sancig@gmail.com
