# global.R

# Load necessary libraries
library(shiny)
library(shinyWidgets)
library(shinyBS)
library(ggplot2)
library(dplyr)
library(purrr)
library(patchwork)
library(visNetwork)
library(igraph)
library(DT)
library(jsonlite)
library(plotly)

# Source the dataset registry first — provides CIRC_DATASETS, CIRC_TISSUE_CHOICES,
# and CIRC_DATASET_CHOICES_MAP used by all modules below.
source("datasets.R")

# Source the functions file
source("functions.R")
source("explorer_module.R")
source("comparison_functions.R")
source("comparison_module.R")