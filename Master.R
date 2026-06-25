# ==============================================================================
# THE WHEEL: CONTROL PANEL & SCRAPER INITIALIZATION
# ==============================================================================

# 1. PATHS & DIRECTORIES ####
trial_folder  <- "C:/Users/james.baker/Forest Research/TW CBC-TBA-NextGenBritishConifers - Share/Sitka/High GCA Fullsib P85-P87 experiments/Brecon 15"
project_name  <- "Brecon_15_S"
run_id        <- "Baseline"

# Path to ASReml Standalone
asreml_path   <- "C:/Program Files/ASReml4/bin/asreml.exe"

# 2. FILE NAMES COLLATING ####
as_file       <- paste0(project_name, ".as")
csv_file      <- paste0(project_name, ".csv")
dms_xml_file  <- paste0(project_name, ".xml") 

# 3. CORE LIBRARIES ####
library(here)
library(tidyverse)
library(xml2)
library(ggplot2)
library(patchwork)
library(svglite)
library(flextable)
library(officer)

# 4. ADJUSTMENT ENGINE CONTROLS ####
# (Only used if running the Spatial Adjustment Engine)
target_trait  <- "Pil_16"    
target_model  <- "Design+"   # "Design", "Design+", or "Spatial AR1"

# ==============================================================================
# 5. EXECUTE THE ENGINES
# Run the engine you want to run for this session.
# ==============================================================================
source(here("ASReml_univariate_loops.R"))


source(here("spatial_adjustments.R"))