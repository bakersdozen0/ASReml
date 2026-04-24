# ==============================================================================
# ASREML DATA ADJUSTER & XML GENERATOR (Standalone)
# ==============================================================================

library(tidyverse)
library(xml2) 

# 0. CONTROL PANEL #### 
trial_folder  <- "C:/Users/james.baker/Forest Research/TW CBC-TBA-NextGenBritishConifers - Share/Sitka/Backwards Selected Fullsib P96-P99 experiments/Kintyre 18"
csv_file      <- "Kintyre_18_S.csv"
dms_xml_file  <- "Kintyre_18_DMS2.xml" # The master parent XML from Dataplan

# ---> SPECIFY WHAT YOU WANT TO ADJUST HERE <---
target_run_id <- "Edge_correction"  # The folder to read solutions from
target_trait  <- "Dbhob_28"         # The specific trait to adjust
target_model  <- "Spatial AR1"      # "Design", "Design+", or "Spatial AR1"
# ----------------------------------------------

# 1. SETUP PATHS & READ RAW DATA ####
run_dir <- file.path(trial_folder, "Analyses", target_run_id)
out_dir <- file.path(trial_folder, "Analyses") # Save to master folder

if(!dir.exists(run_dir)) stop("Cannot find the specified run folder: ", run_dir)
if(!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw_data <- read.csv(file.path(trial_folder, csv_file), stringsAsFactors = FALSE)
raw_data <- raw_data %>% mutate(Record = row_number()) 

# Verify column exists
if (!(target_trait %in% colnames(raw_data))) {
  stop(paste0("Trait '", target_trait, "' not found in CSV. Check casing!"))
}

# Ensure trait is numeric
raw_data[[target_trait]] <- suppressWarnings(as.numeric(as.character(raw_data[[target_trait]])))

# Define naming convention: Adj_[trait] and Surface_Sum_[trait]
adj_col     <- paste0("Adj_", target_trait)
ss_col      <- paste0("Surface_Sum_", target_trait)

# 2. LOCATE MODEL OUTPUTS ####
sln_file <- file.path(run_dir, paste0(target_trait, "_", target_model, ".sln"))
yht_file <- file.path(run_dir, paste0(target_trait, "_", target_model, ".yht"))

if(!file.exists(sln_file)) stop("Cannot find .sln file: ", sln_file)

cat("Processing Adjustments for:", target_trait, "\n")

# 3. EXTRACT SOLUTIONS ####
sln <- read.table(sln_file, skip = 1, fill = TRUE, stringsAsFactors = FALSE, 
                  colClasses = c("character", "character", "numeric", "numeric"))
colnames(sln) <- c("Term", "Level", "Estimate", "SE")

b_eff <- sln %>% filter(Term == "Block") %>% select(Level, Estimate)
r_eff <- sln %>% filter(Term == "Block.Prow") %>% select(Level, Estimate)
p_eff <- sln %>% filter(Term == "Block.Ppos") %>% select(Level, Estimate)

adj_data <- raw_data %>% 
  mutate(
    Block_key = as.character(Block),
    BProw_key = sprintf("%d.%03d", as.integer(Block), as.integer(Prow)),
    BPpos_key = sprintf("%d.%03d", as.integer(Block), as.integer(Ppos))
  ) %>%
  left_join(b_eff, by = c("Block_key" = "Level")) %>% rename(Block_Est = Estimate) %>%
  left_join(r_eff, by = c("BProw_key" = "Level")) %>% rename(Row_Est = Estimate) %>%
  left_join(p_eff, by = c("BPpos_key" = "Level")) %>% rename(Col_Est = Estimate)

# 4. CONDITIONAL MATH BASED ON MODEL ####
if (target_model == "Spatial AR1") {
  if(!file.exists(yht_file)) stop("Spatial adjustment requires .yht file.")
  yht <- read.table(yht_file, skip = 1)
  colnames(yht) <- c("Record", "Yhat", "Residual", "Hat")
  
  adj_data <- adj_data %>%
    left_join(yht, by = "Record") %>%
    mutate(
      Design_Sum = replace_na(Block_Est, 0) + replace_na(Row_Est, 0) + replace_na(Col_Est, 0),
      Local_Trend = replace_na(Residual, 0)
    )
} else {
  adj_data <- adj_data %>%
    mutate(
      Design_Sum = replace_na(Block_Est, 0) + replace_na(Row_Est, 0) + replace_na(Col_Est, 0),
      Local_Trend = 0
    )
}

# --- TRANSFORMATION & ROW OMISSION ---
# 1. Calculate the values
# 2. Filter out any rows where the tree was not measured
# 3. Rename columns to the final Dataplan-compliant names
adj_data <- adj_data %>%
  mutate(
    temp_ss = Design_Sum + Local_Trend,
    temp_adj = .data[[target_trait]] - temp_ss
  ) %>%
  filter(!is.na(.data[[target_trait]])) %>% 
  rename(!!ss_col := temp_ss, !!adj_col := temp_adj)

# 5. EXPORT DATA (Surgical 3-Column Export) ####
final_export <- adj_data %>%
  select(
    Stem_id, 
    all_of(ss_col), 
    all_of(adj_col)
  )

# CSV naming: Adj_DBHOB_28.csv
out_name_base <- paste0("Adj_", target_trait)
write.csv(final_export, file.path(out_dir, paste0(out_name_base, ".csv")), row.names = FALSE, na = "")

cat("SUCCESS! CSV saved with unmeasured lines omitted.\n")

# 6. PARSE MASTER XML & GENERATE NEW XML ####
base_xml_path <- file.path(trial_folder, dms_xml_file)
orig_group <- "Unknown"; orig_units <- ""; orig_desc  <- target_trait

if (file.exists(base_xml_path)) {
  doc <- read_xml(base_xml_path)
  traits <- xml_find_all(doc, "//trait")
  trait_names <- xml_attr(traits, "name")
  match_idx <- which(tolower(trait_names) == tolower(target_trait))
  
  if (length(match_idx) > 0) {
    t_node <- traits[[match_idx[1]]]
    orig_group <- xml_attr(t_node, "group_name")
    orig_units <- xml_attr(t_node, "units")
    orig_desc  <- xml_attr(t_node, "description")
  }
}

adj_desc <- paste0("Surface Sum Adjusted (", target_model, ") ", orig_desc)
ss_desc  <- paste0("Surface Sum (", target_model, ") ", orig_desc)

xml_content <- paste0(
  '<?xml version="1.0" encoding="UTF-8" ?>\n',
  '<trial><traitlist>\n',
  '<trait group_name="', orig_group, '" name="', adj_col, '" trait_type="M" description="', adj_desc, '" data_type="N" units="', orig_units, '" validate="none" is_solver_mappable="1" ></trait>\n',
  '<trait group_name="', orig_group, '" name="', ss_col, '" trait_type="M" description="', ss_desc, '" data_type="N" units="', orig_units, '" validate="none" is_solver_mappable="0" ></trait>\n',
  '</traitlist></trial>'
)

writeLines(xml_content, file.path(out_dir, paste0(out_name_base, ".XML")))

cat("SUCCESS! Minimal CSV and XML saved to master Analyses folder.\n")