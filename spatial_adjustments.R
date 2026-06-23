# ==============================================================================
# ENGINE 2: SPATIAL ADJUSTMENTS & DATAPLAN EXPORTS
# ==============================================================================

# 1. SETUP PATHS & READ RAW DATA ####
run_dir <- file.path(trial_folder, "Analyses", run_id)
out_dir <- file.path(trial_folder, "Analyses") 

if(!dir.exists(run_dir)) stop("Cannot find the specified run folder: ", run_dir)

raw_data <- read.csv(file.path(trial_folder, csv_file), stringsAsFactors = FALSE) %>% mutate(Record = row_number()) 

# Verify column exists
if (!(target_trait %in% colnames(raw_data))) stop(paste0("Trait '", target_trait, "' not found in CSV. Check casing!"))
raw_data[[target_trait]] <- suppressWarnings(as.numeric(as.character(raw_data[[target_trait]])))

# Define naming convention: Adj_[trait] and Surface_Sum_[trait]
adj_col <- paste0("Adj_", target_trait)
ss_col  <- paste0("Surface_Sum_", target_trait)

# 2. LOCATE MODEL OUTPUTS ####
sln_file <- file.path(run_dir, paste0(target_trait, "_", target_model, ".sln"))
yht_file <- file.path(run_dir, paste0(target_trait, "_", target_model, ".yht"))
if(!file.exists(sln_file)) stop("Cannot find .sln file: ", sln_file)

cat("Processing Adjustments for:", target_trait, "\n")

# 3. EXTRACT SOLUTIONS ####
sln <- read.table(sln_file, skip = 1, fill = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric"))
colnames(sln) <- c("Term", "Level", "Estimate", "SE")

# Use tolower() to protect against ASReml's random capitalizations
b_eff <- sln %>% filter(tolower(Term) == "block") %>% select(Level, Estimate)
r_eff <- sln %>% filter(tolower(Term) == "block.prow") %>% select(Level, Estimate)
c_eff <- sln %>% filter(tolower(Term) == "block.ppos") %>% select(Level, Estimate)
p_eff <- sln %>% filter(tolower(Term) %in% c("plot", "block.plot")) %>% select(Level, Estimate)

plot_col_name <- grep("^plot$", colnames(raw_data), ignore.case = TRUE, value = TRUE)

adj_data <- raw_data %>% 
  mutate(
    Block_key = as.character(Block),
    BProw_key = sprintf("%d.%03d", suppressWarnings(as.integer(Block)), suppressWarnings(as.integer(Prow))),
    BPpos_key = sprintf("%d.%03d", suppressWarnings(as.integer(Block)), suppressWarnings(as.integer(Ppos))),
    Plot_key  = if(length(plot_col_name) > 0) as.character(.data[[plot_col_name[1]]]) else NA
  ) %>%
  left_join(b_eff, by = c("Block_key" = "Level")) %>% rename(Block_Est = Estimate) %>%
  left_join(r_eff, by = c("BProw_key" = "Level")) %>% rename(Row_Est = Estimate) %>%
  left_join(c_eff, by = c("BPpos_key" = "Level")) %>% rename(Col_Est = Estimate) %>%
  left_join(p_eff, by = c("Plot_key" = "Level")) %>% rename(Plot_Est = Estimate)

# --- DYNAMIC EDGE & COVARIATE ADDER ---
adj_data$Extra_Env_Sum <- 0 # Initialize a baseline 0 column

for (e_term in c("Edge", "Distright", "Distleft", "Distedge")) {
  term_idx <- tolower(sln$Term) == tolower(e_term)
  
  # Only attempt the math if the term is in the model AND the raw data
  if (any(term_idx) && e_term %in% colnames(adj_data)) {
    t_sln <- sln[term_idx, ]
    
    if (nrow(t_sln) > 1) {
      # Scenario A: It's a Factor (e.g., Edge 1 and Edge 2)
      temp_join <- adj_data %>% mutate(Join_Key = as.character(.data[[e_term]])) %>% left_join(t_sln %>% select(Level, Estimate), by = c("Join_Key" = "Level"))
      adj_data$Extra_Env_Sum <- adj_data$Extra_Env_Sum + replace_na(temp_join$Estimate, 0)
      cat(paste("  -> Included FACTOR edge effect:", e_term, "\n"))
      
    } else if (nrow(t_sln) == 1) {
      # Scenario B: It's a Covariate Slope (e.g., continuous distance integers)
      adj_data$Extra_Env_Sum <- adj_data$Extra_Env_Sum + replace_na(suppressWarnings(as.numeric(as.character(adj_data[[e_term]]))) * t_sln$Estimate[1], 0)
      cat(paste("  -> Included COVARIATE edge effect:", e_term, "\n"))
    }
  }
}

# 4. CONDITIONAL MATH BASED ON MODEL ####
if (target_model == "Spatial AR1") {
  if(!file.exists(yht_file)) stop("Spatial adjustment requires .yht file.")
  yht <- read.table(yht_file, skip = 1); colnames(yht) <- c("Record", "Yhat", "Residual", "Hat")
  
  adj_data <- adj_data %>% 
    left_join(yht, by = "Record") %>% 
    mutate(
      Design_Sum = replace_na(Block_Est, 0) + replace_na(Row_Est, 0) + replace_na(Col_Est, 0) + replace_na(Plot_Est, 0) + Extra_Env_Sum, 
      Local_Trend = replace_na(Residual, 0)
    )
} else {
  adj_data <- adj_data %>% 
    mutate(
      Design_Sum = replace_na(Block_Est, 0) + replace_na(Row_Est, 0) + replace_na(Col_Est, 0) + replace_na(Plot_Est, 0) + Extra_Env_Sum, 
      Local_Trend = 0
    )
}

# --- TRANSFORMATION & ROW OMISSION ---
# 1. Calculate the values
# 2. Filter out any rows where the tree was not measured
# 3. Rename columns to the final Dataplan-compliant names
adj_data <- adj_data %>% 
  mutate(temp_ss = Design_Sum + Local_Trend, temp_adj = .data[[target_trait]] - temp_ss) %>%
  filter(!is.na(.data[[target_trait]])) %>% rename(!!ss_col := temp_ss, !!adj_col := temp_adj)

# 5. EXPORT DATA (Surgical 3-Column Export) ####
out_name_base <- paste0("Adj_", target_trait)
write.csv(adj_data %>% select(Stem_id, all_of(ss_col), all_of(adj_col)), file.path(out_dir, paste0(out_name_base, ".csv")), row.names = FALSE, na = "")

# 6. PARSE MASTER XML & GENERATE NEW XML ####
base_xml_path <- file.path(trial_folder, dms_xml_file)
orig_group <- "Unknown"; orig_units <- ""; orig_desc  <- target_trait

if (file.exists(base_xml_path)) {
  doc <- read_xml(base_xml_path)
  traits <- xml_find_all(doc, "//trait")
  match_idx <- which(tolower(xml_attr(traits, "name")) == tolower(target_trait))
  if (length(match_idx) > 0) { 
    t_node <- traits[[match_idx[1]]]
    orig_group <- xml_attr(t_node, "group_name")
    orig_units <- xml_attr(t_node, "units")
    orig_desc <- xml_attr(t_node, "description") 
  }
}

writeLines(paste0('<?xml version="1.0" encoding="UTF-8" ?>\n<trial><traitlist>\n<trait group_name="', orig_group, '" name="', adj_col, '" trait_type="M" description="Surface Sum Adjusted (', target_model, ') ', orig_desc, '" data_type="N" units="', orig_units, '" validate="none" is_solver_mappable="1" ></trait>\n<trait group_name="', orig_group, '" name="', ss_col, '" trait_type="M" description="Surface Sum (', target_model, ') ', orig_desc, '" data_type="N" units="', orig_units, '" validate="none" is_solver_mappable="0" ></trait>\n</traitlist></trial>'), file.path(out_dir, paste0(out_name_base, ".XML")))

cat("SUCCESS! Minimal CSV and XML saved to master Analyses folder.\n")

# --- 7. AUDIT LOGGER ---
log_file <- file.path(out_dir, "Adjusted_Traits_Manifest.txt")
write(sprintf("Trait: %-15s | Model Used: %-15s | Run ID: %-15s | Timestamp: %s", target_trait, target_model, run_id, format(Sys.time(), "%Y-%m-%d %H:%M")), file = log_file, append = TRUE)

cat("\n========================================\n      CURRENT ADJUSTMENT MANIFEST       \n========================================\n")
cat(readLines(log_file), sep = "\n")
cat("========================================\n")
