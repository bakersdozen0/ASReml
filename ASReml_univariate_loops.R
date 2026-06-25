# ==============================================================================
# ENGINE 1: UNIVARIATE LOOPS & MODEL DIAGNOSTICS
# ==============================================================================

# --- A. INITIAL COMPONENT VERIFICATION ---
asreml_exe <- shQuote(asreml_path) 
raw_data <- read.csv(file.path(trial_folder, csv_file), stringsAsFactors = FALSE)

# Lock in the 'data file order' exactly as the PPGVal Sfile step output it.
# This ensures ASReml's un-indexed .yht residuals will always join flawlessly.
raw_data <- raw_data %>% mutate(Record = row_number())

# Legend-fixing custom helper functions
force_breaks <- function(x) { 
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  if (!is.finite(min_val) || !is.finite(max_val)) return(c(0, 1)) # Failsafe for all NAs
  seq(min_val, max_val, length.out = 5) 
}
format_labels <- function(x) { round(x, 2) }

# --- B. DYNAMIC TRAIT SCRAPER & SURVIVAL BLOCK FILTER ---
as_text <- paste(readLines(file.path(trial_folder, as_file)), collapse = " ")
found_traits <- unique(unlist(str_extract_all(as_text, "\\b[A-Za-z0-9]+_[0-9]+\\b")))
traits_to_test <- found_traits[found_traits %in% colnames(raw_data)]

# FILTER OUT SURVIVAL TRAITS
# The "^" ensures it only looks at the start of the string
traits_to_test <- traits_to_test[!grepl("^Sur_", traits_to_test, ignore.case = TRUE)] 

cat("Automated Discovery: Found", length(found_traits), "potential matches.\n")
cat("Guard Rail: Proceeding with", length(traits_to_test), "traits (Survival excluded).\n")

# --- C. SANDBOX DIRECTORY PREPARATION ---
sandbox_dir <- file.path(trial_folder, "Analyses", run_id)
dir.create(sandbox_dir, recursive = TRUE, showWarnings = FALSE)
cat("Running scenario:", run_id, "| Sandbox:", sandbox_dir, "\n")

# Copy the CSV to the sandbox (Bypasses Fortran string limits)
file.copy(file.path(trial_folder, csv_file), file.path(sandbox_dir, csv_file), overwrite = TRUE)

# Read the master .as file, strip the messy file paths, and save to Sandbox
as_lines <- readLines(file.path(trial_folder, as_file), warn = FALSE)
csv_line_index <- grep("\\.csv", as_lines, ignore.case = TRUE)
if(length(csv_line_index) > 0) {
  qualifiers <- sub(".*\\.csv\\s*", "", as_lines[csv_line_index[1]], ignore.case = TRUE)
  as_lines[csv_line_index[1]] <- paste(csv_file, qualifiers) # Forces local read
}
writeLines(as_lines, file.path(sandbox_dir, as_file))

# Step into the Sandbox for the rest of the script
setwd(sandbox_dir)
out_dir <- sandbox_dir
models_to_run <- c("1" = "Design", "2" = "Design+", "3" = "Spatial AR1")

# --- D. VISUAL LAYOUT MAPPING GEOMETRY ---
# Calculate the geometric center for labels
block_labels <- raw_data %>%
  filter(!is.na(Ppos) & !is.na(Prow) & !is.na(Block) & 
           as.character(Block) != "0" & as.character(Block) != "" & as.character(Block) != ".") %>% 
  group_by(Block) %>%
  summarize(x_mid = mean(as.numeric(Ppos), na.rm = TRUE), y_mid = mean(as.numeric(Prow), na.rm = TRUE))

# Calculate Internal Segments
h_segments <- raw_data %>% arrange(Prow, Ppos) %>% mutate(next_block = lead(Block)) %>%
  filter(Block != next_block & Prow == lead(Prow) & !is.na(Block)) %>%
  mutate(x = Ppos + 0.5, xend = Ppos + 0.5, y = Prow - 0.5, yend = Prow + 0.5)

v_segments <- raw_data %>% arrange(Ppos, Prow) %>% mutate(next_block = lead(Block)) %>%
  filter(Block != next_block & Ppos == lead(Ppos) & !is.na(Block)) %>%
  mutate(x = Ppos - 0.5, xend = Ppos + 0.5, y = Prow + 0.5, yend = Prow + 0.5)

# Calculate Outer Trial Border
trial_border <- data.frame(
  xmin = min(raw_data$Ppos, na.rm = TRUE) - 0.5, xmax = max(raw_data$Ppos, na.rm = TRUE) + 0.5,
  ymin = min(raw_data$Prow, na.rm = TRUE) - 0.5, ymax = max(raw_data$Prow, na.rm = TRUE) + 0.5
)

# Standard Theme
theme_trial <- theme_void() + theme(plot.title = element_text(size = 20, face = "bold", margin = margin(b = 10)), plot.subtitle = element_text(size = 14, color = "grey30"), legend.key.height = unit(1.5, "cm"), legend.text = element_text(size = 12))
shared_layers <- list(geom_segment(data = h_segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 0.8, inherit.aes = FALSE), geom_segment(data = v_segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 0.8, inherit.aes = FALSE), geom_rect(data = trial_border, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=1.5, inherit.aes=FALSE), geom_label(data = block_labels, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 6, fontface = "bold", fill = alpha("white", 0.9), label.size = NA))
ppgmap_colors <- c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred")

master_results_list <- list()
master_fixed_list <- list()

# --- E. EXECUTION TRAIT LOOP RUNNING UNIVARIATES ---
for (trait in traits_to_test) {
  cat("\n========================================\nProcessing:", trait, "\n")
  raw_data[[trait]] <- suppressWarnings(as.numeric(as.character(raw_data[[trait]])))
  
  trait_variance_list <- list()
  
  # Initialize plots for this trait
  blank_plot <- ggplot() + theme_void() + annotate("text", x = 0.5, y = 0.5, label = "Model Did Not Converge", color = "darkred", fontface = "italic", size = 5) + coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) 
  res_plots <- list(blank_plot, blank_plot, blank_plot, blank_plot); sol_plots <- list(blank_plot, blank_plot, blank_plot, blank_plot)
  
  design_Ve <- NA 
  design_logL <- NA 
  designplus_logL <- NA
  
  # --- PANEL 1: RAW DATA ---
  res_plots[[1]] <- ggplot(raw_data, aes(x = Ppos, y = Prow)) + geom_tile(aes(fill = .data[[trait]]), color = "lightgray", size = 0.2) + scale_fill_gradientn(colors = ppgmap_colors, na.value = "white", breaks = force_breaks, labels = format_labels) + shared_layers + theme_trial + scale_y_reverse() + coord_fixed() + labs(title = "1. Raw Data")
  sol_plots[[1]] <- res_plots[[1]]
  
  # --- RUN MODELS ---
  for (part in names(models_to_run)) {
    model_name <- models_to_run[[part]]
    cat("  -> Running", model_name, "... ")
    
    # 1. EXECUTE ASREML (SILENTLY)
    system(paste(asreml_exe, "-n", as_file, part, trait, "1.2"), wait = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE) 
    
    # 2. CHECK CONVERGENCE
    sandbox_asr <- paste0(project_name, ".asr"); sandbox_yht <- paste0(project_name, ".yht"); sandbox_sln <- paste0(project_name, ".sln")
    if(!file.exists(sandbox_asr)) { cat("Failed to generate .asr\n"); next }
    
    asr_lines <- readLines(sandbox_asr, warn = FALSE)
    if(!any(grepl("LogL Converged", asr_lines))) { cat("Converge Fail\n"); file.copy(sandbox_asr, paste0(trait, "_", model_name, "_FAILED.asr"), overwrite = TRUE); next }
    cat("Success\n      [Diagnostics]:")
    
    # NEW: DIAGNOSTIC WARNING SCRAPER
    finish_line <- grep("Finished:", asr_lines, value = TRUE)
    if(length(finish_line) > 0 && grepl("Warning|Error", finish_line[1], ignore.case = TRUE)) cat("\n        -> RUN MSG: ", trimws(sub(".*Finished:.*(Warning:|Error:)", "\\1", finish_line[1], ignore.case=TRUE)))
    sing_lines <- grep("singularities detected", asr_lines, ignore.case = TRUE, value = TRUE)
    if(length(sing_lines) > 0) cat("\n        ->", trimws(sing_lines[1]))
    
    mt_start <- grep("Model_Term", asr_lines, ignore.case = TRUE); wald_start <- grep("Wald F statistics", asr_lines, ignore.case = TRUE)
    if(length(mt_start) > 0) {
      bad_vc <- grep("\\s+[B\\?S]\\s*$", asr_lines[mt_start[1]:ifelse(length(wald_start) > 0, wald_start[1] - 1, length(asr_lines))], value = TRUE)
      if(length(bad_vc) > 0) { cat("\n        -> BOUNDARY/SINGULAR COMPONENTS DETECTED:"); for(b in bad_vc) cat("\n             ", trimws(b)) }
    }
    cat("\n")
    
    # 3. SAFELY RENAME FILES
    out_asr <- paste0(trait, "_", model_name, ".asr"); out_yht <- paste0(trait, "_", model_name, ".yht"); out_sln <- paste0(trait, "_", model_name, ".sln")
    file.copy(sandbox_asr, out_asr, overwrite = TRUE); file.copy(sandbox_yht, out_yht, overwrite = TRUE); if(file.exists(sandbox_sln)) file.copy(sandbox_sln, out_sln, overwrite = TRUE)
    
    # Scrape Variances
    logl_val <- as.numeric(gsub("LogL=\\s*", "", str_extract(tail(grep("LogL=", asr_lines, value = TRUE), 1), "LogL=\\s*[-0-9.]+")))    
    
    terms <- c("Block", "SubBlock", "Block\\.SubBlock", "Block\\.Prow", "Block\\.Ppos", "Prow", "Ppos", "Family_id", "Family_name", "uni\\(Crosstype,2\\)", "units", "Residual", "Plot","Tree") 
    current_model_vars <- list() 
    
    # THE FIX: Isolate the actual Variance Components table!
    model_term_idx <- grep("Model_Term", asr_lines, ignore.case = TRUE)
    asr_lines_vars <- if (length(model_term_idx) > 0) asr_lines[model_term_idx[1]:length(asr_lines)] else asr_lines
    
    for(t in terms) {
      lines <- grep(paste0("^\\s*", t), asr_lines_vars, value = TRUE, ignore.case = TRUE)
      lines <- lines[!grepl("effects", lines, ignore.case = TRUE)]
      if(length(lines) > 0) {
        term_clean <- gsub("\\\\", "", t)
        if (tolower(term_clean) == "units") term_clean <- "Independent Error" else if (tolower(term_clean) == "residual") term_clean <- ifelse(part == "3", "Spatial Variance", "Independent Error")
        df_row <- data.frame(Trait=trait, Model=model_name, Term=term_clean, Variance=as.numeric(unlist(strsplit(trimws(lines[1]), "\\s+"))[5]))
        if (!(term_clean == "Independent Error" && any(sapply(current_model_vars, function(x) x$Term == "Independent Error")))) {
          trait_variance_list[[length(trait_variance_list) + 1]] <- df_row; master_results_list[[length(master_results_list) + 1]] <- df_row %>% mutate(LogL=logl_val); current_model_vars[[length(current_model_vars) + 1]] <- df_row
        }
      }
    }
    current_vars_df <- bind_rows(current_model_vars)
    
    # Assign Delta LogL based on the correct baseline
    if (part == "1") { 
      design_logL <- logl_val; d_logL <- "Base"
      ve_row <- current_vars_df %>% filter(Term %in% c("Independent Error", "Spatial Variance"))
      # Bulletproof fallback: sum the variances instead of defaulting to 1
      design_Ve <- if(nrow(ve_row) > 0) ve_row %>% pull(Variance) %>% .[1] else sum(current_vars_df$Variance, na.rm=TRUE)
      
      design_Ve_current <- design_Ve
      
    } else {
      # Fallback: If Model 1 failed, don't leak previous traits. Use current model's Ve.
      if(is.na(design_Ve)) {
        ve_fallback <- current_vars_df %>% filter(Term %in% c("Independent Error", "Spatial Variance"))
        design_Ve_current <- if(nrow(ve_fallback) > 0) ve_fallback %>% pull(Variance) %>% .[1] else sum(current_vars_df$Variance, na.rm=TRUE)
      } else {
        design_Ve_current <- design_Ve
      }
      
      if (part == "2") { 
        designplus_logL <- logl_val
        d_logL <- if(!is.na(design_logL)) round(logl_val - design_logL, 2) else "NA (Base Failed)"
      } else if (part == "3") { 
        d_logL <- if(!is.na(designplus_logL)) round(logl_val - designplus_logL, 2) else "NA (Design+ Failed)" 
      }
    }
    
    # Calculate the Total Phenotypic Variance (Vp) for the current model
    current_Vp <- sum(current_vars_df$Variance, na.rm = TRUE)
    
    # Calculate percentages relative to the Total Variance (Vp), matching Section G
    var_string <- str_wrap(paste(current_vars_df %>% 
                                   mutate(Text = paste0(Term, ": ", round((Variance / current_Vp) * 100, 1), "%")) %>% 
                                   pull(Text), collapse = " | "), width = 75)
    
    metrics_subtitle <- paste0("dLogL: ", d_logL, " | Raw Ve: ", round(current_vars_df %>% filter(Term == "Independent Error") %>% pull(Variance) %>% .[1], 3), "\n", var_string)
    
    if (part == "3") {
      ar_lines <- grep("(?i)AR_[RC]|AR1", asr_lines, value = TRUE)
      metrics_subtitle <- paste0(metrics_subtitle, "\nAR-Row: ", if(length(ar_lines) >= 1) round(as.numeric(unlist(strsplit(trimws(ar_lines[1]), "\\s+"))[5]), 2) else "NA", "  |  AR-Col: ", if(length(ar_lines) >= 2) round(as.numeric(unlist(strsplit(trimws(ar_lines[2]), "\\s+"))[5]), 2) else "NA")
    }
    
    # B. Process Maps & Fixed Effects
    sln <- read.table(out_sln, skip = 1, fill = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric"))
    colnames(sln) <- c("Term", "Level", "Estimate", "SE")
    
    # EXTRACT ALL FIXED EFFECTS (ORIGIN & SPATIAL COVARIATES) & WALD P-VALUE
    for (f_term in c("Origin", "Edge", "Distright", "Distleft", "Distedge")) {
      term_sln <- sln %>% filter(tolower(Term) == tolower(f_term))
      if (nrow(term_sln) > 0) {
        term_p <- "NA"; wald_idx <- grep("Wald F statistics", asr_lines, ignore.case = TRUE)
        if (length(wald_idx) > 0) { term_line <- grep(paste0("\\b", f_term, "\\b"), asr_lines[wald_idx[1]:length(asr_lines)], ignore.case = TRUE, value = TRUE); if (length(term_line) > 0) term_p <- tail(unlist(strsplit(trimws(term_line[1]), "\\s+")), 1) }
        master_fixed_list[[length(master_fixed_list) + 1]] <- term_sln %>% mutate(Trait=trait, Model=model_name, T_value=Estimate/SE, Wald_P_Value=term_p)
      }
    }
    
    yht <- read.table(out_yht, skip = 1); colnames(yht) <- c("Record", "Yhat", "Residual", "Hat")
    map_data <- raw_data %>% left_join(yht, by = "Record") %>% mutate(Resid = Residual)
    
    # Design and Design+ Map Construction
    if (part == "1" || part == "2") {
      b_eff <- sln %>% filter(tolower(Term) == "block") %>% select(Level, Estimate)
      
      # Grab Row and Col effects for Design+ 
      r_eff <- sln %>% filter(tolower(Term) == "block.prow") %>% select(Level, Estimate)
      c_eff <- sln %>% filter(tolower(Term) == "block.ppos") %>% select(Level, Estimate)
      
      # Use a robust regex to catch "Subblock" or "SubBlock" in the raw data
      sub_col_name <- grep("subblock", colnames(raw_data), ignore.case = TRUE, value = TRUE)
      plot_col_name <- grep("^plot$", colnames(raw_data), ignore.case = TRUE, value = TRUE)
      sub_eff <- sln %>% filter(tolower(Term) %in% c("block.subblock", "subblock")) %>% select(Level, Estimate)
      p_eff <- sln %>% filter(tolower(Term) %in% c("plot", "block.plot")) %>% select(Level, Estimate)
      
      map_data <- map_data %>% 
        mutate(Block_key = as.character(Block),
               Sub_key = if(length(sub_col_name) > 0) paste(Block, .data[[sub_col_name[1]]], sep=".") else NA,
               Plot_key = if(length(plot_col_name) > 0) as.character(.data[[plot_col_name[1]]]) else NA,
               BProw_key = sprintf("%d.%03d", suppressWarnings(as.integer(Block)), suppressWarnings(as.integer(Prow))),
               BPpos_key = sprintf("%d.%03d", suppressWarnings(as.integer(Block)), suppressWarnings(as.integer(Ppos)))) %>%
        left_join(b_eff, by = c("Block_key" = "Level")) %>% rename(B_est = Estimate) %>%
        left_join(sub_eff, by = c("Sub_key" = "Level")) %>% rename(Sub_est = Estimate) %>%
        left_join(p_eff, by = c("Plot_key" = "Level")) %>% rename(P_est = Estimate) %>%
        left_join(r_eff, by = c("BProw_key" = "Level")) %>% rename(Row_Est = Estimate) %>%
        left_join(c_eff, by = c("BPpos_key" = "Level")) %>% rename(Col_Est = Estimate) %>%
        mutate(D_Sum = replace_na(B_est,0) + replace_na(Sub_est,0) + replace_na(P_est,0) + replace_na(Row_Est,0) + replace_na(Col_Est,0))
      
      # DYNAMIC EDGE SCRAPER FOR VISUAL MAPS
      for (e_term in c("Edge", "Distright", "Distleft", "Distedge")) {
        if (any(tolower(sln$Term) == tolower(e_term)) && e_term %in% colnames(map_data)) {
          t_sln <- sln[tolower(sln$Term) == tolower(e_term), ]
          if (nrow(t_sln) > 1) { 
            # Scenario A: It's a Factor
            map_data$D_Sum <- map_data$D_Sum + replace_na((map_data %>% mutate(Join_Key = as.character(.data[[e_term]])) %>% left_join(t_sln %>% select(Level, Estimate), by = c("Join_Key" = "Level")))$Estimate, 0)
            cat(paste("    -> Plotted FACTOR edge effect:", e_term, "\n"))
          } else if (nrow(t_sln) == 1) { 
            # Scenario B: It's a Covariate Slope
            map_data$D_Sum <- map_data$D_Sum + replace_na(suppressWarnings(as.numeric(as.character(map_data[[e_term]]))) * t_sln$Estimate[1], 0) 
            cat(paste("    -> Plotted COVARIATE edge effect:", e_term, "\n"))
          }
        }
      }
      map_data <- map_data %>% mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, D_Sum))
    } else { 
      # Spatial Correction Map Construction
      map_data <- map_data %>% mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, Resid), Resid = ifelse(is.na(.data[[trait]]), NA, yht$Yhat - mean(yht$Yhat, na.rm=TRUE))) 
    }
    
    # Define the individual plots
    sol_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = Ppos, y = Prow)) + geom_tile(aes(fill = PlotVal), color = "lightgray", size = 0.2) + scale_fill_gradientn(colors = ppgmap_colors, na.value = "white", breaks = force_breaks, labels = format_labels) + shared_layers + theme_trial + scale_y_reverse() + coord_fixed() + labs(title = model_name, subtitle = metrics_subtitle)
    res_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = Ppos, y = Prow)) + geom_tile(aes(fill = Resid), color = "lightgray", size = 0.2) + scale_fill_gradientn(colors = ppgmap_colors, na.value = "white", breaks = force_breaks, labels = format_labels) + shared_layers + theme_trial + scale_y_reverse() + coord_fixed() + labs(title = paste(model_name, "Resid"), subtitle = metrics_subtitle)
    
    # CLEAN THE SANDBOX FOR THE NEXT LOOP
    file.remove(sandbox_asr); if(file.exists(sandbox_yht)) file.remove(sandbox_yht); if(file.exists(sandbox_sln)) file.remove(sandbox_sln)
  }
  
  # --- F. CANVAS SAVING ---
  # Dynamic 4-Panel Plot Sizing (Tile-Based Scaling)
  if(length(res_plots) == 4) {
    dyn_width  <- max(((max(raw_data$Ppos, na.rm = TRUE) - min(raw_data$Ppos, na.rm = TRUE) + 1) * 0.15 * 2) + 2.5, 10)
    dyn_height <- max(((max(raw_data$Prow, na.rm = TRUE) - min(raw_data$Prow, na.rm = TRUE) + 1) * 0.15 * 2) + 2.5, 7)
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_RESIDUALS.svg")), (res_plots[[1]] | res_plots[[2]]) / (res_plots[[3]] | res_plots[[4]]) + plot_annotation(title = paste(project_name, "- Trait:", trait), theme = theme(plot.title = element_text(size = 18, face = "bold"))), width = dyn_width, height = dyn_height, dpi = 600)
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_SOLUTIONS.svg")), (sol_plots[[1]] | sol_plots[[2]]) / (sol_plots[[3]] | sol_plots[[4]]) + plot_annotation(title = paste(project_name, "- Trait:", trait), theme = theme(plot.title = element_text(size = 18, face = "bold"))), width = dyn_width, height = dyn_height, dpi = 600)
  }
}

# --- G. COMPILING TABLES & REPORTING OUTPUTS ---
if(length(master_results_list) > 0) {
  df_base <- bind_rows(master_results_list) %>% group_by(Trait, Model, Term) %>% slice_tail(n = 1) %>% ungroup()
  df_wide <- bind_rows(df_base %>% distinct(Trait, Model, LogL) %>% mutate(Term = "LogL", Variance = LogL, Pct = NA) %>% select(Trait, Model, Term, Variance, Pct), df_base %>% group_by(Trait, Model) %>% mutate(Total_Var = sum(Variance, na.rm = TRUE), Pct = Variance / Total_Var) %>% ungroup() %>% select(Trait, Model, Term, Variance, Pct)) %>% pivot_wider(names_from = Model, values_from = c(Variance, Pct), names_sep = "_")
  for(col in c("Variance_Design", "Variance_Design+", "Variance_Spatial AR1", "Pct_Design", "Pct_Design+", "Pct_Spatial AR1")) { if(!col %in% colnames(df_wide)) df_wide[[col]] <- NA }
  
  df_calc <- df_wide %>% mutate(`Delta_Design+` = `Variance_Design+` - Variance_Design, `Delta_Spatial` = `Variance_Spatial AR1` - `Variance_Design+`, `Pct_Delta_Design+` = `Pct_Design+` - Pct_Design, `Pct_Delta_Spatial` = `Pct_Spatial AR1` - `Pct_Design+`) %>% select(Trait, Term, Variance_Design, `Variance_Design+`, `Variance_Spatial AR1`, `Delta_Design+`, `Delta_Spatial`, Pct_Design, `Pct_Design+`, `Pct_Spatial AR1`, `Pct_Delta_Design+`, `Pct_Delta_Spatial`) %>% arrange(Trait, desc(Term == "LogL"), Term)
  
  write.csv(df_base %>% mutate(Trial = project_name) %>% select(Trial, Trait, Model, Term, Variance, LogL), file.path(out_dir, "All_Traits_Variance_Summary_Long.csv"), row.names = FALSE, na = "")
  write.csv(df_calc, file.path(out_dir, "All_Traits_Variance_Summary_Wide_Deltas.csv"), row.names = FALSE, na = "")
  
  # Flextable Generation
  doc <- read_docx(); trait_list <- split(df_calc, df_calc$Trait)
  for (tr in names(trait_list)) {
    
    # 1. Prepare the data by dropping 'Trait' first so column counts are true
    table_data <- trait_list[[tr]] %>% select(-Trait)
    num_cols <- ncol(table_data)
    
    # 2. Build the table using the safe dynamic boundary
    ft <- flextable(table_data %>% 
                      mutate(across(starts_with("Pct"), ~ ifelse(is.na(.), "", sprintf("%.1f%%", . * 100)))) %>% 
                      mutate(across(where(is.numeric), ~ ifelse(is.na(.), "", sprintf("%.3f", .))))) %>% 
      theme_booktabs() %>% 
      fontsize(size = 9, part = "all") %>% 
      padding(padding = 3, part = "all") %>% 
      set_header_labels(Variance_Design = "Var\nDesign", `Variance_Design+` = "Var\nDesign+", `Variance_Spatial AR1` = "Var\nSpatial", `Delta_Design+` = "Delta\nDesign+", Delta_Spatial = "Delta\nSpatial", Pct_Design = "%\nDesign", `Pct_Design+` = "%\nDesign+", `Pct_Spatial AR1` = "%\nSpatial", `Pct_Delta_Design+` = "% Delta\nDesign+", Pct_Delta_Spatial = "% Delta\nSpatial") %>% 
      autofit() %>% 
      align(j = "Term", align = "left", part = "all") %>% 
      align(j = 2:num_cols, align = "center", part = "all") # FIX: Safely targets exactly what's available
    
    doc <- doc %>% body_add_par(paste("Trait:", tr), style = "heading 2") %>% body_add_flextable(value = ft, align = "left") %>% body_add_par("", style = "Normal") 
  }
  
  print(doc %>% body_end_section_landscape(), target = file.path(out_dir, "All_Traits_Variance_Summary_Formatted.docx"))
}

# --- H. EDGE COVARIATE FORMAL TESTING ---
base_path <- file.path(trial_folder, "Analyses", "Baseline", "All_Traits_Variance_Summary_Long.csv")
edge_path <- file.path(trial_folder, "Analyses", "DistEdge", "All_Traits_Variance_Summary_Long.csv")

if(file.exists(base_path) && file.exists(edge_path)) {
  comparison_results <- read_csv(base_path) %>% select(Trait, Model, LogL) %>% distinct() %>% rename(LogL_Baseline = LogL) %>%
    inner_join(read_csv(edge_path) %>% select(Trait, Model, LogL) %>% distinct() %>% rename(LogL_DistEdge = LogL), by = c("Trait", "Model")) %>%
    mutate(Delta_LogL = LogL_DistEdge - LogL_Baseline, LRT_Stat = 2 * Delta_LogL, Significant = LRT_Stat > 3.84)
  write.csv(comparison_results, file.path(trial_folder, "Analyses", "Model_Comparison_Results.csv"), row.names = FALSE)
}
cat("\nPipeline Process Completed Successfully.\n")