# 
# 0. CONTROL PANEL (Change these for your specific project) #### 
# 

trial_folder  <- "C:/Users/james.baker/Forest Research/TW CBC-TBA-NextGenBritishConifers - Share/Sitka/Backwards Selected Fullsib P96-P99 experiments/Kintyre 18"
project_name  <- "Kintyre_18_S"
as_file       <- paste0(project_name, ".as")
csv_file      <- paste0(project_name, ".csv")

# Path to ASReml Standalone (Usually standard across company machines)
asreml_path   <- "C:/Program Files/ASReml4/bin/asreml.exe"

# 
# 1. SETUP & LIBRARIES (Do not edit below this line)
# 
library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)

asreml_exe <- shQuote(asreml_path) 
raw_data <- read.csv(file.path(trial_folder,csv_file), stringsAsFactors = FALSE)

# Lock in the 'data file order' exactly as the PPGVal Sfile step output it.
# This ensures ASReml's un-indexed .yht residuals will always join flawlessly.
raw_data <- raw_data %>%
  mutate(Record = row_number())

# --- DYNAMIC REGEX TRAIT SCRAPER WITH GUARD RAIL ---
as_text <- paste(readLines(file.path(trial_folder, as_file)), collapse = " ")
found_traits <- unique(unlist(str_extract_all(as_text, "\\b[A-Za-z0-9]+_[0-9]+\\b")))
traits_to_test <- found_traits[found_traits %in% colnames(raw_data)]

cat("Automated Discovery: Found", length(found_traits), "potential matches.\n")
cat("Guard Rail: Proceeding with", length(traits_to_test), "traits found in CSV.\n")

out_dir <- file.path(trial_folder, "Analyses")
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
  cat("Created 'Analyses' folder for outputs in the trial directory.\n")
}

# Ensure R looks for the outputs inside the trial folder
base_name <- file.path(trial_folder, project_name) 
out_asr <- paste0(base_name, ".asr")
out_yht <- paste0(base_name, ".yht")
out_sln <- paste0(base_name, ".sln")

models_to_run <- c("1" = "Design", "2" = "Design+", "3" = "Spatial AR1")
master_results_list <- list()

# Calculate Block Boundaries for the Black Outlines & Text Labels
block_bounds <- raw_data %>%
  filter(!is.na(Ppos) & !is.na(Prow) & !is.na(Block)) %>%
  group_by(Block) %>%
  summarize(
    xmin = min(as.numeric(Ppos), na.rm = TRUE) - 0.5,
    xmax = max(as.numeric(Ppos), na.rm = TRUE) + 0.5,
    ymin = min(as.numeric(Prow), na.rm = TRUE) - 0.5,
    ymax = max(as.numeric(Prow), na.rm = TRUE) + 0.5,
    # NEW: Calculate the geometric center for the labels
    x_mid = (xmin + xmax) / 2,
    y_mid = (ymin + ymax) / 2
  )

ppgmap_colors <- c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred")

# 
# 2. Execution section ####
# 
for (trait in traits_to_test) {
  cat("\n========================================\nProcessing:", trait, "\n")
  
  raw_data[[trait]] <- suppressWarnings(as.numeric(as.character(raw_data[[trait]])))
  res_plots <- list(); sol_plots <- list(); trait_variance_list <- list()
  
  # --- PANEL 1: RAW DATA ---
  raw_map <- ggplot(raw_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = .data[[trait]])) +
    geom_tile(color = "black", size = 0.05) + 
    geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
    scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90", n.breaks=10,guide = guide_colorbar(barheight = 15)) +
    geom_label(data = block_bounds, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 4, fontface = "bold", fill = alpha("white", 0.6), label.size = NA) +
    scale_y_reverse() + theme_void() + coord_fixed() + labs(title = paste("1. Raw Data:", trait))
  res_plots[[1]] <- raw_map; sol_plots[[1]] <- raw_map
  
  # --- RUN MODELS ---
  for (part in names(models_to_run)) {
    model_name <- models_to_run[[part]]
    cat("  -> Running", model_name, "... ")
    
    suppressWarnings(file.remove(out_asr, out_yht, out_sln))
    # Point list.files specifically to the trial folder to clean up junk
    suppressWarnings(file.remove(list.files(path = trial_folder, pattern = paste0("^", project_name, "\\.(msv|veo|ask|tmp|tsv)$"), full.names = TRUE)))
    
    # 1. Save our current location, then step into the Shared Drive
    original_wd <- getwd()
    setwd(trial_folder)
   
    # 2. Build the command WITH the silent flag (> NUL 2>&1)
    command <- paste(asreml_exe, "-n", paste0(project_name, ".as"), part, trait, "1.2 > NUL 2>&1")
    
    # Print the command to the console for tracking
    cat("  -> Command:", command, "\n")
    
    # 3. Run ASReml silently
    shell(command, wait = TRUE)
    # 3. Step back to our local R project directory
    setwd(original_wd)
    
    if(!file.exists(out_asr)) { cat("FAILED\n"); next }
    asr_lines <- readLines(out_asr)
    if(!any(grepl("LogL Converged", asr_lines))) { cat("Converge Fail\n"); next }
    cat("Success\n")
    
    # A. Scrape Variances
    logl_val <- as.numeric(str_extract(tail(grep("LogL=", asr_lines, value = TRUE), 1), "(?<=LogL=\\s)[-0-9.]+"))
    
    # UPGRADE: Added regex to catch AR1/COR spatial parameters
    terms <- c("Block", "SubBlock", "Block\\.SubBlock", "Block\\.Prow", "Block\\.Ppos", "Family_id", "uni\\(Crosstype,2\\)", "units", "Tree", "Prow.*ar1", "Ppos.*ar1", "Prow.*cor", "Ppos.*cor")
    
    current_model_vars <- list() # Store just this model's terms temporarily
    
    for(t in terms) {
      line <- grep(paste0("^\\s*", t, "\\s+"), asr_lines, value = TRUE)
      if(length(line) > 0) {
        val <- as.numeric(unlist(strsplit(trimws(line[1]), "\\s+"))[5])
        df_row <- data.frame(Trait=trait, Model=model_name, Term=gsub("\\\\", "", t), Variance=val)
        trait_variance_list[[length(trait_variance_list) + 1]] <- df_row
        master_results_list[[length(master_results_list) + 1]] <- df_row %>% mutate(LogL=logl_val)
        current_model_vars[[length(current_model_vars) + 1]] <- df_row
      }
    }
    
    current_vars_df <- bind_rows(current_model_vars)
    
    # CAPTURE THE BASELINES FROM THE DESIGN MODEL (Part 1)
    if (part == "1") {
      design_logL <- logl_val
      if ("units" %in% current_vars_df$Term) {
        design_Ve <- current_vars_df %>% filter(Term == "units") %>% pull(Variance) %>% .[1]
      } else {
        design_Ve <- 1 # Failsafe to prevent division by zero
      }
    }
    
    # BUILD THE METRICS SUBTITLE FOR THE MAPS
    d_logL <- round(logl_val - design_logL, 2)
    curr_Ve <- current_vars_df %>% filter(Term == "units") %>% pull(Variance) %>% round(3) %>% .[1]
    
    metrics_subtitle <- paste0("dLogL: ", d_logL, " | Ve: ", curr_Ve)
    
    # Add spatial auto-correlations to the subtitle if it's the spatial model
    if (part == "3") {
      ar_row <- current_vars_df %>% filter(grepl("Prow.*(ar1|cor)", Term)) %>% pull(Variance) %>% round(3) %>% .[1]
      ar_col <- current_vars_df %>% filter(grepl("Ppos.*(ar1|cor)", Term)) %>% pull(Variance) %>% round(3) %>% .[1]
      metrics_subtitle <- paste0(metrics_subtitle, " | AR-Row: ", ifelse(is.na(ar_row), "NA", ar_row), " | AR-Col: ", ifelse(is.na(ar_col), "NA", ar_col))
    }
    
    # B. Process Maps (CRITICAL FIX: colClasses forces 'Level' to remain text)
    sln <- read.table(out_sln, skip = 1, fill = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric"))
    colnames(sln) <- c("Term", "Level", "Estimate", "SE")
    
    yht <- read.table(out_yht, skip = 1)
    colnames(yht) <- c("Record", "Yhat", "Residual", "Hat")
    
    # JOIN perfectly based on the 'data file order' Record index!
    map_data <- raw_data %>% 
      left_join(yht, by = "Record") %>% 
      mutate(Resid = ifelse(is.na(.data[[trait]]), NA, Residual))
    
    if (part == "1") {
      block_eff <- sln %>% filter(Term == "Block") %>% select(Level, Estimate)
      map_data <- map_data %>% 
        mutate(Block_key = as.character(Block)) %>%
        left_join(block_eff, by = c("Block_key" = "Level")) %>%
        mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, Estimate))
      sol_title <- paste("2. Design (Block Solutions):",trait)
      
    } else if (part == "2") {
      b_eff <- sln %>% filter(Term == "Block") %>% select(Level, Estimate)
      r_eff <- sln %>% filter(Term == "Block.Prow") %>% select(Level, Estimate)
      p_eff <- sln %>% filter(Term == "Block.Ppos") %>% select(Level, Estimate)
      
      # FIX: Build exact matching strings (e.g. 1.001) using sprintf
      map_data <- map_data %>% 
        mutate(
          Block_key = as.character(Block),
          BProw_key = sprintf("%d.%03d", as.integer(Block), as.integer(Prow)),
          BPpos_key = sprintf("%d.%03d", as.integer(Block), as.integer(Ppos))
        ) %>%
        left_join(b_eff, by = c("Block_key" = "Level")) %>% rename(B_est = Estimate) %>%
        left_join(r_eff, by = c("BProw_key" = "Level")) %>% rename(R_est = Estimate) %>%
        left_join(p_eff, by = c("BPpos_key" = "Level")) %>% rename(P_est = Estimate) %>%
        mutate(D_Sum = replace_na(B_est,0) + replace_na(R_est,0) + replace_na(P_est,0)) %>%
        mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, D_Sum))
      sol_title <- paste("3. Design+ (Blk+Row+Col Solutions):", trait)
      
    } else {
      map_data <- map_data %>% 
        mutate(
          # Map the smooth spatial surface to the Solutions plot
          PlotVal = ifelse(is.na(.data[[trait]]), NA, Resid), 
          
          # Map the independent genetic/design noise to the Residuals plot 
          # (Subtracting the mean centers it around 0 for the color scale)
          Resid = ifelse(is.na(.data[[trait]]), NA, yht$Yhat - mean(yht$Yhat, na.rm=TRUE)) 
        )
      sol_title <- paste("4. Spatial Effects Correction:", trait)
    }
    
    sol_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = PlotVal)) +
      geom_tile(color = "black", size = 0.05) + 
      geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90", n.breaks=10, guide = guide_colorbar(barheight = 15)) +
      geom_label(data = block_bounds, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 4, fontface = "bold", fill = alpha("white", 0.6), label.size = NA) +
      scale_y_reverse() + theme_void() + coord_fixed() + labs(title = sol_title, subtitle = metrics_subtitle)
    
    res_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = Resid)) +
      geom_tile(color = "black", size = 0.05) + 
      geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90", n.breaks=10,guide = guide_colorbar(barheight = 15)) +
      geom_label(data = block_bounds, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 4, fontface = "bold", fill = alpha("white", 0.6), label.size = NA) +
      scale_y_reverse() + theme_void() + coord_fixed() + labs(title = paste0(as.numeric(part)+1, ". ", model_name, " Res: ",trait), subtitle = metrics_subtitle)    
    
    # --- ARCHIVE ALL FILES ---
    file.copy(out_asr, file.path(out_dir, paste0(trait, "_", model_name, ".asr")), overwrite = TRUE)
    if(file.exists(out_sln)) file.copy(out_sln, file.path(out_dir, paste0(trait, "_", model_name, ".sln")), overwrite = TRUE)
    if(file.exists(out_yht)) file.copy(out_yht, file.path(out_dir, paste0(trait, "_", model_name, ".yht")), overwrite = TRUE)
  }
  
  if(length(res_plots) == 4) {
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_RESIDUALS.png")), (res_plots[[1]] | res_plots[[2]]) / (res_plots[[3]] | res_plots[[4]]), width = 12, height = 10)
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_SOLUTIONS.png")), (sol_plots[[1]] | sol_plots[[2]]) / (sol_plots[[3]] | sol_plots[[4]]), width = 12, height = 10)
  }
  
  if(length(trait_variance_list) > 0) {
    trait_df <- bind_rows(trait_variance_list) %>% 
      # UPGRADE: Scale all variances relative to the Design Model's Ve
      mutate(Scaled_Var = Variance / design_Ve) %>% 
      ungroup()
    
    trait_df$Model <- factor(trait_df$Model, levels = c("Design", "Design+", "Spatial AR1"))
    
    bp <- ggplot(trait_df, aes(x = Model, y = Scaled_Var, fill = Term)) +
      geom_col(color = "black") + theme_minimal() + scale_fill_brewer(palette = "Set3") +
      labs(title = paste("Variance Components:", trait), 
           subtitle = paste("Scaled relative to Design Ve (", round(design_Ve, 3), ")"),
           y = "Variance (Proportion of Design Ve)") 
    ggsave(file.path(out_dir, paste0(trait, "_VC_Barplot.png")), bp, width = 7, height = 6)
  }
}

# 3. EXPORT FINAL MASTER TABLE (Wide-by-Model with Deltas)

if(length(master_results_list) > 0) {
  
  # 1. Clean duplicates and base data
  df_base <- bind_rows(master_results_list) %>% 
    group_by(Trait, Model, Term) %>% slice_tail(n = 1) %>% ungroup()
  
  # 2. Calculate Percentages (Variance / Total Variance per Model)
  df_var <- df_base %>%
    group_by(Trait, Model) %>%
    mutate(
      Total_Var = sum(Variance, na.rm = TRUE),
      Pct = Variance / Total_Var # Keep as decimal (e.g. 0.31 for 31%) so Excel reads it well
    ) %>% 
    ungroup()
  
  # 3. Format LogL as a "Term" row so it stacks nicely
  df_logl <- df_base %>%
    distinct(Trait, Model, LogL) %>%
    mutate(Term = "LogL", Variance = LogL, Pct = NA) %>%
    select(Trait, Model, Term, Variance, Pct)
  
  # 4. Combine and Pivot Wide by Model
  df_combined <- bind_rows(df_logl, df_var %>% select(Trait, Model, Term, Variance, Pct))
  
  # We use backticks to handle the `+` and spaces in model names safely
  df_wide <- df_combined %>%
    pivot_wider(
      names_from = Model, 
      values_from = c(Variance, Pct),
      names_sep = "_"
    )
  
  # 5. Safety Check: Ensure columns exist even if a model failed to converge
  expected_cols <- c("Variance_Design", "Variance_Design+", "Variance_Spatial AR1",
                     "Pct_Design", "Pct_Design+", "Pct_Spatial AR1")
  for(col in expected_cols) {
    if(!col %in% colnames(df_wide)) df_wide[[col]] <- NA
  }
  
  # 6. Calculate Deltas
  df_calc <- df_wide %>%
    mutate(
      `Delta_Design+` = `Variance_Design+` - Variance_Design,
      `Delta_Spatial` = `Variance_Spatial AR1` - `Variance_Design+`,
      `Pct_Delta_Design+` = `Pct_Design+` - Pct_Design,
      `Pct_Delta_Spatial` = `Pct_Spatial AR1` - `Pct_Design+`
    ) %>%
    # Order the columns exactly like the Excel image
    select(
      Trait, Term,
      Variance_Design, `Variance_Design+`, `Variance_Spatial AR1`,
      `Delta_Design+`, `Delta_Spatial`,
      Pct_Design, `Pct_Design+`, `Pct_Spatial AR1`,
      `Pct_Delta_Design+`, `Pct_Delta_Spatial`
    ) %>%
    # Sort so LogL is at the top of each trait block
    arrange(Trait, desc(Term == "LogL"), Term)
  
  # 7. VISUAL SEPARATION: Insert a blank header row between traits
  trait_list <- split(df_calc, df_calc$Trait)
  spaced_list <- lapply(names(trait_list), function(tr) {
    # Create an empty row
    blank_row <- df_calc[1, ]
    blank_row[1, ] <- NA
    # Turn the Trait column into a visual header
    blank_row$Trait <- paste(">>> TRAIT:", tr, "<<<")
    
    bind_rows(blank_row, trait_list[[tr]])
  })
  
  final_table <- bind_rows(spaced_list)
  
  # 8. Export!
  write.csv(final_table, file.path(out_dir, "All_Traits_Variance_Summary.csv"), row.names = FALSE, na = "")
  cat("\nSUCCESS! Formatted CSV with Deltas and Percentages saved to 'Analyses'.\n")
} else {
  cat("\n[!] WARNING: No results were captured to save.\n")
}
