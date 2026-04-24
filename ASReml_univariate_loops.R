# 
# 0. CONTROL PANEL (Change these for your specific project) #### 
# 

trial_folder  <- "C:/Users/james.baker/Forest Research/TW CBC-TBA-NextGenBritishConifers - Share/Sitka/Backwards Selected Fullsib P96-P99 experiments/Kintyre 18"
project_name  <- "Kintyre_18_S"
as_file       <- paste0(project_name, ".as")
csv_file      <- paste0(project_name, ".csv")

# NEW: Change this when running different model testing scenarios!
run_id        <- "Baseline"

# Path to ASReml Standalone (Usually standard across company machines)
asreml_path   <- "C:/Program Files/ASReml4/bin/asreml.exe"

# 
# 1. SETUP & LIBRARIES (Do not edit below this line)
# 
library(here)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(svglite)
library(flextable)
library(officer)

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

traits_to_test <-c("Dre_28","Drs_28")

# --- 1. SETUP THE SANDBOX DIRECTORY ---
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
out_dir <- sandbox_dir # Keep out_dir variable so your plot/word exports still work

models_to_run <- c("1" = "Design", "2" = "Design+", "3" = "Spatial AR1")

master_results_list <- list()
master_fixed_list <- list() # NEW: To hold Origin and other fixed effects

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
  
  # NEW: Create a blank failsafe plot that holds its rigid shape!
  blank_plot <- ggplot() + 
    theme_void() + 
    annotate("text", x = 0.5, y = 0.5, label = "Model Did Not Converge", color = "darkred", fontface = "italic", size = 5) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) 
  
  # Pre-fill the lists with the blanks so the 4-panel grid never collapses
  res_plots <- list(blank_plot, blank_plot, blank_plot, blank_plot)
  sol_plots <- list(blank_plot, blank_plot, blank_plot, blank_plot)
  trait_variance_list <- list()
  
  design_Ve <- NA 
  design_logL <- NA
  
  # --- PANEL 1: RAW DATA ---
  force_breaks <- function(x) { 
    min_val <- min(x, na.rm = TRUE)
    max_val <- max(x, na.rm = TRUE)
    if (!is.finite(min_val) || !is.finite(max_val)) return(c(0, 1)) # Failsafe for all NAs
    seq(min_val, max_val, length.out = 5) 
  }
  
  format_labels <- function(x) sprintf("%.2f", x)
  
  raw_map <- ggplot(raw_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = .data[[trait]])) +
    geom_tile(color = "black", size = 0.05) + 
    geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
    scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90", breaks = force_breaks, labels = format_labels) +
    theme_void() + 
    theme(legend.key.height = unit(1.5, "cm")) +
    geom_label(data = block_bounds, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 4, fontface = "bold", fill = alpha("white", 0.6), label.size = NA) +
    scale_y_reverse() + coord_fixed() + labs(title = "1. Raw Data")
  res_plots[[1]] <- raw_map; sol_plots[[1]] <- raw_map
  
  # --- RUN MODELS ---
  for (part in names(models_to_run)) {
    model_name <- models_to_run[[part]]
    cat("  -> Running", model_name, "... ")
    
    # --- 1. EXECUTE ASREML (SILENTLY) ---
    command <- paste(asreml_exe, "-n", as_file, part, trait, "1.2")
    system(command, wait = TRUE, ignore.stdout = TRUE, ignore.stderr = TRUE) 
    
    # --- 2. CHECK CONVERGENCE ---
    sandbox_asr <- paste0(project_name, ".asr")
    sandbox_yht <- paste0(project_name, ".yht")
    sandbox_sln <- paste0(project_name, ".sln")
    
    if(!file.exists(sandbox_asr)) { cat("Failed to generate .asr\n"); next }
    
    asr_lines <- readLines(sandbox_asr, warn = FALSE)
    
    if(!any(grepl("LogL Converged", asr_lines))) { 
      cat("Converge Fail\n")
      file.copy(sandbox_asr, paste0(trait, "_", model_name, "_FAILED.asr"), overwrite = TRUE)
      next 
    }
    cat("Success\n")
    
    # --- 3. SAFELY RENAME FILES ---
    out_asr <- paste0(trait, "_", model_name, ".asr")
    out_yht <- paste0(trait, "_", model_name, ".yht")
    out_sln <- paste0(trait, "_", model_name, ".sln")
    
    file.copy(sandbox_asr, out_asr, overwrite = TRUE)
    file.copy(sandbox_yht, out_yht, overwrite = TRUE)
    if(file.exists(sandbox_sln)) file.copy(sandbox_sln, out_sln, overwrite = TRUE)
    
    # A. Scrape Variances
    logl_line <- tail(grep("LogL=", asr_lines, value = TRUE), 1)
    logl_match <- str_extract(logl_line, "LogL=\\s*[-0-9.]+")
    logl_val <- as.numeric(gsub("LogL=\\s*", "", logl_match))    
    
    terms <- c("Block", "SubBlock", "Block\\.SubBlock", "Block\\.Prow", "Block\\.Ppos", "Prow", "Ppos", "Family_id", "Family_name", "uni\\(Crosstype,2\\)", "units", "Residual", "Tree") 
    
    current_model_vars <- list() 
    
    # --- THE FIX: Isolate the actual Variance Components table! ---
    model_term_idx <- grep("Model_Term", asr_lines, ignore.case = TRUE)
    asr_lines_vars <- if (length(model_term_idx) > 0) asr_lines[model_term_idx[1]:length(asr_lines)] else asr_lines
    
    for(t in terms) {
      # Search only the isolated bottom table
      lines <- grep(paste0("^\\s*", t), asr_lines_vars, value = TRUE, ignore.case = TRUE)
      lines <- lines[!grepl("effects", lines, ignore.case = TRUE)]
      
      if(length(lines) > 0) {
        parts <- unlist(strsplit(trimws(lines[1]), "\\s+"))
        val <- as.numeric(parts[5]) 
        term_clean <- gsub("\\\\", "", t)
        
        if (tolower(term_clean) == "units") {
          term_clean <- "Independent Error"
        } else if (tolower(term_clean) == "residual") {
          if (part == "3") {
            term_clean <- "Spatial Variance"
          } else {
            term_clean <- "Independent Error"
          }
        }
        
        df_row <- data.frame(Trait=trait, Model=model_name, Term=term_clean, Variance=val)
        
        if (!(term_clean == "Independent Error" && any(sapply(current_model_vars, function(x) x$Term == "Independent Error")))) {
          trait_variance_list[[length(trait_variance_list) + 1]] <- df_row
          master_results_list[[length(master_results_list) + 1]] <- df_row %>% mutate(LogL=logl_val)
          current_model_vars[[length(current_model_vars) + 1]] <- df_row
        }
      }
    }
    
    current_vars_df <- bind_rows(current_model_vars)
    
    if (part == "1") {
      design_logL <- logl_val
      ve_row <- current_vars_df %>% filter(Term == "Independent Error")
      design_Ve <- if(nrow(ve_row) > 0) ve_row %>% pull(Variance) %>% .[1] else 1
    }
    
    d_logL <- round(logl_val - design_logL, 2)
    curr_Ve_raw <- current_vars_df %>% filter(Term == "Independent Error") %>% pull(Variance) %>% .[1]
    
    var_text_list <- current_vars_df %>%
      mutate(Pct = round((Variance / design_Ve) * 100, 1)) %>%
      mutate(Text = paste0(Term, ": ", Pct, "%")) %>%
      pull(Text)
    
    var_string <- str_wrap(paste(var_text_list, collapse = " | "), width = 75)
    metrics_subtitle <- paste0("dLogL: ", d_logL, " | Raw Ve: ", round(curr_Ve_raw, 3), "\n", var_string)
    
    if (part == "3") {
      ar_lines <- grep("(?i)AR_[RC]|AR1", asr_lines, value = TRUE)
      
      ar_row <- if(length(ar_lines) >= 1) {
        parts <- unlist(strsplit(trimws(ar_lines[1]), "\\s+"))
        round(as.numeric(parts[5]), 2) 
      } else "NA"
      
      ar_col <- if(length(ar_lines) >= 2) {
        parts <- unlist(strsplit(trimws(ar_lines[2]), "\\s+"))
        round(as.numeric(parts[5]), 2) 
      } else "NA"
      
      metrics_subtitle <- paste0(metrics_subtitle, "\nAR-Row: ", ar_row, "  |  AR-Col: ", ar_col)
    }
    
    # B. Process Maps & Fixed Effects
    sln <- read.table(out_sln, skip = 1, fill = TRUE, stringsAsFactors = FALSE, colClasses = c("character", "character", "numeric", "numeric"))
    colnames(sln) <- c("Term", "Level", "Estimate", "SE")
    
    
    # --- NEW: EXTRACT ALL FIXED EFFECTS (ORIGIN & SPATIAL COVARIATES) & WALD P-VALUE ---
    fixed_terms_to_scrape <- c("Origin", "Edge", "Distright", "Distleft", "Distedge")
    
    for (f_term in fixed_terms_to_scrape) {
      # Use case-insensitive matching to protect against ASReml's weird capitalization
      term_sln <- sln %>% filter(tolower(Term) == tolower(f_term))
      
      if (nrow(term_sln) > 0) {
        term_p <- "NA"
        
        # 1. Find exactly where the Wald table starts
        wald_idx <- grep("Wald F statistics", asr_lines, ignore.case = TRUE)
        
        if (length(wald_idx) > 0) {
          # 2. Only search the text from the Wald table downwards
          wald_section <- asr_lines[wald_idx[1]:length(asr_lines)]
          
          # 3. Find the specific line for this term
          term_line <- grep(paste0("\\b", f_term, "\\b"), wald_section, ignore.case = TRUE, value = TRUE)
          
          if (length(term_line) > 0) {
            # 4. Split by spaces and grab the very last string (the p-value)
            parts <- unlist(strsplit(trimws(term_line[1]), "\\s+"))
            term_p <- tail(parts, 1)
          }
        }
        
        term_df <- term_sln %>%
          mutate(
            Trait = trait,
            Model = model_name,
            T_value = Estimate / SE,
            Wald_P_Value = term_p
          )
        
        master_fixed_list[[length(master_fixed_list) + 1]] <- term_df
      }
    }
    
    yht <- read.table(out_yht, skip = 1)
    colnames(yht) <- c("Record", "Yhat", "Residual", "Hat")
    
    map_data <- raw_data %>% 
      left_join(yht, by = "Record") %>% 
      mutate(Resid = ifelse(is.na(.data[[trait]]), NA, Residual))
    
    if (part == "1") {
      block_eff <- sln %>% filter(Term == "Block") %>% select(Level, Estimate)
      map_data <- map_data %>% 
        mutate(Block_key = as.character(Block)) %>%
        left_join(block_eff, by = c("Block_key" = "Level")) %>%
        mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, Estimate))
      sol_title <- "2. Design (Block Solutions)"
      
    } else if (part == "2") {
      b_eff <- sln %>% filter(Term == "Block") %>% select(Level, Estimate)
      r_eff <- sln %>% filter(Term == "Block.Prow") %>% select(Level, Estimate)
      p_eff <- sln %>% filter(Term == "Block.Ppos") %>% select(Level, Estimate)
      
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
      sol_title <- "3. Design+ (Blk+Row+Col Solutions)"
      
    } else {
      map_data <- map_data %>% 
        mutate(
          PlotVal = ifelse(is.na(.data[[trait]]), NA, Resid), 
          Resid = ifelse(is.na(.data[[trait]]), NA, yht$Yhat - mean(yht$Yhat, na.rm=TRUE)) 
        )
      sol_title <- "4. Spatial Effects Correction"
    }
    
    sol_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = PlotVal)) +
      geom_tile(color = "black", size = 0.05) + 
      geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90", breaks = force_breaks, labels = format_labels) +
      theme_void() + 
      theme(legend.key.height = unit(1.5, "cm"), plot.subtitle = element_text(size = 8)) +
      geom_label(data = block_bounds, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 4, fontface = "bold", fill = alpha("white", 0.6), label.size = NA) +
      scale_y_reverse() + coord_fixed() + labs(title = sol_title, subtitle = metrics_subtitle)
    
    res_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = Resid)) +
      geom_tile(color = "black", size = 0.05) + 
      geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90", breaks = force_breaks, labels = format_labels) +
      theme_void() + 
      theme(legend.key.height = unit(1.5, "cm"), plot.subtitle = element_text(size = 8)) +
      geom_label(data = block_bounds, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 4, fontface = "bold", fill = alpha("white", 0.6), label.size = NA) +
      scale_y_reverse() + coord_fixed() + labs(title = paste0(as.numeric(part)+1, ". ", model_name, " Res"), subtitle = metrics_subtitle)
    
    # --- 4. CLEAN THE SANDBOX FOR THE NEXT LOOP ---
    file.remove(sandbox_asr)
    if(file.exists(sandbox_yht)) file.remove(sandbox_yht)
    if(file.exists(sandbox_sln)) file.remove(sandbox_sln)
  }
  
  if(length(res_plots) == 4) {
    master_title <- paste(project_name, "- Trait:", trait)
    
    res_assembled <- (res_plots[[1]] | res_plots[[2]]) / (res_plots[[3]] | res_plots[[4]]) +
      plot_annotation(title = master_title, theme = theme(plot.title = element_text(size = 18, face = "bold")))
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_RESIDUALS.svg")), res_assembled, width = 12, height = 10,dpi = 600)
    
    sol_assembled <- (sol_plots[[1]] | sol_plots[[2]]) / (sol_plots[[3]] | sol_plots[[4]]) +
      plot_annotation(title = master_title, theme = theme(plot.title = element_text(size = 18, face = "bold")))
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_SOLUTIONS.svg")), sol_assembled, width = 12, height = 10,dpi = 600)
  }
  
  if(length(trait_variance_list) > 0) {
    trait_df <- bind_rows(trait_variance_list) %>% 
      mutate(Pct_Var = (Variance / design_Ve) * 100) %>% 
      ungroup()
    
    trait_df$Model <- factor(trait_df$Model, levels = c("Design", "Design+", "Spatial AR1"))
    
    all_terms <- unique(trait_df$Term)
    other_terms <- setdiff(all_terms, c("Independent Error", "Spatial Variance"))
    
    ordered_levels <- c(sort(other_terms), "Spatial Variance", "Independent Error")
    trait_df$Term <- factor(trait_df$Term, levels = ordered_levels)
    
    bp <- ggplot(trait_df, aes(x = Model, y = Pct_Var, fill = Term)) +
      geom_col(color = "black") + 
      geom_hline(yintercept = 100, linetype = "dashed", color = "red", linewidth = 1) +
      theme_minimal() + 
      scale_fill_brewer(palette = "Set3") +
      labs(title = paste("Variance Components:", trait), 
           y = paste0("% of Design Ve (", round(design_Ve, 3), ")")) +
      theme(legend.position = "right")
    
    ggsave(file.path(out_dir, paste0(trait, "_VC_Barplot.svg")), bp, width = 7, height = 6,dpi = 600)
  }
  
  if(length(master_fixed_list) > 0) {
    origin_mapping <- c(
      "1" = "Ro_North_HG", "2" = "Ro_North_Unk", "3" = "Ro_North_WC",
      "4" = "Ro_South_Unk", "5" = "Ro_HG", "6" = "Ro_Ledmore",
      "7" = "Ro_WC", "8" = "Ro_Filler"
    )
    
    # NEW: Plot ONLY Origin so Edge doesn't crash the graph
    trait_origin_df <- bind_rows(master_fixed_list) %>% 
      filter(Trait == trait & Model == "Spatial AR1" & Term == "Origin") %>%
      mutate(Origin_Name = factor(origin_mapping[as.character(Level)], levels = origin_mapping))
    
    if(nrow(trait_origin_df) > 0) {
      origin_p <- trait_origin_df$Wald_P_Value[1]
      
      op_plot <- ggplot(trait_origin_df, aes(x = Origin_Name, y = Estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2, color = "darkblue") +
        geom_point(size = 3, color = "darkblue") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        theme_minimal() +
        labs(
          title = paste("Origin Solutions (Spatial AR1):", trait),
          subtitle = paste("Wald P-value:", origin_p),
          x = "Origin",
          y = "Solution Estimate (+/- 1 SE)"
        ) +
        scale_x_discrete(drop = TRUE) 
      
      ggsave(file.path(out_dir, paste0(trait, "_Origin_Plot.svg")), op_plot, width = 7, height = 5,dpi = 600)
    }
  }
}

# 3. EXPORT FINAL MASTER TABLE (Wide-by-Model with Deltas) AND Formatted word document. 

if(length(master_results_list) > 0) {
  
  df_base <- bind_rows(master_results_list) %>% 
    group_by(Trait, Model, Term) %>% slice_tail(n = 1) %>% ungroup()
  
  df_var <- df_base %>%
    group_by(Trait, Model) %>%
    mutate(
      Total_Var = sum(Variance, na.rm = TRUE),
      Pct = Variance / Total_Var 
    ) %>% 
    ungroup()
  
  df_logl <- df_base %>%
    distinct(Trait, Model, LogL) %>%
    mutate(Term = "LogL", Variance = LogL, Pct = NA) %>%
    select(Trait, Model, Term, Variance, Pct)
  
  df_combined <- bind_rows(df_logl, df_var %>% select(Trait, Model, Term, Variance, Pct))
  
  df_wide <- df_combined %>%
    pivot_wider(
      names_from = Model, 
      values_from = c(Variance, Pct),
      names_sep = "_"
    )
  
  expected_cols <- c("Variance_Design", "Variance_Design+", "Variance_Spatial AR1",
                     "Pct_Design", "Pct_Design+", "Pct_Spatial AR1")
  for(col in expected_cols) {
    if(!col %in% colnames(df_wide)) df_wide[[col]] <- NA
  }
  
  df_calc <- df_wide %>%
    mutate(
      `Delta_Design+` = `Variance_Design+` - Variance_Design,
      `Delta_Spatial` = `Variance_Spatial AR1` - `Variance_Design+`,
      `Pct_Delta_Design+` = `Pct_Design+` - Pct_Design,
      `Pct_Delta_Spatial` = `Pct_Spatial AR1` - `Pct_Design+`
    ) %>%
    select(
      Trait, Term,
      Variance_Design, `Variance_Design+`, `Variance_Spatial AR1`,
      `Delta_Design+`, `Delta_Spatial`,
      Pct_Design, `Pct_Design+`, `Pct_Spatial AR1`,
      `Pct_Delta_Design+`, `Pct_Delta_Spatial`
    ) %>%
    arrange(Trait, desc(Term == "LogL"), Term)
  
  trait_list <- split(df_calc, df_calc$Trait)
  spaced_list <- lapply(names(trait_list), function(tr) {
    blank_row <- df_calc[1, ]
    blank_row[1, ] <- NA
    blank_row$Trait <- paste(">>> TRAIT:", tr, "<<<")
    bind_rows(blank_row, trait_list[[tr]])
  })
  
  final_table <- bind_rows(spaced_list)
  
  df_long_export <- df_base %>%
    mutate(Trial = project_name) %>%
    select(Trial, Trait, Model, Term, Variance, LogL) 
  
  write.csv(df_long_export, file.path(out_dir, "All_Traits_Variance_Summary_Long.csv"), row.names = FALSE, na = "")
  cat("\nSUCCESS! Machine-readable Long CSV saved to 'Analyses'.\n")
  
  doc <- read_docx()
  trait_list <- split(df_calc, df_calc$Trait)
  
  for (tr in names(trait_list)) {
    df_formatted <- trait_list[[tr]] %>%
      select(-Trait) %>% 
      mutate(across(starts_with("Pct"), ~ ifelse(is.na(.), "", sprintf("%.1f%%", . * 100)))) %>%
      mutate(across(where(is.numeric), ~ ifelse(is.na(.), "", sprintf("%.3f", .))))
    
    ft <- flextable(df_formatted) %>%
      theme_booktabs() %>%
      fontsize(size = 9, part = "all") %>%
      padding(padding = 3, part = "all") %>%
      set_header_labels(
        Variance_Design = "Var\nDesign",
        `Variance_Design+` = "Var\nDesign+",
        `Variance_Spatial AR1` = "Var\nSpatial",
        `Delta_Design+` = "Delta\nDesign+",
        Delta_Spatial = "Delta\nSpatial",
        Pct_Design = "%\nDesign",
        `Pct_Design+` = "%\nDesign+",
        `Pct_Spatial AR1` = "%\nSpatial",
        `Pct_Delta_Design+` = "% Delta\nDesign+",
        Pct_Delta_Spatial = "% Delta\nSpatial"
      ) %>%
      autofit() %>%
      align(j = "Term", align = "left", part = "all") %>%
      align(j = 2:ncol(df_formatted), align = "center", part = "all")
    
    doc <- doc %>%
      body_add_par(paste("Trait:", tr), style = "heading 2") %>% 
      body_add_flextable(value = ft, align = "left") %>%
      body_add_par("", style = "Normal") 
  }
  
  doc <- doc %>% body_end_section_landscape()
  print(doc, target = file.path(out_dir, "All_Traits_Variance_Summary_Formatted.docx"))
  cat("SUCCESS! Formatted Landscape Word tables (split by trait) saved to 'Analyses'.\n")
  
} else {
  cat("\n[!] WARNING: No results were captured to save.\n")
}

### 4. EXPORT MASTER FIXED EFFECTS TABLES #### 
if(length(master_fixed_list) > 0) {
  
  origin_mapping <- c(
    "1"="Ro_North_HG", "2"="Ro_North_Unk", "3"="Ro_North_WC", 
    "4"="Ro_South_Unk", "5"="Ro_HG", "6"="Ro_Ledmore", 
    "7"="Ro_WC", "8"="Ro_Filler"
  )
  
  fixed_export <- bind_rows(master_fixed_list) %>%
    mutate(
      Trial = project_name,
      # NEW: Safely label the levels based on the specific Term
      Level_Name = case_when(
        tolower(Term) == "origin" ~ as.character(origin_mapping[as.character(Level)]),
        tolower(Term) == "edge" ~ paste("Edge", Level),
        tolower(Term) %in% c("distright", "distleft", "distedge") ~ "Slope",
        TRUE ~ as.character(Level)
      )
    ) %>%
    select(Trial, Trait, Model, Term, Level, Level_Name, Estimate, SE, T_value, Wald_P_Value) %>%
    arrange(Trait, Model, Term, Level)
  
  write.csv(fixed_export, file.path(out_dir, "All_Traits_Fixed_Effects.csv"), row.names = FALSE, na = "")
  cat("SUCCESS! Exported Fixed Effects (CSV) to 'Analyses'.\n")
  
  # --- NEW: PUBLICATION-READY WORD TABLE FOR FIXED EFFECTS ---
  doc_fixed <- read_docx()
  fixed_trait_list <- split(fixed_export, fixed_export$Trait)
  
  for (tr in names(fixed_trait_list)) {
    df_formatted_fixed <- fixed_trait_list[[tr]] %>%
      select(Model, Term, Level_Name, Estimate, SE, T_value, Wald_P_Value) %>%
      mutate(across(c(Estimate, SE, T_value), ~ ifelse(is.na(.), "", sprintf("%.3f", .))))
    
    ft_fixed <- flextable(df_formatted_fixed) %>%
      theme_booktabs() %>%
      fontsize(size = 9, part = "all") %>%
      padding(padding = 3, part = "all") %>%
      set_header_labels(
        Level_Name = "Level",
        T_value = "T-Value",
        Wald_P_Value = "Wald p-value"
      ) %>%
      autofit() %>%
      merge_v(j = c("Model", "Term")) %>% # Merges duplicate names vertically for a super clean look
      align(j = 1:3, align = "left", part = "all") %>%
      align(j = 4:ncol(df_formatted_fixed), align = "center", part = "all")
    
    doc_fixed <- doc_fixed %>%
      body_add_par(paste("Trait:", tr, "- Fixed Effects"), style = "heading 2") %>% 
      body_add_flextable(value = ft_fixed, align = "left") %>%
      body_add_par("", style = "Normal") 
  }
  
  doc_fixed <- doc_fixed %>% body_end_section_portrait()
  print(doc_fixed, target = file.path(out_dir, "All_Traits_Fixed_Effects_Formatted.docx"))
  cat("SUCCESS! Formatted Portrait Word tables for Fixed Effects saved to 'Analyses'.\n")
}

