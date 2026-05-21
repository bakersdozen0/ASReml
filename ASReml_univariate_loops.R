# 
# 0. CONTROL PANEL (Change these for your specific project) #### 
# 

trial_folder  <- "C:/Users/james.baker/Forest Research/TW CBC-TBA-NextGenBritishConifers - Share/Sitka/High GCA Fullsib P85-P87 experiments/Brecon 8"
project_name  <- "Brecon_8_S"
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

# traits_to_test <-c("Ht_06")

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

# --- MAPPING CALCULATIONS & THEMES ---
## calculate min and max values for legend: 
force_breaks <- function(x) { 
  min_val <- min(x, na.rm = TRUE)
  max_val <- max(x, na.rm = TRUE)
  if (!is.finite(min_val) || !is.finite(max_val)) return(c(0, 1)) # Failsafe for all NAs
  seq(min_val, max_val, length.out = 5) 
}
format_labels <- function(x) {
  round(x, 2) # Rounds the exact min/max to 2 decimal places so the legend stays clean
}

# Calculate the geometric center for labels
block_labels <- raw_data %>%
  filter(!is.na(Ppos) & !is.na(Prow) & !is.na(Block) & 
           as.character(Block) != "0" & 
           as.character(Block) != "" & 
           as.character(Block) != ".") %>% 
  group_by(Block) %>%
  summarize(
    x_mid = mean(as.numeric(Ppos), na.rm = TRUE),
    y_mid = mean(as.numeric(Prow), na.rm = TRUE)
  )

# Calculate Internal Segments
h_segments <- raw_data %>%
  arrange(Prow, Ppos) %>%
  mutate(next_block = lead(Block)) %>%
  filter(Block != next_block & Prow == lead(Prow) & !is.na(Block)) %>%
  mutate(x = Ppos + 0.5, xend = Ppos + 0.5, y = Prow - 0.5, yend = Prow + 0.5)

v_segments <- raw_data %>%
  arrange(Ppos, Prow) %>%
  mutate(next_block = lead(Block)) %>%
  filter(Block != next_block & Ppos == lead(Ppos) & !is.na(Block)) %>%
  mutate(x = Ppos - 0.5, xend = Ppos + 0.5, y = Prow + 0.5, yend = Prow + 0.5)

# Calculate Outer Trial Border
trial_border <- data.frame(
  xmin = min(raw_data$Ppos, na.rm = TRUE) - 0.5,
  xmax = max(raw_data$Ppos, na.rm = TRUE) + 0.5,
  ymin = min(raw_data$Prow, na.rm = TRUE) - 0.5,
  ymax = max(raw_data$Prow, na.rm = TRUE) + 0.5
)

# Standard Theme
theme_trial <- theme_void() + 
  theme(
    plot.title = element_text(size = 20, face = "bold", margin = margin(b = 10)),
    plot.subtitle = element_text(size = 14, color = "grey30"),
    legend.key.height = unit(1.5, "cm"),
    legend.text = element_text(size = 12)
  )

# Shared Mapping Layers
shared_layers <- list(
  geom_segment(data = h_segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 0.8, inherit.aes = FALSE),
  geom_segment(data = v_segments, aes(x = x, xend = xend, y = y, yend = yend), color = "black", size = 0.8, inherit.aes = FALSE),
  geom_rect(data = trial_border, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=1.5, inherit.aes=FALSE),
  geom_label(data = block_labels, aes(x = x_mid, y = y_mid, label = Block), inherit.aes = FALSE, size = 6, fontface = "bold", fill = alpha("white", 0.9), label.size = NA)
)

ppgmap_colors <- c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred")

# 
# 2. Execution section ####
# 
for (trait in traits_to_test) {
  cat("\n========================================\nProcessing:", trait, "\n")
  raw_data[[trait]] <- suppressWarnings(as.numeric(as.character(raw_data[[trait]])))
  
  trait_variance_list <- list()
  master_fixed_list <- list()
  
  # Initialize plots for this trait
  blank_plot <- ggplot() + theme_void() + 
    annotate("text", x = 0.5, y = 0.5, label = "Model Did Not Converge", color = "darkred", fontface = "italic", size = 5) +
    coord_fixed(xlim = c(0, 1), ylim = c(0, 1)) 
  
  res_plots <- list(blank_plot, blank_plot, blank_plot, blank_plot)
  sol_plots <- list(blank_plot, blank_plot, blank_plot, blank_plot)
  
  design_Ve <- NA 
  design_logL <- NA
  
   # --- PANEL 1: RAW DATA ---
  res_plots[[1]] <- ggplot(raw_data, aes(x = Ppos, y = Prow)) +
    geom_tile(aes(fill = .data[[trait]]), color = "lightgray", size = 0.2) + 
    scale_fill_gradientn(colors = ppgmap_colors, na.value = "white", breaks = force_breaks, labels = format_labels) +
    shared_layers + # Now included in every plot
    theme_trial + scale_y_reverse() + coord_fixed() + labs(title = "1. Raw Data")
  sol_plots[[1]] <- res_plots[[1]]
  
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
    
    # --- NEW: DIAGNOSTIC WARNING SCRAPER ---
    cat("      [Diagnostics]:")
    
    # 1. End-of-file warnings (e.g., Parameters Not Converged)
    finish_line <- grep("Finished:", asr_lines, value = TRUE)
    if(length(finish_line) > 0 && grepl("Warning|Error", finish_line[1], ignore.case = TRUE)) {
      warn_msg <- trimws(sub(".*Finished:.*(Warning:|Error:)", "\\1", finish_line[1], ignore.case=TRUE))
      cat("\n        -> RUN MSG: ", warn_msg)
    }
    
    # 2. Design Singularities
    sing_lines <- grep("singularities detected", asr_lines, ignore.case = TRUE, value = TRUE)
    if(length(sing_lines) > 0) {
      cat("\n        ->", trimws(sing_lines[1]))
    }
    
    # 3. Variance Component Boundaries (B, ?, S)
    mt_start <- grep("Model_Term", asr_lines, ignore.case = TRUE)
    wald_start <- grep("Wald F statistics", asr_lines, ignore.case = TRUE)
    
    if(length(mt_start) > 0) {
      end_idx <- ifelse(length(wald_start) > 0, wald_start[1] - 1, length(asr_lines))
      vc_table <- asr_lines[mt_start[1]:end_idx]
      
      # Regex: Look for rows ending with B, ?, or S (accounting for trailing spaces)
      bad_vc <- grep("\\s+[B\\?S]\\s*$", vc_table, value = TRUE)
      
      if(length(bad_vc) > 0) {
        cat("\n        -> BOUNDARY/SINGULAR COMPONENTS DETECTED:")
        for(b in bad_vc) cat("\n             ", trimws(b))
      }
    }
    cat("\n")
    
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
    
    terms <- c("Block", "SubBlock", "Block\\.SubBlock", "Block\\.Prow", "Block\\.Ppos", "Prow", "Ppos", "Family_id", "Family_name", "uni\\(Crosstype,2\\)", "units", "Residual", "Plot","Tree") 
    
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
      # 1. Baseline Model
      design_logL <- logl_val
      d_logL <- "Base" # Set text for the subtitle instead of a number
      ve_row <- current_vars_df %>% filter(Term == "Independent Error")
      design_Ve <- if(nrow(ve_row) > 0) ve_row %>% pull(Variance) %>% .[1] else 1
      
    } else if (part == "2") {
      # 2. Design+ Model (Compared to Baseline Design)
      designplus_logL <- logl_val
      d_logL <- round(logl_val - design_logL, 2)
      
    } else if (part == "3") {
      # 3. Spatial AR1 Model (Compared to Design+)
      d_logL <- round(logl_val - designplus_logL, 2)
    }
    
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
    
    
    # --- EXTRACT ALL FIXED EFFECTS (ORIGIN & SPATIAL COVARIATES) & WALD P-VALUE ---
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
      mutate(Resid = Residual)
    
    if (part == "1" || part == "2") {
      b_eff <- sln %>% filter(tolower(Term) == "block") %>% select(Level, Estimate)
      # Use a robust regex to catch "Subblock" or "SubBlock" in the raw data
      sub_col_name <- grep("subblock", colnames(raw_data), ignore.case = TRUE, value = TRUE)
      
      sub_eff <- sln %>% filter(tolower(Term) %in% c("block.subblock", "subblock")) %>% select(Level, Estimate)
      p_eff <- sln %>% filter(tolower(Term) %in% c("plot", "block.plot")) %>% select(Level, Estimate)
      
      map_data <- map_data %>% 
        mutate(
          Block_key = as.character(Block),
          # Safely build the Sub_key only if the column actually exists
          Sub_key = if(length(sub_col_name) > 0) paste(Block, .data[[sub_col_name[1]]], sep=".") else NA,
          Plot_key = as.character(Plot)
        ) %>%
        left_join(b_eff, by = c("Block_key" = "Level")) %>% rename(B_est = Estimate) %>%
        left_join(sub_eff, by = c("Sub_key" = "Level")) %>% rename(Sub_est = Estimate) %>%
        left_join(p_eff, by = c("Plot_key" = "Level")) %>% rename(P_est = Estimate) %>%
        mutate(D_Sum = replace_na(B_est,0) + replace_na(Sub_est,0) + replace_na(P_est,0))
      
      
      # ---  DYNAMIC EDGE SCRAPER FOR VISUAL MAPS ---
      env_terms <- c("Edge", "Distright", "Distleft", "Distedge")
      for (e_term in env_terms) {
        term_idx <- tolower(sln$Term) == tolower(e_term)
        if (any(term_idx) && e_term %in% colnames(map_data)) {
          t_sln <- sln[term_idx, ]
          if (nrow(t_sln) > 1) { # Factor logic
            temp_join <- map_data %>%
              mutate(Join_Key = as.character(.data[[e_term]])) %>%
              left_join(t_sln %>% select(Level, Estimate), by = c("Join_Key" = "Level"))
            map_data$D_Sum <- map_data$D_Sum + replace_na(temp_join$Estimate, 0)
          } else if (nrow(t_sln) == 1) { # Covariate logic
            slope <- t_sln$Estimate[1]
            raw_val <- suppressWarnings(as.numeric(as.character(map_data[[e_term]])))
            map_data$D_Sum <- map_data$D_Sum + replace_na(raw_val * slope, 0)
          }
        }
      }

      
      map_data <- map_data %>% 
        mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, D_Sum))
      sol_title <- "3. Design+ (Blk+Row+Col+Edge Solutions)"
      
    } else {
      map_data <- map_data %>% 
        mutate(
          PlotVal = ifelse(is.na(.data[[trait]]), NA, Resid), 
          Resid = ifelse(is.na(.data[[trait]]), NA, yht$Yhat - mean(yht$Yhat, na.rm=TRUE)) 
        )
      sol_title <- "4. Spatial Effects Correction"
    }
    
   # Define the plot
    # Notice the explicit inherit.aes = FALSE in the border layer
    sol_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = Ppos, y = Prow)) +
      geom_tile(aes(fill = PlotVal), color = "lightgray", size = 0.2) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "white", breaks = force_breaks, labels = format_labels) +
      shared_layers + # Standard template
      theme_trial + scale_y_reverse() + coord_fixed() + labs(title = model_name, subtitle = metrics_subtitle)
    
    res_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = Ppos, y = Prow)) +
      geom_tile(aes(fill = Resid), color = "lightgray", size = 0.2) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "white", breaks = force_breaks, labels = format_labels) +
      shared_layers + # Standard template
      theme_trial + scale_y_reverse() + coord_fixed() + labs(title = paste(model_name, "Resid"), subtitle = metrics_subtitle)
    
    # --- 4. CLEAN THE SANDBOX FOR THE NEXT LOOP ---
    file.remove(sandbox_asr)
    if(file.exists(sandbox_yht)) file.remove(sandbox_yht)
    if(file.exists(sandbox_sln)) file.remove(sandbox_sln)
  }
  
  # NEW: Dynamic 4-Panel Plot Sizing (Tile-Based Scaling)
  if(length(res_plots) == 4) {
    master_title <- paste(project_name, "- Trait:", trait)
    
    # --- DYNAMIC ASPECT RATIO CALCULATION ---
    trial_cols <- max(raw_data$Ppos, na.rm = TRUE) - min(raw_data$Ppos, na.rm = TRUE) + 1
    trial_rows <- max(raw_data$Prow, na.rm = TRUE) - min(raw_data$Prow, na.rm = TRUE) + 1
    
    # 1. Define physical size per tree (0.15 inches gives a great text-to-plot balance)
    inches_per_tree <- 0.15 
    
    # 2. Calculate dimensions for ONE plot
    single_plot_width <- trial_cols * inches_per_tree
    single_plot_height <- trial_rows * inches_per_tree
    
    # 3. Calculate Total Canvas (2x2 grid + padding for titles and legends)
    # Add ~2.5 inches width for the legends, and ~2.5 inches height for main titles
    dyn_width <- (single_plot_width * 2) + 2.5
    dyn_height <- (single_plot_height * 2) + 2.5
    
    # 4. Failsafe: Prevent the canvas from getting so small that standard text breaks
    dyn_width <- max(dyn_width, 10)
    dyn_height <- max(dyn_height, 7)
    # -----------------------------------------
    
    res_assembled <- (res_plots[[1]] | res_plots[[2]]) / (res_plots[[3]] | res_plots[[4]]) +
      plot_annotation(title = master_title, theme = theme(plot.title = element_text(size = 18, face = "bold")))
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_RESIDUALS.svg")), res_assembled, width = dyn_width, height = dyn_height, dpi = 600)
    
    sol_assembled <- (sol_plots[[1]] | sol_plots[[2]]) / (sol_plots[[3]] | sol_plots[[4]]) +
      plot_annotation(title = master_title, theme = theme(plot.title = element_text(size = 18, face = "bold")))
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_SOLUTIONS.svg")), sol_assembled, width = dyn_width, height = dyn_height, dpi = 600)
  }
  
  # --- BARPLOT SECTION ---
  if(length(trait_variance_list) > 0) {
    
    # Official terms list. 
    expected_terms <- c("Block", "SubBlock", "Block.SubBlock", "Block.Prow", 
                        "Block.Ppos", "Prow", "Ppos", "Plot", "Tree", "Family_id", 
                        "Family_name", "uni(Crosstype,2)", "Spatial Variance", "Independent Error")
    
    trait_df <- bind_rows(trait_variance_list) %>% 
      group_by(Model, Term) %>% 
      slice_tail(n = 1) %>%  # Removes the duplicates
      ungroup() %>%
      filter(Term %in% expected_terms) %>% # Uses the exact terms the scraper generated
      group_by(Model) %>%
      mutate(
        Total_Var_In_Model = sum(Variance, na.rm = TRUE),
        Pct_Var = (Variance / Total_Var_In_Model) * 100
      ) %>% 
      ungroup()
    
    # Force the ordering so colors stay mapped logically
    trait_df$Model <- factor(trait_df$Model, levels = c("Design", "Design+", "Spatial AR1"))
    trait_df$Term <- factor(trait_df$Term, levels = expected_terms)
    
    bp <- ggplot(trait_df, aes(x = Model, y = Pct_Var, fill = Term)) +
      geom_col(color = "black", size = 0.2) + 
      geom_hline(yintercept = 100, linetype = "dashed", color = "black", size = 0.5) +
      theme_minimal() + 
      scale_fill_brewer(palette = "Set3") +
      labs(title = paste("Variance Components:", trait), 
           y = "% of Total Variance") +
      theme(legend.position = "right")
    
    ggsave(file.path(out_dir, paste0(trait, "_VC_Barplot.png")), bp, width = 8, height = 6, dpi = 600)
  }
  # --- ORIGIN SOLUTIONS PLOT (PER TRAIT) ---
  if(length(master_fixed_list) > 0) {
    # 1. Map labels
    origin_mapping <- c(
      "1"="Ro_North_HG", "2"="Ro_North_Unk", "3"="Ro_North_WC", 
      "4"="Ro_South_Unk", "5"="Ro_HG", "6"="Ro_Ledmore", 
      "7"="Ro_WC", "8"="Ro_Filler"
    )
    
    # 2. Filter for ONLY the current trait and Spatial AR1 results
    trait_origin_df <- bind_rows(master_fixed_list) %>% 
      filter(Trait == trait & Model == "Spatial AR1" & Term == "Origin") %>%
      mutate(Origin_Name = factor(origin_mapping[as.character(Level)], levels = origin_mapping))
    
    # 3. Only plot if we have data for this trait
    if(nrow(trait_origin_df) > 0) {
      origin_p <- trait_origin_df$Wald_P_Value[1]
      
      op_plot <- ggplot(trait_origin_df, aes(x = Origin_Name, y = Estimate)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        geom_errorbar(aes(ymin = Estimate - SE, ymax = Estimate + SE), width = 0.2, color = "darkblue") +
        geom_point(size = 3, color = "darkblue") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(
          title = paste("Origin Solutions (Spatial AR1):", trait),
          subtitle = paste("Wald P-value:", origin_p),
          x = "Origin",
          y = "Solution Estimate (+/- 1 SE)"
        )
      
      ggsave(file.path(out_dir, paste0(trait, "_Origin_Plot.png")), op_plot, width = 7, height = 5, dpi = 600)
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

#
#### formally test for whether or not to include edge effects: ####
# 
library(tidyverse)

# 1. READ DATA DIRECTLY FROM THE FOLDERS
base_path <- file.path(trial_folder, "Analyses", "Baseline", "All_Traits_Variance_Summary_Long.csv")
edge_path <- file.path(trial_folder, "Analyses", "DistEdge", "All_Traits_Variance_Summary_Long.csv")

if(!file.exists(base_path) || !file.exists(edge_path)) stop("Baseline or DistEdge CSVs not found. Check your folder paths.")

base_df <- read.csv(base_path) %>% select(Trait, Model, LogL) %>% distinct()
edge_df <- read.csv(edge_path) %>% select(Trait, Model, LogL) %>% distinct()

# 2. STRICT JOIN (Pairs Baseline with DistEdge for the SAME Model)
comparison_results <- base_df %>%
  rename(LogL_Baseline = LogL) %>%
  inner_join(
    edge_df %>% rename(LogL_DistEdge = LogL), 
    by = c("Trait", "Model") # <--- THIS IS THE FIX. It only joins exact matches.
  ) %>%
  mutate(
    Delta_LogL = LogL_DistEdge - LogL_Baseline,
    LRT_Stat = 2 * Delta_LogL,
    Significant = LRT_Stat > 3.84 # Critical value for df=1 at p=0.05
  )

# 3. PRINT RESULTS
print(comparison_results)
write.csv(comparison_results, file.path(trial_folder, "Analyses", "Model_Comparison_Results.csv"), row.names = FALSE)
