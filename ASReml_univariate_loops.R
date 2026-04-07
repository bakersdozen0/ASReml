######## ASReml R is not cooporating in the license activation arena; let's try using R as the text parser
## and asreml SA as the engine: 

library(here)
library(asreml)
library(tidyverse)
library(ggplot2)
library(patchwork)

# 1. Setup
# Read raw data - ensure trait columns are treated as numeric
raw_data <- read.csv("Brecon_59_S.csv", stringsAsFactors = FALSE)

traits_to_test <- c("Dbhob_28")
# other testers; "Dbhob_28", "Dbhub_28", "Drs_28", "Bte_28", "Dre_28"
as_file <- "Brecon_59_S.as"

asreml_path <- "C:/Program Files/ASReml4/bin/asreml.exe"
asreml_exe <- shQuote(asreml_path) # Properly handles the space in "Program Files"

# --- DYNAMIC REGEX TRAIT SCRAPER WITH GUARD RAIL ---
as_text <- paste(readLines(as_file), collapse = " ")
# Find all potential traits (AlphaNumeric_Digits)
found_traits <- unique(unlist(str_extract_all(as_text, "\\b[A-Za-z0-9]+_[0-9]+\\b")))

# GUARD RAIL: Only keep things that are actually columns in your CSV
traits_to_test <- found_traits[found_traits %in% colnames(raw_data)]

cat("Automated Discovery: Found", length(found_traits), "potential matches.\n")
cat("Guard Rail: Proceeding with", length(traits_to_test), "traits found in CSV.\n")

#Save all ASReml outputs to separate directory; 
out_dir <- "Analyses"
if(!dir.exists(out_dir)) {
  dir.create(out_dir)
  cat("Created 'Analyses' folder for outputs.\n")
}

base_name <- gsub("\\.as$", "", as_file) 
out_asr <- paste0(base_name, ".asr")
out_yht <- paste0(base_name, ".yht")
out_sln <- paste0(base_name, ".sln") # Added this just in case you want the solutions too!

# Models to loop through (matches your !PART numbers)
models_to_run <- c("1" = "Design", "2" = "Design+", "3" = "Spatial AR1")

master_results <- list()

# --- NEW: Calculate Block Boundaries for the Black Outlines ---
block_bounds <- raw_data %>%
  filter(!is.na(Ppos) & !is.na(Prow) & !is.na(Block)) %>%
  group_by(Block) %>%
  summarize(
    xmin = min(as.numeric(Ppos), na.rm = TRUE) - 0.5,
    xmax = max(as.numeric(Ppos), na.rm = TRUE) + 0.5,
    ymin = min(as.numeric(Prow), na.rm = TRUE) - 0.5,
    ymax = max(as.numeric(Prow), na.rm = TRUE) + 0.5
  )

ppgmap_colors <- c("darkblue", "blue", "cyan", "green", "yellow", "orange", "red", "darkred")



### EXECUTE::: ####
for (trait in traits_to_test) {
  cat("\n========================================\nProcessing:", trait, "\n")
  
  raw_data[[trait]] <- suppressWarnings(as.numeric(as.character(raw_data[[trait]])))
  res_plots <- list(); sol_plots <- list(); trait_variance_list <- list()
  
  # --- PANEL 1: RAW DATA ---
  raw_map <- ggplot(raw_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = .data[[trait]])) +
    geom_tile(color = "black", size = 0.05) + 
    geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
    scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90") +
    scale_y_reverse() + theme_void() + coord_fixed() + labs(title = "1. Raw Data")
  res_plots[[1]] <- raw_map; sol_plots[[1]] <- raw_map
  
  # --- RUN MODELS ---
  for (part in names(models_to_run)) {
    model_name <- models_to_run[[part]]
    cat("  -> Running", model_name, "... ")
    
    suppressWarnings(file.remove(out_asr, out_yht, out_sln))
    suppressWarnings(file.remove(list.files(pattern = paste0("^", base_name, "\\.(msv|veo|ask|tmp|tsv)$"))))
    
    # SILENCED COMMAND
    command <- paste(asreml_exe, "-n", as_file, part, trait, "1.2 > NUL 2>&1")
    shell(command, wait = TRUE)
    
    if(!file.exists(out_asr)) { cat("FAILED\n"); next }
    asr_lines <- readLines(out_asr)
    if(!any(grepl("LogL Converged", asr_lines))) { cat("Converge Fail\n"); next }
    cat("Success\n")
    
    # A. Scrape Variances
    logl_val <- as.numeric(str_extract(tail(grep("LogL=", asr_lines, value = TRUE), 1), "(?<=LogL=\\s)[-0-9.]+"))
    terms <- c("Block", "SubBlock", "Block\\.SubBlock", "Block\\.Prow", "Block\\.Ppos", "Family_id", "uni\\(Crosstype,2\\)", "units")
    for(t in terms) {
      line <- grep(paste0("^\\s*", t, "\\s+"), asr_lines, value = TRUE)
      if(length(line) > 0) {
        val <- as.numeric(unlist(strsplit(trimws(line[1]), "\\s+"))[5])
        df_row <- data.frame(Trait=trait, Model=model_name, Term=gsub("\\\\", "", t), Variance=val)
        trait_variance_list[[length(trait_variance_list) + 1]] <- df_row
        master_results_list[[length(master_results_list) + 1]] <- df_row %>% mutate(LogL=logl_val)
      }
    }
    
    # B. Process Maps
    sln <- read.table(out_sln, skip = 1, fill = TRUE, stringsAsFactors = FALSE)
    colnames(sln) <- c("Term", "Level", "Estimate", "SE")
    yht <- read.table(out_yht, skip = 1)
    colnames(yht) <- c("Record", "Yhat", "Residual", "Hat")
    
    map_data <- raw_data %>% mutate(Resid = ifelse(is.na(.data[[trait]]), NA, yht$Residual))
    
    if (part == "1") {
      block_eff <- sln %>% filter(Term == "Block") %>% select(Level, Estimate)
      map_data <- map_data %>% left_join(block_eff, by = c("Block" = "Level")) %>%
        mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, Estimate))
      sol_title <- "2. Design (Block Solutions)"
    } else if (part == "2") {
      b_eff <- sln %>% filter(Term == "Block") %>% select(Level, Estimate); r_eff <- sln %>% filter(Term == "Block.Prow") %>% select(Level, Estimate); p_eff <- sln %>% filter(Term == "Block.Ppos") %>% select(Level, Estimate)
      map_data <- map_data %>% 
        left_join(b_eff, by = c("Block" = "Level")) %>% rename(B_est = Estimate) %>%
        left_join(r_eff, by = c("Prow" = "Level")) %>% rename(R_est = Estimate) %>%
        left_join(p_eff, by = c("Ppos" = "Level")) %>% rename(P_est = Estimate) %>%
        mutate(D_Sum = replace_na(B_est,0) + replace_na(R_est,0) + replace_na(P_est,0)) %>%
        mutate(PlotVal = ifelse(is.na(.data[[trait]]), NA, D_Sum))
      sol_title <- "3. Design+ (Blk+Row+Col Solutions)"
    } else {
      map_data <- map_data %>% 
        mutate(Trend = ifelse(is.na(.data[[trait]]), NA, yht$Yhat - mean(yht$Yhat, na.rm=TRUE))) %>%
        rename(PlotVal = Trend)
      sol_title <- "4. Total Spatial Trend Correction"
    }
    
    sol_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = PlotVal)) +
      geom_tile(color = "black", size = 0.05) + 
      geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90") +
      scale_y_reverse() + theme_void() + coord_fixed() + labs(title = sol_title)
    
    res_plots[[as.numeric(part) + 1]] <- ggplot(map_data, aes(x = as.numeric(Ppos), y = as.numeric(Prow), fill = Resid)) +
      geom_tile(color = "black", size = 0.05) + 
      geom_rect(data = block_bounds, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill=NA, color="black", size=0.8, inherit.aes=FALSE) +
      scale_fill_gradientn(colors = ppgmap_colors, na.value = "grey90") +
      scale_y_reverse() + theme_void() + coord_fixed() + labs(title = paste0(as.numeric(part)+1, ". ", model_name, " Res"))
    
    file.copy(out_asr, file.path(out_dir, paste0(trait, "_", model_name, ".asr")), overwrite = TRUE)
  }
  
  if(length(res_plots) == 4) {
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_RESIDUALS.png")), (res_plots[[1]] | res_plots[[2]]) / (res_plots[[3]] | res_plots[[4]]), width = 12, height = 10)
    ggsave(file.path(out_dir, paste0(trait, "_4Panel_SOLUTIONS.png")), (sol_plots[[1]] | sol_plots[[2]]) / (sol_plots[[3]] | sol_plots[[4]]), width = 12, height = 10)
  }
  
  if(length(trait_variance_list) > 0) {
    trait_df <- bind_rows(trait_variance_list) %>% 
      group_by(Model) %>% mutate(Pct = Variance / sum(Variance) * 100) %>% ungroup()
    trait_df$Model <- factor(trait_df$Model, levels = c("Design", "Design+", "Spatial AR1"))
    bp <- ggplot(trait_df, aes(x = Model, y = Pct, fill = Term)) +
      geom_col(color = "black") + theme_minimal() + scale_fill_brewer(palette = "Set3") +
      labs(title = paste("Variance Components:", trait), y = "% Variance") 
    ggsave(file.path(out_dir, paste0(trait, "_VC_Barplot.png")), bp, width = 7, height = 6)
  }
}

if(length(master_results_list) > 0) {
  final_table <- bind_rows(master_results_list) %>% 
    # 1. Clean up duplicates so pivot_wider doesn't create 'list-cols'
    # We group by the key columns and take only the last (final) value found
    group_by(Trait, Model, Term) %>%
    slice_tail(n = 1) %>% 
    ungroup() %>%
    # 2. Pivot wide so terms are columns and traits/models are rows
    pivot_wider(names_from = Term, values_from = Variance) %>% 
    # 3. Reorder for the 'DesignSpatial' layout
    select(Trait, Model, LogL, everything()) %>%
    arrange(Trait, Model)
  
  # Write the file to your 'Analyses' folder
  write.csv(final_table, file.path(out_dir, "All_Traits_Variance_Summary.csv"), row.names = FALSE)
  
  cat("\nSUCCESS! All data and plots for the discovered traits are in the 'Analyses' folder.\n")
} else {
  cat("\n[!] WARNING: No results were captured to save.\n")
}
