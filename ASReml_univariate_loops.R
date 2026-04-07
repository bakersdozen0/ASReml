

#### install; ####
#install.packages("biometryassist")
#library(biometryassist)

install_asreml()

packages <- c("data.table", "ggplot2", "jsonlite") 
install.packages(packages) 


### alternative; download zip file and point this command at it; 
# NB: I had to do it this way ; 

install.packages("C:/Users/james.baker/R projects/asreml_4.2.0.392.zip", repos=NULL,type="win.binary")
# activate
library(asreml)
asreml.license.activate()

######## ASReml R is not cooporating in the license activation arena; let's try using R as the text parser
## and asreml SA as the engine: 

# 1. Setup
traits_to_test <- c("Dre_28", "Dbhob_28", "Dbhub_28", "Drs_28", "Bte_28")
as_file <- "Brecon_59_S.as"

base_name <- gsub("\\.as$", "", as_file) 
out_asr <- paste0(base_name, ".asr")
out_yht <- paste0(base_name, ".yht")
out_sln <- paste0(base_name, ".sln") # Added this just in case you want the solutions too!

# 2. The exact executable path with escaped double quotes
asreml_exe <- "\"C:/Program Files/ASReml4/bin/asreml.exe\"" 

# 3. Run the loop
for (trait in traits_to_test) {
  
  cat("\n----------------------------------------\n")
  cat("Running ASReml for trait:", trait, "\n")
  
  # Construct the command
  command <- paste(asreml_exe, "-n", as_file, "3", trait, "1.2")
  
  # Execute the command (silently wait for it to finish)
  system(command, wait = TRUE)
  
  # Check for and rename the results files
  if(file.exists(out_asr)) {
    new_asr_name <- paste0(trait, "_results.asr")
    file.rename(out_asr, new_asr_name)
    cat("Success! Saved as", new_asr_name, "\n")
  } else {
    cat("ERROR: Did not create the .asr file for", trait, "\n")
  }
  
  if(file.exists(out_yht)) {
    file.rename(out_yht, paste0(trait, "_predictions.yht"))
  }
  
  if(file.exists(out_sln)) {
    file.rename(out_sln, paste0(trait, "_solutions.sln"))
  }
}

cat("\nAll traits finished running!\n")