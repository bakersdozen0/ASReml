# 🌲 ASReml-R Spatial Pipeline: The "Distrust & Verify" Protocol

**Overview**
This pipeline automates the sequential execution, extraction, and visualization of ASReml 4 spatial models (`Design`, `Design+`, and `Spatial AR1`) for forestry progeny trials. 

Because spatial matrices are highly sensitive to row-order corruption, and because ASReml output logs are notoriously inconsistent, this script was engineered with a "Zero-Trust" architecture. It does not assume data order will be preserved, it does not assume ASReml will name terms consistently, and it calculates exact component scaling natively.

Here is the exact, step-by-step breakdown of how the data is ingested, scraped, and reassembled.

### 1. Data Ingestion & The "Spatial Lock-In"
* **The Premise:** ASReml `.yht` (residuals/fitted values) files are *blind*. They do not contain `Tree_ID` or `Ppos`/`Prow` identifiers; they only output sequential row records (1 to N) based on the order the data was fed into the `.as` file.
* **The Verification:** To guarantee the `.yht` spatial residuals perfectly map back to the correct physical trees, this script **hard-locks the data index upon ingestion**. It reads the pre-sorted CSV (sorted by `Prow` then `Ppos` via the `sfile` protocol) and instantly stamps every row with an immutable `Record` ID matching that exact data-file sequence. No alphabetical sorting or filtering is permitted prior to spatial joining.

### 2. Execution & Environmental Quarantine
* **The Premise:** Running ASReml from a shared network drive via R often causes pathing errors or dumps temporary files (`.msv`, `.veo`) onto the user's local desktop.
* **The Verification:** The script uses `setwd()` to physically step into the shared network trial folder, executes the `.exe` via a silent Windows shell command (`> NUL 2>&1`), and steps back out. This guarantees all `.asr`, `.yht`, and `.sln` files are written directly to the target network directory with zero local bleed-over. 

### 3. The Scraper Engine (Zero-Trust Text Parsing)
ASReml 4 output tables frequently shift column counts, insert dummy headers, and change terminology across model steps. The scraper engine handles this surgically:
* **Targeting Sigma:** The script explicitly splits the `asr` tables and hard-targets **Column 5 (`Sigma`)** to capture the true variance component. It ignores Column 4 (`Gamma` / variance ratios) entirely, preferring to recalculate percentages natively.
* **Filtering Dummy Text:** ASReml prints a blank text line summarizing the residual degrees of freedom (e.g., `Residual 1680 effects`). The script uses a strict Regex filter (`!grepl("effects")`) to destroy these text artifacts before they can corrupt the math.
* **Unifying Error Terminology:** ASReml names standard errors `Residual`, but renames them `units` when an AR1 model explicitly defines them as an independent random effect. The scraper intercepts both terms and forces them into a unified **`Independent Error`** category. This guarantees that baseline model error (Ve) is mathematically comparable across all models.
* **Capturing AR Parameters:** The spatial autocorrelations are not printed as standard variances. The scraper explicitly hunts for ASReml 4's internal AR codes (`AR_R` and `AR_C`) to safely extract the spatial correlation parameters.

### 4. Spatial Reassembly & Visualization
* **The Premise:** We need to view the raw data, block solutions, design+ solutions, and the final AR1 spatial residuals and fitted surfaces.
* **The Verification:** The script reads the blind `.yht` file and performs a strict relational `left_join` using the immutable `Record` ID established in Step 1. The spatial grid is then drawn using the true `Prow` and `Ppos` coordinates. 

### 5. Baseline Scaling & Master Outputs
To ensure changes in genetic variance and error unmasking are visually and mathematically obvious:
* **The Baseline:** The script captures the `Independent Error` (Ve) of the base **Design Model**. 
* **The Scaling:** Every subsequent variance component in the advanced models (including the AR1 model's inflated total variance due to genetic unmasking) is divided by the Design Model's Ve and multiplied by 100. 
* **The Visuals:** 1. **4-Panel Maps:** Includes delta-LogL relative to the base model, current Ve, and AR correlations directly in the plot subtitles.
    2. **Scaled VC Bar Plots:** Displayed as stacked bars so the expansion of genetic signal (`Family_id`) and the shrinkage of `Independent Error` are instantly quantifiable.
    3. **Master CSV:** Automatically pivoted into a wide-format Excel-ready table, complete with calculated delta-LogL, delta-Variances, and calculated percentages.
