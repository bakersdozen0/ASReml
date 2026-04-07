# ASReml Spatial Pipeline

This repository contains an automated R pipeline for executing and visualizing univariate spatial analyses using ASReml Standalone. 

The script dynamically reads an ASReml `.as` file to identify traits, runs a tiered model progression (Design, Design+, and Spatial AR1), and extracts the relevant variance components and spatial mapping files for visualization.

## Prerequisites

To run this pipeline, you will need:
* **ASReml Standalone 4.2+** installed locally with an active license.
* **R and RStudio**.
* The following R packages: `tidyverse`, `ggplot2`, `patchwork`, `here`.

## File Structure & `.gitignore`

**Important:** Do not commit proprietary data or massive output files to GitHub. Make sure you have a `.gitignore` file in this directory containing the following:
```text
*.csv
*.as
*.asr
*.yht
*.sln
Analyses/
