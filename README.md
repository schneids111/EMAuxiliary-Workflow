# EMAuxiliary Workflow: Preparing, Evaluating, and Integrating Auxiliary Variables in EMA Data Analyses

This repository accompanies the paper:

> **Schneider, S., Walentynowicz, M., Toledo, M., Hernandez, R., Junghaenel, D.U., Smyth, J.M., Stone, A.A. (under review). Using auxiliary variables to reduce bias from missing EMA prompts: a practical guide. Advances in Methods and Practices in Psychological Science. [Full citation will be updated upon publication.]**  
> The workflow provides a three-step pipeline for preparing ecological momentary assessment (EMA) datasets, evaluating potential auxiliary variables, and incorporating them into multivariate multilevel models using the **EMAuxiliary** R function (which generates **Blimp** input code for auxiliary-variable modeling).

---

## 📘 Overview of the Three-Step Workflow

| Step | Purpose | Output / Next Input |
|------|----------|---------------------|
| **Step 1 – Prepare EMA Data** (`Step1_EMA_DataPrep.Rmd`) | Helps create a complete long-format EMA dataset where each scheduled prompt—completed or missed—is represented. Add rows for missed prompts, lagged variables, and person-level data. | EMA dataset with auxiliary variables |
| **Step 2 – Evaluate Auxiliary Variables** (`Step2_AuxEval.Rmd`) | Evaluate potential auxiliary variables based on prediction of missingness and correlations with focal variables. | Auxiliary variable diagnostics |
| **Step 3 – Fit Models with EMAuxiliary** (`Step3_ExamplesEMAuxiliary.Rmd`) | Incorporate selected auxiliaries into multilevel models in Blimp. | Full Blimp output and model estimates |

Each component is fully stand-alone and can be run independently, but the workflow is designed to flow from one step to the next.

---

## ⚙️ The EMAuxiliary Function

The core functionality of this workflow is implemented in the EMAuxiliary() function, located in the EMAuxiliary/ folder.

  📄 EMAuxiliary.R -- Generates Blimp syntax for multivariate multilevel models with auxiliary variables
  
  📄 EMAuxiliary Reference Manual -- Provides detailed documentation and trouble shooting tips

The function allows users to specify models using familiar lme4-style syntax and automatically translates these into Blimp input files for estimation.

---

## 🚀 Getting Started

### 1. Installation and Requirements
- **R (≥ 4.2)**  
- **Blimp**: [Download Blimp](https://www.appliedmissingdata.com/blimp)  
- **R interface for Blimp (rblimp)**:
  ```r
  remotes::install_github("blimp-stats/rblimp")
- **EMAuxiliary function**:
  ```r
  source("https://raw.githubusercontent.com/schneids111/EMAuxiliary-Workflow/main/EMAuxiliary/EMAuxiliary.R")


### 2. Using the Tutorials
   
Step 1: Run relevant pieces of the notebook Step1_EMA_DataPrep.Rmd to reconstruct or simulate your EMA dataset.

Step 2: Run Step2_AuxEval.Rmd using your prepared data as input. This will generate detailed information and heatmaps showing how candidate auxiliaries relate to missingness and focal model variables.

Step 3: Use the EMAuxiliary() function (located in EMAuxiliary/EMAuxiliary.R) to generate a Blimp model that includes your selected auxiliaries. See Step3_ExamplesEMAuxiliary.Rmd for worked examples. 
      **Use the EMAuxiliary Reference Manual for guidance.**

---

## 📊 Supplemental Materials (Reproducibility)

The Supplemental_materials/ folder contains materials used to reproduce the empirical example and simulation study reported in the manuscript.

These materials are not required for typical use of the workflow, but are provided for transparency and reproducibility.

  📄 Empirical_example/ — Data and scripts for the applied example
  
  📄 simulation_study/ — Code and results for the Monte Carlo simulation

Most users can focus on Step 1–3 without using these materials.

---

## 📦 Repository Contents

* 📁 **EMAuxiliary-Workflow/**
  * 📄 `README.md`
  * 📁 `EMAuxiliary/`
    * 📄 `EMAuxiliary.R`
    * 📄 `EMAuxiliary Reference Manual v0.91.pdf`
  * 📁 `Step1/`
    * 📄 `Step1_EMA_DataPrep.html`
    * 📄 `Step1_EMA_DataPrep.Rmd`
  * 📁 `Step2/`
    * 📄 `Step2_AuxEval.html`
    * 📄 `Step2_AuxEval.Rmd`
  * 📁 `Step3/`
    * 📄 `Step3_ExamplesEMAuxiliary.html`
    * 📄 `Step3_ExamplesEMAuxiliary.Rmd`
  * 📁 `Supplemental_materials/`
    * 📁 `Empirical_example/`
      * 📄 `ema_study_data.csv`
      * 📄 `ema_study_variable_codebook.csv`
      * 📁 `Empirical_example_step2/`
        * 📄 `Step2_empirical.Rmd`
        * 📄 `Step2_empirical.html`
      * 📁 `Empirical_example_step3/`
        * 📄 `EMAux_noaux_focal_output.txt`
        * 📄 `EMAux_noaux_full_output.txt`
        * 📄 `EMAux_withaux_focal_output.txt`
        * 📄 `EMAux_withaux_full_output.txt`
        * 📄 `Step3_empirical.R`
    * 📁 `simulation_study/`
      * 📄 `01_simulate_complete_dataset.R`
      * 📄 `02_Monte_Carlo_Simulation.R`
      * 📄 `blimp_diagnostics_details.csv`
      * 📄 `diagnostic_summary.csv`
      * 📄 `simulation_results_details.csv`
      * 📄 `simulation_summary.csv`
 
## 📎 Citation and Archival Links
GitHub repository (latest version): https://github.com/schneids111/EMAuxiliary-Workflow

## 🧠 Contact
For questions, please contact:
Stefan Schneider, PhD, University of Southern California
Email: schneids@usc.edu
