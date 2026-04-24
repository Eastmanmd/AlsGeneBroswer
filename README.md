# ALS Gene Browser

> An interactive Shiny application for exploring transcriptomic data from the NYGC ALS Consortium allowing researchers to visualize differential gene expression across tissues, genetic subtypes, and clinical variables in ALS and control post-mortem samples.

---

## Background

The Center of Genomics of Neurodegenerative Diseases (CGND) flagship ALS manuscripts present a comprehensive analysis of the full consortium cohort, focusing on **RNA-seq** and **whole genome sequencing (WGS)** data generated from post-mortem tissue samples collected from ALS patients and neurologically normal controls. 

The ALS Gene Browser provides an accessible, point-and-click interface for researchers to query gene-level results without requiring direct interaction with raw data or analysis pipelines.

---

## Dataset Overview

| Feature | Details |
|---|---|
| Data types | RNA-seq,  clinical metadata|
| Sample types | Post-mortem human tissue |
| Groups | ALS cases, Neurologically normal controls |
| Contributing sites | 41 sites across the NYGC ALS Consortium |
| Tissues profiled | 5 major CNS and non-CNS tissues (see below) |

### Tissues

1. Motor Cortex (MCX)
2. Frontal Cortex (FCX)
3. Cervical Spinal Cord (CSC)
4. Lumbar Spinal Cord (LSC)
5. Cerebellum (CBL)

---

## Application Overview

The browser is organized into **two main tabs**

---

## Tab 1 — ALS vs. Control Differential Expression

### Purpose
Explore differential gene expression (DEG) results comparing gene expression data from **ALS samples to neurologically normal controls** across all five tissues.

### How to Use

1. **Enter a gene of interest** in the gene input field (HGNC symbol, e.g., `TARDBP`, `FUS`, `SOD1`).
2. The app will return results across **all five tissues simultaneously**.

### Outputs

#### Boxplot
- Displays normalized expression of the selected gene across ALS and control samples, faceted by tissue.

#### DEG Results Table
- Summary of differential expression statistics for the queried gene in each tissue.
- Includes fold change (LFC), p-value, adjusted p-value (FDR), and direction of effect.
- DEG includes covariate selection (see paper)
- A positive LFC shows higher gene expression in ALS cases compared to controls.
- A negative LFC shows higher genne expression in controls compared to ALS cases. 

### Statistical Approach
- DEG analysis was performed using **robust covariates** to account for sample collection across multiple contributing sites. (see paper)
- Covariate correction mitigates batch effects introduced by cross-site variability in sample handling, RNA extraction, and sequencing. (see paper)

---

## Tab 2 — Case-Only Analysis (C9orf72 & Clinical Variables)

### Purpose
Explore gene expression patterns **within ALS cases **, stratified by *C9orf72* repeat expansion status and correlated with key clinical variables including age at symptom onset and disease duration. 

### How to Use

1. **Enter a gene of interest** in the gene input field.
2. **Select a tissue** from the five available tissues using the dropdown menu.

### Outputs

#### C9orf72 Boxplot
- Visualizes expression of the selected gene comparing **C9orf72-positive vs. C9orf72-negative** ALS cases in the chosen tissue.

#### Clinical Variable Scatter Plots
- **Age at Death vs. Gene Expression** — scatter plot showing the relationship between expression and age at death across ALS cases.
- **Disease Duration vs. Gene Expression** — scatter plot showing the relationship between expression and disease duration (time from symptom onset to death).

#### DEG Results Table
- Results from three complementary analyses:
  - Summary of differential expression statistics for the queried gene in each tissue.
  - Includes fold change, p-value, adjusted p-value (FDR), and direction of effect.
  - **C9orf72 status** (categorical): DEG between C9orf72-positive and C9orf72-negative cases.

#### Note:
  - **Age at Death** (continuous): Association between gene expression and age at death.
  - **Disease Duration** (continuous): Association between gene expression and disease duration.

### Statistical Approach
- All analyses use **DESeq2** for count-based differential expression modeling.
- For **Motor Cortex only**, analyses were additionally run using **DREAM** (differential expression for repeated measures).

---

## Repository Structure

```
ALS_Gene_Browser/
├── app.R                  # Main Shiny application (or ui.R + server.R)
├── data/
│   ├── deg_als_vs_ctrl.rds   # DEG results: ALS vs. Control (Tab 1)
│   ├── deg_case_only_results.rds     # DEG results: C9orf72 pos vs. neg (Tab 2)
│   └── tpm_data.rds # Normalized expression across all samples
└── README.md
```

---

## Citation

If you use this browser or the underlying data in your research, please cite the NYGC ALS Consortium flagship manuscript:

> *[Full citation to be added upon publication]*

For questions about the browser or the underlying data, please contact the **Center of Genomics of Neurodegenerative Diseases (CGND)** at the New York Genome Center.

---
