# ALS Gene Browser

> An interactive Shiny application for exploring transcriptomic data from the NYGC ALS Consortium, allowing researchers to visualize differential gene expression across tissues, genetic subtypes, and clinical variables in ALS and control post-mortem samples.

---

## Background

The Center for Genomics of Neurodegenerative Diseases (CGND) flagship ALS manuscript presents a comprehensive analysis of the full consortium cohort, focusing on **RNA-seq** and **whole genome sequencing (WGS)** data from post-mortem tissue collected from ALS patients and neurologically normal controls.

The ALS Gene Browser provides a point-and-click interface for querying gene-level results without requiring direct interaction with raw data or analysis pipelines.

---

## Dataset Overview

| Feature | Details |
|---|---|
| Data types | RNA-seq, clinical metadata |
| Sample types | Post-mortem human tissue |
| Groups | ALS cases, neurologically normal controls |
| Contributing sites | 41 sites across the NYGC ALS Consortium |
| Tissues profiled | Motor Cortex (MCX), Frontal Cortex (FCX), Cervical Spinal Cord (CSC), Lumbar Spinal Cord (LSC), Cerebellum (CBL) |

---

## Application Overview

The browser has **two main tabs**:

### Tab 1 — ALS vs. Control Differential Expression

Explore DEG results comparing **ALS vs. neurologically normal controls** across all five tissues.

- Enter a gene (HGNC symbol, e.g. `TARDBP`, `FUS`, `SOD1`) to return results across all five tissues simultaneously.
- **Boxplot:** normalized expression by group, faceted by tissue.
- **DEG table:** fold change (LFC), p-value, FDR, and direction of effect. Positive LFC = higher expression in ALS; negative LFC = higher expression in controls.
- DEG analysis uses **robust covariates** to correct for cross-site batch effects in sample handling, RNA extraction, and sequencing (see paper).

### Tab 2 — Case-Only Analysis (C9orf72 & Clinical Variables)

Explore expression **within ALS cases**, stratified by *C9orf72* repeat expansion status and correlated with clinical variables.

- Enter a gene and select a tissue.
- **C9orf72 boxplot:** C9orf72-positive vs. -negative cases in the chosen tissue.
- **Clinical scatter plots:** gene expression vs. age at death, and vs. disease duration (symptom onset to death).
- **DEG table:** results from three analyses — C9orf72 status (categorical), age at death (continuous), and disease duration (continuous).
- All analyses use **DESeq2**; Motor Cortex additionally uses **DREAM** (differential expression for repeated measures).

---

## Citation

If you use this browser or the underlying data in your research, please cite the NYGC ALS Consortium flagship manuscript:

> NYGC ALS Consortium. The New York Genome Center ALS Consortium resource integrates postmortem tissue transcriptomics and whole genome sequencing to empower biological discovery. *medRxiv* (2026). https://www.medrxiv.org/content/10.64898/2026.04.29.26350889v1

> Humphrey, J., Oku, A., Byrska-Bishop, M., Basile, A. O., Evani, U. S., Corvelo, A., Tokolyi, A., Bp, K., Réal, A., Kim, Y., Bond, M. L., Clarke, W. E., Fu, R., Geiger, H., Chang, S., Naito, T., Jang, B., Musunuri, R., Dredge, W. H., Al-Abri, R., … Phatnani, H. (2026). The New York Genome Center ALS Consortium resource integrates postmortem tissue transcriptomics and whole genome sequencing to empower biological discovery. *medRxiv: the preprint server for health sciences*, 2026.04.29.26350889. https://doi.org/10.64898/2026.04.29.26350889



For questions about the browser or underlying data, contact the **Center for Genomics of Neurodegenerative Diseases (CGND)** at the New York Genome Center.

---