# Ketamine-Depressive-Circuits-Analysis
This repository contains the R-based analytical pipeline for Mendelian Randomization (MR) and genomic visualization for the study: "Ketamine rewires depressive circuits via coordinated multi-cell signaling".
## 1.Overview
This project integrates human brain cell-type-specific eQTL data with mouse single-nucleus transcriptomics to identify key causal genes involved in major depressive disorder (MDD) and their rescue by ketamine treatment. By evaluating 8 key cell subtypes, we identified 8 core causal genes as central nodes in ketamine's rapid antidepressant effect.
## 2. Computational Workflow & Script Details
### Step 0: Genomic Coordinate Alignment (0.eqtl_tob37.r)
Function: Converts human eQTL association data from hg38 to the b37 (GRCh37) system to match MDD GWAS statistics.

Logic: Implements "sign flipping" for beta values to maintain effect direction relative to the dbSNP (v153) reference.
