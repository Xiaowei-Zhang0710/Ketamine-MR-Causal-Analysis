# Ketamine-Depressive-Circuits-Analysis
This repository contains the R-based analytical pipeline for Mendelian Randomization (MR) and genomic visualization for the study: "Ketamine rewires depressive circuits via coordinated multi-cell signaling".
## 1.Overview
This project integrates human brain cell-type-specific eQTL data with mouse single-nucleus transcriptomics to identify key causal genes involved in major depressive disorder (MDD) and their rescue by ketamine treatment. By evaluating 8 key cell subtypes, we identified 8 core causal genes as central nodes in ketamine's rapid antidepressant effect.
## 2. Computational Workflow & Script Details
### Step 0: Genomic Coordinate Alignment (0.eqtl_tob37.r)
Function: Converts human eQTL association data from hg38 to the b37 (GRCh37) system to match MDD GWAS statistics.  
Logic: Implements "sign flipping" for beta values to maintain effect direction relative to the dbSNP (v153) reference.
### Step 1: Causal Inference & Local Mapping (1.run.smr_all_process.r)
Two-Sample MR: Uses the TwoSampleMR package ($P \le 10^{-5}$, $r^2 < 0.01$) to evaluate the causality of mouse DEGs in human MDD risk.  
Gviz Plotting: Visualizes specific risk loci  relative to gene structures using biomaRt.
### Step 2: Result Aggregation and Leading SNP Identification (2.plot0.get_cytoband_full.r) 
Purpose: To serve as a bridge between MR statistical results and spatial visualization.   
Core Logic: Aggregates multi-omics causal results across cell types. It identifies a leading SNP for each causal target gene.   
Priority Logic Logic: Identifies the most representative leading SNP by prioritizing the SNP with the lowest $P$-value from the GWAS results; if none is available, defaults to the SNP with the lowest $P$-value from the exposure eQTLs.   
Output: Generates a standardized coordinate table for genome-wide mapping.
### Step 3: Genomic Landscape Visualization (2.plot1_cytoband_with_degSTAR.r)
Figure 3q: Displays the complete genomic distribution of the 8 causal genes across human chromosomes 1–22.  
Visual Mapping: Integrates cytoband ideograms (using hg19_cytoBand.txt), SNP markers , cell-type associations (dots), and ketamine-rescued DEG  (stars) .
## 3. Requirements
Language: R (v4.3.2)  
Key Packages: TwoSampleMR (v0.5.6), Seurat (v3.1.4), Gviz.
## 4. Reference Data
Cytobands: hg19_cytoBand.txt (included in this repository) provides the chromosomal structure for hg19/b37 assembly.
## 5. Note on Spatial Analysis
The Python-based spatial deconvolution and mapping pipeline (Tangram) is available in our companion repository: [插入第二个仓库链接].
