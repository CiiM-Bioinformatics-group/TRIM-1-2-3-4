# TRIM-1-2-3-4

This repository contains custom R, Python and command-line scripts used for single-cell RNA-seq, single-nucleus ATAC-seq and downstream integrative analyses, including differential gene expression, chromatin accessibility, regulon inference, cell-cell communication and disease enrichment analyses.

## Overview

Peripheral blood mononuclear cells (PBMCs) from healthy donors were trained in vitro and profiled using single-cell RNA sequencing and single-nucleus ATAC sequencing. The analysis workflow includes:

- scRNA-seq quality control, doublet removal, normalization, integration and cell-type annotation.
- Differential gene expression analysis across trained immunity conditions.
- snATAC-seq quality control and integration.
- Gene activity scoring and transfer of cell-type labels from the annotated scRNA-seq reference.
- Regulon and transcription factor activity analysis using pySCENIC.
- Cell-cell communication analysis using CellChat.
- Gene set enrichment analyses in external disease datasets.

## Repository Contents

| File | Description |
| --- | --- |
| `trim_qc_harmony_pipeline.R` | scRNA-seq sample loading, QC filtering, scDblFinder doublet removal, Harmony integration, clustering and broad cell-type annotation. |
| `DEG_analysis.R` | Differential gene expression analysis for Day 6, Day 0 + 4 h and Day 6 + 4 h comparisons across major immune cell types. |
| `TIGs_accross4types.R` | Heatmap for comparison the effect size of TIGs across 4 TRIMs. Enrichment for the module genes. |
| `cell_cell_interaction.R` | CellChat analysis comparing trained conditions with matched RPMI controls. |
| `Regulon_scenic_pipeline.R` | R workflow for regulon activity analysis, differential regulon activity testing and TIG-associated TF prioritization. |
| `ATAC_step1.R`, `ATAC_step2.R`, `ATAC_step3.R` | Steps for snATAC-seq data pre-processiong and integration. |
| `dapeak.R`,`daTF.R` | Differential accessibility peaks and differential acitivity TFs. |
| `Regulon_scenic_commands.md` | Python and shell commands used to create loom files and run pySCENIC. |
| `trim_external_disease_gsea_pipeline.R` | Gene set enrichment analysis of TRIM-induced gene signatures in external disease datasets. |

## Main Software and Packages

The analyses were performed using established open-source software and R/Python packages, including:

- **Single-cell RNA-seq:** Cell Ranger, Seurat, scDblFinder, Harmony, Azimuth.
- **Single-nucleus ATAC-seq:** cellranger-arc, Cell Ranger ATAC, Signac, Seurat, GenomicRanges, EnsDb.Hsapiens.v86, Harmony, MACS2.
- **Demultiplexing and genotype matching:** Souporcell and BCFtools.
- **Differential analysis and visualization:** Seurat, ggplot2, dplyr and related tidyverse packages.
- **Pathway and disease enrichment:** clusterProfiler, MSigDB Hallmark gene sets and Reactome gene sets.
- **Regulon analysis:** pySCENIC, GRNBoost2, cisTarget databases, AUCell and SCopeLoomR.
- **Cell-cell communication:** CellChat.

## Data Availability

The scRNA-seq and snATAC-seq data generated in this study will be made available through the European Genome-phenome Archive (EGA). Accession codes are pending and will be provided upon publication. 


## Reproducibility

Custom scripts used for data preprocessing, statistical analysis and visualization are provided in this repository.

# sessionInfo()
JASPAR2020_0.99.10 
BiocParallel_1.30.3 
scDblFinder_1.10.0
ggsignif_0.6.4                   
harmony_0.1.1                    
Rcpp_1.0.10                      
glmGamPoi_1.8.0                   
DoubletFinder_2.0.3              
readxl_1.4.0                      
circlize_0.4.15                  
ComplexHeatmap_2.12.1             
ReactomePA_1.40.0                
TFBSTools_1.34.0                  
DOSE_3.22.0                      
org.Hs.eg.db_3.15.0               
clusterProfiler_4.4.4            
EnhancedVolcano_1.14.0            
BSgenome.Hsapiens.UCSC.hg38_1.4.4
BSgenome_1.64.0                   
rtracklayer_1.56.1               
Biostrings_2.64.0                 
XVector_0.36.0                   
EnsDb.Hsapiens.v86_2.99.0         
ensembldb_2.20.2                 
AnnotationFilter_1.20.0           
GenomicFeatures_1.48.3           
AnnotationDbi_1.58.0              
forcats_0.5.1                    
stringr_1.4.1                     
purrr_1.0.2                      
readr_2.1.2                       
tidyr_1.2.0                      
tibble_3.2.1                      
tidyverse_1.3.2                  
ggrepel_0.9.1 
harmony_0.1.1
writexl_1.4.1
Signac_1.10.0                     
SingleCellExperiment_1.18.1      
SummarizedExperiment_1.26.1       
Biobase_2.56.0                   
GenomicRanges_1.48.0              
GenomeInfoDb_1.35.15             
IRanges_2.30.0                    
S4Vectors_0.34.0                 
BiocGenerics_0.42.0               
MatrixGenerics_1.8.1             
matrixStats_0.62.0                
ggplot2_3.4.0                    
dplyr_1.1.4                       
sp_1.5-0                         
SeuratObject_4.1.0                
Seurat_4.1.0 


## Citation

If you use this code, please cite the associated manuscript once available.
