# TRIM-1-2-3-4

Code accompanying the study of four trained immunity programs induced by BCG, *Aspergillus*, R848 and IVIG in human peripheral blood immune cells.

This repository contains custom R, Python and command-line scripts used for single-cell RNA-seq, single-nucleus ATAC-seq and downstream integrative analyses, including differential gene expression, chromatin accessibility, regulon inference, cell-cell communication and disease enrichment analyses.

## Overview

Peripheral blood mononuclear cells (PBMCs) from healthy donors were trained in vitro and profiled using single-cell RNA sequencing and single-nucleus ATAC sequencing. The analysis workflow includes:

- scRNA-seq quality control, doublet removal, normalization, integration and cell-type annotation.
- Differential gene expression analysis across trained immunity conditions.
- snATAC-seq quality control, TF-IDF normalization, LSI dimensional reduction and Harmony integration.
- Gene activity scoring and transfer of cell-type labels from the annotated scRNA-seq reference.
- Peak-level chromatin accessibility analysis and peak-gene linkage analysis.
- Regulon and transcription factor activity analysis using pySCENIC.
- Cell-cell communication analysis using CellChat.
- Gene set enrichment analyses in external disease datasets.

## Repository Contents

| File | Description |
| --- | --- |
| `trim_qc_harmony_pipeline.R` | scRNA-seq sample loading, QC filtering, scDblFinder doublet removal, Harmony integration, clustering and broad cell-type annotation. |
| `trim_deg_pipeline.R` | Differential gene expression analysis for Day 6, Day 0 + 4 h and Day 6 + 4 h comparisons across major immune cell types. |
| `trim_cellchat_pipeline.R` | CellChat analysis comparing trained conditions with matched RPMI controls. |
| `Regulon_scenic_pipeline.R` | R workflow for regulon activity analysis, differential regulon activity testing and TIG-associated TF prioritization. |
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

The scRNA-seq and snATAC-seq data generated in this study will be made available through the European Genome-phenome Archive (EGA). Accession codes are pending and will be provided upon publication. Access to controlled human sequencing data will be managed through EGA in accordance with applicable ethical and data protection regulations.

Large sequencing files, processed Seurat objects and controlled-access human genomic data are not included in this repository.

## Usage

The scripts are intended to be run on a high-performance computing environment with access to the relevant raw sequencing outputs and processed intermediate objects.

Before running the scripts, update the input paths at the top of each file to match your local directory structure. Several paths in the scripts point to project-specific server locations and are provided to document the analysis workflow.

Example:

```r
combined_rds <- "../combined_integrated_harmony_sampleNpool_annot.rds"
```

The typical analysis order is:

1. Process scRNA-seq data and generate the integrated annotated object:

```bash
Rscript trim_qc_harmony_pipeline.R
```

2. Run differential gene expression analyses:

```bash
Rscript trim_deg_pipeline.R
```

3. Run CellChat communication analysis:

```bash
Rscript trim_cellchat_pipeline.R
```

4. Run pySCENIC commands described in:

```text
Regulon_scenic_commands.md
```

5. Analyze regulon activity and TF enrichment:

```bash
Rscript Regulon_scenic_pipeline.R
```

6. Run external disease gene set enrichment analyses:

```bash
Rscript trim_external_disease_gsea_pipeline.R
```

## Notes on snATAC-seq Integration

The snATAC-seq workflow used TF-IDF normalization followed by SVD/LSI dimensional reduction. The first LSI component was excluded because of its correlation with sequencing depth. Batch effects across samples and sequencing pools were corrected using Harmony on LSI dimensions 2-30. UMAP, nearest-neighbor graph construction and clustering were performed using the Harmony embeddings. Gene activity scores were computed with `GeneActivity` and log-normalized before label transfer from the annotated scRNA-seq reference.

## Reproducibility

Custom scripts used for data preprocessing, statistical analysis and visualization are provided in this repository. Analyses rely on external data files and software environments that should be configured according to the Methods section of the manuscript.

For reproducibility, users should record package versions with:

```r
sessionInfo()
```

## Citation

If you use this code, please cite the associated manuscript once available.
