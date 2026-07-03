# pySCENIC Commands for Regulon Analysis

This file documents the non-R commands used in the regulon analysis workflow.
Run these commands after exporting `mono_d6h4_count.csv` from
`trim_regulon_scenic_pipeline.R`.

## 1. Create loom input

Run in Python, for example inside the `dynamic` conda environment.

```python
import loompy as lp
import numpy as np
import scanpy as sc

x = sc.read_csv("mono_d6h4_count.csv")

row_attrs = {"Gene": np.array(x.var_names)}
col_attrs = {"CellID": np.array(x.obs_names)}

lp.create("mono.loom", x.X.transpose(), row_attrs, col_attrs)
```

## 2. Activate pySCENIC environment

```bash
source /vol/projects/qzhan/pyscenic_env/bin/activate
```

## 3. Infer co-expression network

```bash
pyscenic grn \
  --num_workers 20 \
  --output adj.sample.tsv \
  --method grnboost2 \
  mono.loom \
  ./allTFs_hg38.txt \
  --seed 777
```

## 4. Download cisTarget resources

```bash
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather

wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl
```

## 5. Motif enrichment and regulon pruning

```bash
pyscenic ctx \
  adj.sample.tsv \
  ./hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  ./hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather \
  --annotations_fname ./motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
  --expression_mtx_fname mono.loom \
  --mode dask_multiprocessing \
  --output reg.csv \
  --num_workers 20
```

## 6. Score regulon activity per cell

```bash
pyscenic aucell \
  mono.loom \
  reg.csv \
  --output out_SCENIC.loom \
  --num_workers 20
```

After `out_SCENIC.loom` is generated, continue running the R workflow in
`trim_regulon_scenic_pipeline.R`.
