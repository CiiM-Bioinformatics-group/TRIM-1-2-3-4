#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(future)
  library(BiocParallel)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(ggplot2)
  library(reshape2)
  library(harmony)
})

# -----------------------------
# Configuration
# -----------------------------

plan("multicore", workers = 4)
set.seed(1234)

pool_dirs <- list(
  P1 = "/vol/projects/CIIM/TRIM/202402845_p1_count/outs/per_sample_outs",
  P2 = "/vol/projects/CIIM/TRIM/202402845_p2_count/outs/per_sample_outs",
  P3 = "/vol/projects/CIIM/TRIM/202402845_p3_count/outs/per_sample_outs",
  P4 = "/vol/projects/CIIM/TRIM/202402845_p4_count/outs/per_sample_outs",
  P5 = "/vol/projects/CIIM/TRIM/202402845_p5_count/outs/per_sample_outs"
)

output_dir <- "."
qc_plot_file <- file.path(output_dir, "qc.pdf")
qc_table_file <- file.path(output_dir, "qc_cell_count_summary.csv")
combined_rds_file <- file.path(output_dir, "combined_integrated_harmony_sampleNpool_annot.rds")

min_features <- 300
max_features <- 4000
max_percent_mt <- 10

doublet_workers <- 3
doublet_seed <- 1234

n_variable_features <- 3000
n_pcs <- 50
harmony_dims <- 1:30
cluster_resolution <- 0.5

# -----------------------------
# Helper functions
# -----------------------------

parse_sample_info <- function(sample_name) {
  donor <- sub(".*_Donor", "", sample_name)
  treatment <- sub("_Donor.*", "", sample_name)

  list(
    donor = donor,
    treatment = treatment
  )
}

extract_timepoint <- function(sample_name) {
  ifelse(
    grepl("^Day[0-9]+", sample_name),
    sub("^(Day[0-9]+).*", "\\1", sample_name),
    NA_character_
  )
}

read_filter_and_remove_doublets <- function(pool, pool_path, sample) {
  message("Reading: ", sample, " from ", pool)

  matrix_path <- file.path(pool_path, sample, "count", "sample_filtered_feature_bc_matrix")
  if (!dir.exists(matrix_path)) {
    warning("Skipping missing matrix path: ", matrix_path)
    return(NULL)
  }

  counts <- Read10X(matrix_path)
  obj <- CreateSeuratObject(counts = counts, project = sample, min.cells = 3)
  n_total <- ncol(obj)

  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(
    obj,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_percent_mt
  )
  n_after_filter <- ncol(obj)

  sample_info <- parse_sample_info(sample)
  obj$sample <- sample
  obj$pool <- pool
  obj$donor <- sample_info$donor
  obj$treatment <- sample_info$treatment
  obj$timepoint <- extract_timepoint(sample)

  bp <- MulticoreParam(doublet_workers, RNGseed = doublet_seed)
  sce <- as.SingleCellExperiment(obj, assay = "RNA")
  sce <- scDblFinder(sce, samples = "sample", BPPARAM = bp)

  obj$scDblFinder.class <- sce$scDblFinder.class
  obj$scDblFinder.score <- sce$scDblFinder.score

  n_after_doublet <- sum(obj$scDblFinder.class == "singlet")
  obj <- subset(obj, subset = scDblFinder.class == "singlet")

  qc_row <- data.frame(
    sample = sample,
    pool = pool,
    donor = sample_info$donor,
    treatment = sample_info$treatment,
    n_total = n_total,
    n_after_filter = n_after_filter,
    n_after_doublet = n_after_doublet,
    stringsAsFactors = FALSE
  )

  list(object = obj, qc = qc_row)
}

assign_celltypes <- function(obj) {
  obj$celltype <- "undefined"
  obj$celltype[obj$seurat_clusters %in% c("0", "1")] <- "CD4_T"
  obj$celltype[obj$seurat_clusters %in% c("2", "6")] <- "CD8_T"
  obj$celltype[obj$seurat_clusters %in% c("3")] <- "NK"
  obj$celltype[obj$seurat_clusters %in% c("4")] <- "B"
  obj$celltype[obj$seurat_clusters %in% c("5", "11", "14")] <- "Mono/Macro"
  obj$celltype[obj$seurat_clusters %in% c("12")] <- "Platelet"
  obj$celltype[obj$seurat_clusters %in% c("13")] <- "DC"
  obj$celltype[obj$seurat_clusters %in% c("9")] <- "Proliferating"
  obj$celltype[obj$seurat_clusters %in% c("7")] <- "Mast"
  obj$celltype[obj$seurat_clusters %in% c("10")] <- "Plasma"
  obj
}

# -----------------------------
# Step 1: Load samples and QC
# -----------------------------

seurat_list <- list()
qc_stats <- data.frame()

for (pool in names(pool_dirs)) {
  pool_path <- pool_dirs[[pool]]
  sample_names <- list.dirs(pool_path, recursive = FALSE, full.names = FALSE)

  for (sample in sample_names) {
    result <- read_filter_and_remove_doublets(pool, pool_path, sample)
    if (is.null(result)) {
      next
    }

    key <- paste(pool, sample, sep = "_")
    seurat_list[[key]] <- result$object
    qc_stats <- rbind(qc_stats, result$qc)
  }
}

if (length(seurat_list) == 0) {
  stop("No valid samples were loaded. Please check pool_dirs and matrix paths.")
}

# -----------------------------
# Step 2: Save QC summary
# -----------------------------

qc_stats_long <- melt(
  qc_stats,
  id.vars = c("sample", "pool", "donor", "treatment"),
  measure.vars = c("n_total", "n_after_filter", "n_after_doublet"),
  variable.name = "Stage",
  value.name = "Cell_Count"
)

qc_plot <- ggplot(qc_stats_long, aes(x = sample, y = Cell_Count, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(
    title = "Cell Count at Each Filtering Stage",
    y = "Number of Cells",
    x = "Sample"
  )

ggsave(filename = qc_plot_file, plot = qc_plot, width = 14, height = 4)
write.csv(qc_stats, file = qc_table_file, row.names = FALSE)

# -----------------------------
# Step 3: Merge and integrate by Harmony
# -----------------------------

if (length(seurat_list) == 1) {
  combined <- seurat_list[[1]]
} else {
  combined <- merge(
    seurat_list[[1]],
    y = seurat_list[-1],
    add.cell.ids = names(seurat_list)
  )
}

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, nfeatures = n_variable_features)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = n_pcs)

combined <- RunHarmony(
  object = combined,
  group.by.vars = "treatment",
  reduction = "pca",
  dims.use = harmony_dims
)

combined <- RunUMAP(combined, reduction = "harmony", dims = harmony_dims)
combined <- FindNeighbors(combined, reduction = "harmony", dims = harmony_dims)
combined <- FindClusters(combined, resolution = cluster_resolution)

# -----------------------------
# Step 4: Annotate and save
# -----------------------------

combined <- assign_celltypes(combined)
saveRDS(combined, file = combined_rds_file)

message("Done.")
message("QC plot: ", qc_plot_file)
message("QC table: ", qc_table_file)
message("Combined object: ", combined_rds_file)
