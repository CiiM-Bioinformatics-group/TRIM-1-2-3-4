#!/usr/bin/env Rscript

# TRIM regulon analysis workflow
#
# This script contains the R parts of the regulon analysis:
#   1. Export Mono/Macro Day6_4h count matrix and metadata for pySCENIC
#   2. Load pySCENIC regulon AUC results
#   3. Test differential regulon activity across stimulation conditions
#   4. Identify TIG-associated TFs from GRNBoost2 links
#   5. Plot top regulons and condition-specific TFs
#
# The Python and shell commands required to run pySCENIC are documented in:
#   trim_regulon_scenic_commands.md

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(SCopeLoomR)
  library(AUCell)
  library(PRROC)
  library(parallel)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

options(stringsAsFactors = FALSE)
set.seed(42)

# -----------------------------
# Input files
# -----------------------------

combined_rds <- "../combined_integrated_harmony_sampleNpool_annot.rds"
scenic_loom <- "out_SCENIC.loom"
adj_file <- "adj.sample.tsv"
deg_file <- "/vol/projects/qzhan/TRIM/RNA/all_de_day6_4h.csv"

# -----------------------------
# Output files
# -----------------------------

mono_count_file <- "mono_d6h4_count.csv"
mono_meta_file <- "mono_meta.csv"
da_regulon_rds <- "da_regulons_paired.rds"
top_regulon_pdf <- "top10_da_regulons_paired.pdf"
netea_rds <- "netea_TF_results.rds"
shared_specific_file <- "netea_shared_specific_TFs.txt"
tf_heatmap_pdf <- "TIG_TF_binary_heatmap.pdf"

# -----------------------------
# Step 1: Export Mono/Macro data for pySCENIC
# -----------------------------

combined <- readRDS(combined_rds)

submono <- subset(
  combined,
  subset = day_group == "Day6_4h" & celltype == "Mono/Macro"
)

write.csv(
  t(as.matrix(submono@assays$RNA@counts)),
  file = mono_count_file
)

write.csv(
  submono@meta.data,
  file = mono_meta_file
)

# After this step, run pySCENIC outside R.
# See trim_regulon_scenic_commands.md for the Python and shell commands.

# -----------------------------
# Step 2: Load pySCENIC regulon activity
# -----------------------------

loom <- open_loom(scenic_loom)
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
close_loom(loom)

submono$group <- paste0(submono$stim1, "_", submono$stim2)

cellInfo <- data.frame(
  condition = submono$group,
  row.names = colnames(submono)
)

auc_mat <- getAUC(regulonAUC)
auc_mat <- auc_mat[, rownames(cellInfo)]
stopifnot(identical(colnames(auc_mat), rownames(cellInfo)))

pairs <- list(
  Asp = list(target = "Asp_Candida", ref = "RPMI_Candida"),
  BCG = list(target = "BCG_Staph", ref = "RPMI_Staph"),
  IVIG = list(target = "IVIG_LPS", ref = "RPMI_LPS"),
  R848 = list(target = "R848_PolyIC", ref = "RPMI_PolyIC")
)

# -----------------------------
# Step 3: Differential regulon activity
# -----------------------------

test_paired <- function(auc_mat, cellInfo, target_cond, ref_cond) {
  cells_target <- rownames(cellInfo)[cellInfo$condition == target_cond]
  cells_ref <- rownames(cellInfo)[cellInfo$condition == ref_cond]

  cat(
    "  ", target_cond, " (n=", length(cells_target), ") vs ",
    ref_cond, " (n=", length(cells_ref), ")\n",
    sep = ""
  )

  results <- apply(auc_mat, 1, function(x) {
    v_target <- x[cells_target]
    v_ref <- x[cells_ref]

    if (all(v_target == 0) && all(v_ref == 0)) {
      return(c(mean_target = 0, mean_ref = 0, delta = 0, pvalue = 1))
    }

    wt <- suppressWarnings(wilcox.test(v_target, v_ref))

    c(
      mean_target = mean(v_target),
      mean_ref = mean(v_ref),
      delta = mean(v_target) - mean(v_ref),
      pvalue = wt$p.value
    )
  })

  results <- as.data.frame(t(results))
  results$regulon <- rownames(results)
  results$padj <- p.adjust(results$pvalue, method = "BH")
  results <- results[order(-abs(results$delta)), ]

  results
}

da_regulons <- list()

for (training in names(pairs)) {
  cat("\n======== ", training, " training ========\n", sep = "")

  da_regulons[[training]] <- test_paired(
    auc_mat = auc_mat,
    cellInfo = cellInfo,
    target_cond = pairs[[training]]$target,
    ref_cond = pairs[[training]]$ref
  )

  n_up <- sum(
    da_regulons[[training]]$padj < 0.05 &
      da_regulons[[training]]$delta > 0
  )
  n_down <- sum(
    da_regulons[[training]]$padj < 0.05 &
      da_regulons[[training]]$delta < 0
  )

  cat("  Sig UP: ", n_up, " | Sig DOWN: ", n_down, "\n", sep = "")
}

saveRDS(da_regulons, da_regulon_rds)

# -----------------------------
# Step 4: Plot top regulons
# -----------------------------

top_per_pair <- list()

for (training in names(da_regulons)) {
  df <- da_regulons[[training]]
  sig <- df[df$padj < 0.05, ]

  top_up <- head(
    sig[sig$delta > 0, ][order(-sig[sig$delta > 0, ]$delta), ],
    10
  )
  top_down <- head(
    sig[sig$delta < 0, ][order(sig[sig$delta < 0, ]$delta), ],
    10
  )

  top_per_pair[[training]] <- list(up = top_up, down = top_down)

  cat("\n=== ", training, " TOP 10 UP ===\n", sep = "")
  print(top_up[, c("regulon", "delta", "padj", "mean_target", "mean_ref")])

  cat("\n=== ", training, " TOP 10 DOWN ===\n", sep = "")
  print(top_down[, c("regulon", "delta", "padj", "mean_target", "mean_ref")])
}

top_union <- unique(unlist(lapply(top_per_pair, function(x) {
  c(x$up$regulon, x$down$regulon)
})))

cat("\nTotal unique regulons in panel: ", length(top_union), "\n", sep = "")

delta_mat <- matrix(
  0,
  nrow = length(top_union),
  ncol = length(pairs),
  dimnames = list(top_union, names(pairs))
)

padj_mat <- matrix(
  1,
  nrow = length(top_union),
  ncol = length(pairs),
  dimnames = dimnames(delta_mat)
)

for (training in names(pairs)) {
  df <- da_regulons[[training]]
  idx <- match(top_union, df$regulon)
  delta_mat[, training] <- df$delta[idx]
  padj_mat[, training] <- df$padj[idx]
}

sig_label <- matrix(
  "",
  nrow = nrow(padj_mat),
  ncol = ncol(padj_mat),
  dimnames = dimnames(padj_mat)
)

sig_label[padj_mat < 0.05] <- "*"
sig_label[padj_mat < 0.01] <- "**"
sig_label[padj_mat < 0.001] <- "***"

desired_order <- c("BCG", "Asp", "R848", "IVIG")
delta_mat <- delta_mat[, desired_order]
sig_label <- sig_label[, desired_order]

rng <- max(abs(delta_mat), na.rm = TRUE)

pdf(top_regulon_pdf, width = 5, height = max(8, length(top_union) * 0.22))
pheatmap(
  delta_mat,
  color = colorRampPalette(c("#3a3a98", "white", "#cc0000"))(50),
  breaks = seq(-rng, rng, length.out = 51),
  display_numbers = sig_label,
  fontsize_number = 8,
  number_color = "black",
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  clustering_method = "ward.D2",
  main = "Trained immunity regulon programs\n(AUC: trained vs RPMI control)",
  fontsize_row = 8,
  border_color = NA
)
dev.off()

for (training in names(da_regulons)) {
  sig <- da_regulons[[training]][da_regulons[[training]]$padj < 0.05, ]
  sig <- sig[order(-abs(sig$delta)), ]
  sig$training <- training
  sig$target_cond <- pairs[[training]]$target
  sig$ref_cond <- pairs[[training]]$ref

  write.csv(
    sig[, c(
      "regulon", "training", "target_cond", "ref_cond",
      "delta", "mean_target", "mean_ref", "pvalue", "padj"
    )],
    file = paste0("DA_regulons_", training, "_paired.csv"),
    row.names = FALSE
  )
}

cat("\nDifferential regulon outputs:\n")
cat("  - ", da_regulon_rds, "\n", sep = "")
cat("  - DA_regulons_<training>_paired.csv\n")
cat("  - ", top_regulon_pdf, "\n", sep = "")

# -----------------------------
# Step 5: TIG-TF enrichment
# -----------------------------

adj <- read.table(adj_file, header = TRUE, sep = "\t")

cat("Total TF-target pairs: ", nrow(adj), "\n", sep = "")
cat("Unique TFs: ", length(unique(adj$TF)), "\n", sep = "")
cat("Unique targets: ", length(unique(adj$target)), "\n", sep = "")

D6h4 <- read.csv(deg_file)
D6h4_fil <- D6h4 %>%
  filter(celltype == "Mono/Macro")

D6h4_fil$group[D6h4_fil$group == "Asper"] <- "Asp"

deg_list <- list()
for (cond in c("Asp", "BCG", "IVIG", "R848")) {
  sub <- D6h4_fil[D6h4_fil$group == cond, ]
  sub <- sub[!duplicated(sub$gene), ]
  deg_list[[cond]] <- sub
}

tig_list <- list()
for (cond in names(deg_list)) {
  df <- deg_list[[cond]]
  tig <- df$gene[df$p_val_adj < 0.05]
  tig_list[[cond]] <- tig
  cat(cond, " TIG genes: ", length(tig), "\n", sep = "")
}

calc_tf_score <- function(tf_links, tig_set, n_perm = 500, seed = 42) {
  labels <- as.integer(tf_links$target %in% tig_set)
  scores <- tf_links$importance

  n_pos <- sum(labels)
  n_total <- length(labels)

  if (n_pos < 3 || n_pos == n_total) {
    return(list(
      n_target_in_TIG = n_pos,
      AUC = NA,
      AUPR = NA,
      AUC_p = NA,
      AUPR_p = NA
    ))
  }

  roc_real <- roc.curve(
    scores.class0 = scores[labels == 1],
    scores.class1 = scores[labels == 0],
    curve = FALSE
  )

  pr_real <- pr.curve(
    scores.class0 = scores[labels == 1],
    scores.class1 = scores[labels == 0],
    curve = FALSE
  )

  real_auc <- roc_real$auc
  real_aupr <- pr_real$auc.integral

  set.seed(seed)
  perm_aucs <- numeric(n_perm)
  perm_auprs <- numeric(n_perm)

  for (i in seq_len(n_perm)) {
    perm_labels <- sample(labels)

    perm_aucs[i] <- tryCatch(
      roc.curve(
        scores.class0 = scores[perm_labels == 1],
        scores.class1 = scores[perm_labels == 0],
        curve = FALSE
      )$auc,
      error = function(e) NA
    )

    perm_auprs[i] <- tryCatch(
      pr.curve(
        scores.class0 = scores[perm_labels == 1],
        scores.class1 = scores[perm_labels == 0],
        curve = FALSE
      )$auc.integral,
      error = function(e) NA
    )
  }

  auc_p <- (sum(perm_aucs >= real_auc, na.rm = TRUE) + 1) / (n_perm + 1)
  aupr_p <- (sum(perm_auprs >= real_aupr, na.rm = TRUE) + 1) / (n_perm + 1)

  list(
    n_target_in_TIG = n_pos,
    AUC = real_auc,
    AUPR = real_aupr,
    AUC_p = auc_p,
    AUPR_p = aupr_p
  )
}

run_netea_for_condition <- function(adj, tig_set, n_perm = 500, n_cores = 10) {
  all_TFs <- unique(adj$TF)
  cat("Testing ", length(all_TFs), " TFs against gene set of size ", length(tig_set), "\n", sep = "")

  cat("Pre-filtering TFs with >=3 targets in TIG...\n")
  tfs_to_test <- character()

  for (tf in all_TFs) {
    tf_targets <- adj$target[adj$TF == tf]
    if (sum(tf_targets %in% tig_set) >= 3) {
      tfs_to_test <- c(tfs_to_test, tf)
    }
  }

  cat("TFs after pre-filter: ", length(tfs_to_test), "\n", sep = "")
  cat("Running with ", n_cores, " cores...\n", sep = "")

  results <- mclapply(tfs_to_test, function(tf) {
    tf_links <- adj[adj$TF == tf, ]
    tf_links <- tf_links[order(-tf_links$importance), ]

    res <- calc_tf_score(tf_links, tig_set, n_perm = n_perm)
    res$TF <- tf
    res$total_targets <- nrow(tf_links)
    res
  }, mc.cores = n_cores)

  df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      TF = x$TF,
      total_targets = x$total_targets,
      n_target_in_TIG = x$n_target_in_TIG,
      AUC = x$AUC,
      AUPR = x$AUPR,
      AUC_p = x$AUC_p,
      AUPR_p = x$AUPR_p
    )
  }))

  df[order(df$AUC_p, df$AUPR_p), ]
}

netea_results <- list()

for (cond in c("Asp", "BCG", "IVIG", "R848")) {
  cat("\n=========== Running ", cond, " ===========\n", sep = "")

  netea_results[[cond]] <- run_netea_for_condition(
    adj = adj,
    tig_set = tig_list[[cond]],
    n_perm = 500,
    n_cores = 10
  )

  sig <- netea_results[[cond]][
    netea_results[[cond]]$AUC_p < 0.05 &
      netea_results[[cond]]$AUPR_p < 0.05,
  ]

  cat("Significant TFs: ", nrow(sig), "\n", sep = "")

  if (nrow(sig) > 0) {
    sig <- sig[order(sig$AUC_p), ]
    print(head(sig, 15))
  }
}

saveRDS(netea_results, netea_rds)

# -----------------------------
# Step 6: Shared and condition-specific TFs
# -----------------------------

sig_tfs_per_condition <- lapply(netea_results, function(df) {
  df$TF[df$AUC_p < 0.05 & df$AUPR_p < 0.05]
})

cat("\n=== Significant TFs per condition ===\n")
print(sapply(sig_tfs_per_condition, length))

master_tfs <- Reduce(intersect, sig_tfs_per_condition)
cat("\n=== Master TFs shared across all 4 conditions ===\n")
print(master_tfs)

specific_tfs <- list()
for (cond in names(sig_tfs_per_condition)) {
  others_union <- Reduce(
    union,
    sig_tfs_per_condition[setdiff(names(sig_tfs_per_condition), cond)]
  )
  specific_tfs[[cond]] <- setdiff(sig_tfs_per_condition[[cond]], others_union)
}

cat("\n=== Condition-specific TFs ===\n")
for (cond in names(specific_tfs)) {
  cat(
    cond, " (", length(specific_tfs[[cond]]), " specific): ",
    paste(head(specific_tfs[[cond]], 20), collapse = ", "),
    "\n",
    sep = ""
  )
}

for (cond in names(netea_results)) {
  sig <- netea_results[[cond]][
    netea_results[[cond]]$AUC_p < 0.05 &
      netea_results[[cond]]$AUPR_p < 0.05,
  ]

  sig <- sig[order(sig$AUC_p), ]
  write.csv(sig, paste0("netea_sigTFs_", cond, ".csv"), row.names = FALSE)
}

sink(shared_specific_file)
cat("=== Master TFs shared across all 4 conditions ===\n")
cat(paste(master_tfs, collapse = ", "), "\n\n")

for (cond in names(specific_tfs)) {
  cat("=== ", cond, " specific (n=", length(specific_tfs[[cond]]), ") ===\n", sep = "")
  cat(paste(specific_tfs[[cond]], collapse = ", "), "\n\n")
}
sink()

# -----------------------------
# Step 7: Binary TF heatmap
# -----------------------------

all_tfs <- unique(unlist(sig_tfs_per_condition))
cat("Total unique significant TFs: ", length(all_tfs), "\n", sep = "")

binary_mat <- sapply(sig_tfs_per_condition, function(tfs) {
  all_tfs %in% tfs
})

rownames(binary_mat) <- all_tfs
mode(binary_mat) <- "integer"

binary_mat <- binary_mat[, desired_order]
n_sig <- rowSums(binary_mat)
binary_mat <- binary_mat[order(-n_sig), ]

col_fun <- c("0" = "white", "1" = "#dcaa8a")

pdf(tf_heatmap_pdf, width = 3, height = 5)
Heatmap(
  binary_mat,
  name = "Identified",
  col = col_fun,
  cluster_columns = FALSE,
  cluster_rows = FALSE,
  show_column_names = TRUE,
  column_names_side = "top",
  column_names_rot = 0,
  column_names_centered = TRUE,
  column_names_gp = gpar(fontsize = 10, fontface = "bold"),
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  row_names_side = "right",
  border = TRUE,
  rect_gp = gpar(col = "grey80", lwd = 0.3),
  heatmap_legend_param = list(
    title = "TF identified",
    at = c(0, 1),
    labels = c("No", "Yes"),
    title_gp = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
)
dev.off()

cat("\nDone. Outputs:\n")
cat("  - ", da_regulon_rds, "\n", sep = "")
cat("  - DA_regulons_<training>_paired.csv\n")
cat("  - ", top_regulon_pdf, "\n", sep = "")
cat("  - ", netea_rds, "\n", sep = "")
cat("  - netea_sigTFs_<condition>.csv\n")
cat("  - ", shared_specific_file, "\n", sep = "")
cat("  - ", tf_heatmap_pdf, "\n", sep = "")
