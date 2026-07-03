#integrate ATAC
.libPaths(c("/vol/projects/zzhang/tools/R/4.2.0/library/","/vol/projects/qzhan/miniconda/lib/R/library","/vol/projects/qzhan/miniconda3/envs/dynamic/lib/R/library/","/vol/projects/qzhan/miniconda3/lib/R/library"))
library(future)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(SeuratData)

plan("multiprocess", workers = 10)
options(future.globals.maxSize = 80 * 1024^10)
combined = readRDS("./merged_unprocessedObj.rds")

combined_filtered = subset(combined,subset = nCount_ATAC < 17000 & nCount_ATAC > 1500 & TSS.enrichment < 10 & TSS.enrichment > 3 & nucleosome_signal < 3)

DefaultAssay(combined_filtered) = "ATAC"
combined_filtered = FindTopFeatures(combined_filtered, min.cutoff = 50)
combined_filtered =  RunTFIDF(combined_filtered)
combined_filtered <- RunSVD(combined_filtered)
#change batch effect 
seur_list = SplitObject(combined_filtered,split.by = "sample")
print("Finding integration anchors")
integration.anchors <- FindIntegrationAnchors(object.list = seur_list,anchor.features = rownames(seur_list[['202402854b_count']]),reduction = "rlsi",dims = 2:30)
print("Integrate Embeddings")
integrated_atac = IntegrateEmbeddings(anchorset = integration.anchors,reductions = combined_filtered[["lsi"]], new.reduction.name = "integrated_lsi",dims.to.integrate=1:30)
print("RunUMAP, find neighbours and clusters")
integrated_atac = RunUMAP(integrated_atac, reduction = "integrated_lsi", dims = 2:30)
integrated_atac = FindNeighbors(integrated_atac,reduction = "integrated_lsi",dims = 2:30)
integrated_atac = FindClusters(integrated_atac,verbose = FALSE, algorithm = 3,resolution = 0.2)
#Idents(integrated_atac) = "ATAC_snn_res.0.2"

print("finished integration, now gene activity assay creation")
gene.activities <- GeneActivity(integrated_atac)
integrated_atac[["RNA"]]  = CreateAssayObject(counts = gene.activities)
integrated_atac <- NormalizeData(
  object = integrated_atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(integrated_atac$nCount_RNA)
)
saveRDS(integrated_atac,file="ATACIntegrated.rds")
