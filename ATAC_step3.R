#Integrate ATAC and RNA
library(future)
library(Seurat)
library(Signac)
library(SeuratDisk)
library(SeuratData)


pbmc_rna <- readRDS("../../combined_integrated_harmony_sampleNpool_annot.rds")
pbmc_atac <- readRDS("./ATACIntegrated.rds")
#pbmc_atac <- readRDS("./ATACIntegrated_harmony.rds")
DefaultAssay(pbmc_atac) <- 'RNA'
pbmc_rna <- UpdateSeuratObject(pbmc_rna)
transfer.anchors <- FindTransferAnchors(
  reference = pbmc_rna,
  query = pbmc_atac,
  reduction = 'cca'
)

predicted.labels <- TransferData(
  anchorset = transfer.anchors,
  refdata = pbmc_rna$celltype,
  weight.reduction = pbmc_atac[['lsi']],
  dims = 2:30
)

pbmc_atac <- AddMetaData(object = pbmc_atac, metadata = predicted.labels)

saveRDS(pbmc_atac,file='./ATAC_annot.rds')

plot1 <- DimPlot(
  object = pbmc_rna,
  group.by = 'celltype',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scRNA-seq')

plot2 <- DimPlot(
  object = pbmc_atac,
  group.by = 'predicted.id',
  label = TRUE,
  repel = TRUE) + NoLegend() + ggtitle('scATAC-seq')


pdf('anno_ATAC.pdf',width=10,height=5)
plot1 + plot2
dev.off()
