library(Seurat)
library(dplyr)
library(ggplot2)
library(SingleCellExperiment)
library(Signac) 
library(writexl)
library(ggrepel)
library(tidyverse)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v86)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnhancedVolcano)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(JASPAR2020)
library(TFBSTools)
library(ReactomePA)
library(ComplexHeatmap)
library(circlize)
library(readxl)
library(DoubletFinder)
library(glmGamPoi)
library(harmony)  
library(DoubletFinder)
library(ggsignif)
library(scDblFinder)
library(BiocParallel)
options(stringsAsFactors = FALSE)
rm(list=ls())


pool_dirs <- list(
  P1 = "./TRIM/atac_count/202402849_p1_count/outs",
  P2 = "./TRIM/atac_count/202402849_p2_count/outs",
  P3 = "./TRIM/atac_count/202402849_p3_count/outs",
  P4 = "./TRIM/atac_count/202402849_p4_count/outs",
  P5 = "./TRIM/atac_count/202402849_p5_count/outs",
  pilot = "./TRIM/202402854b_count/outs/outs"
)

batches=c('202402849_p1_count','202402849_p2_count','202402849_p3_count','202402849_p4_count','202402849_p5_count','202402854b_count')

gr_list = list()
for (i in batches) {
  if(i != '202402854b_count'){
   peaks_bed = paste0("./TRIM/atac_count/",i,"/outs/peaks.bed") 
  }
    else{
    peaks_bed = paste0("./TRIM/202402854b_count/outs/peaks.bed") 
    }
  peaksFile = read.table(file = peaks_bed,col.names = c("chr", "start", "end"))
  grFile = makeGRangesFromDataFrame(peaksFile)
  gr_list[[i]] = grFile
}

combined.peaks = GenomicRanges::reduce(x = c(gr_list[['202402849_p1_count']],
                              gr_list[['202402849_p2_count']],
                              gr_list[['202402849_p3_count']],
                              gr_list[['202402849_p4_count']],
                              gr_list[['202402849_p5_count']],
                              gr_list[['202402854b_count']]))


peakwidths = width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
seur_list = list()

for (i in batches) {

  print(paste0('Processing pool: ', i))

  # Reading in singlets pre-processed file#
  singlets <- read.csv(paste0('./TRIM/demultiplex/', i, '/singlet_meta.csv'))
  rownames(singlets) = singlets$barcode
  # Reading in 10X data ATAC
    if ( i != '202402854b_count') {
        data <- Read10X_h5(paste0("./TRIM/atac_count/", i, "/outs/filtered_peak_bc_matrix.h5"))
        }
    else{
        data <- Read10X_h5('./TRIM/202402854b_count/outs/filtered_peak_bc_matrix.h5')
    }
      
  atac = data[,colnames(data) %in% singlets$barcode]

  ########creating atac assay #####
    if ( i == '202402854b_count') {
        frags <- CreateFragmentObject(paste0("./TRIM/202402854b_count/outs/fragments.tsv.gz"),cells = colnames(atac))
        }
    else {
        frags <- CreateFragmentObject(paste0("./TRIM/atac_count/", i,"/outs/fragments.tsv.gz"),cells = colnames(atac))
    }
  
  counts <- FeatureMatrix(fragments = frags,features = combined.peaks,cells = colnames(atac))
  chromAssay <- CreateChromatinAssay(counts, fragments = frags, min.cells=10)
  
              
 seur <- CreateSeuratObject(
    counts = chromAssay,
    assay = "ATAC",
    meta.data = singlets
    )
  seur$poolID=i
  Annotation(seur) = annotations
  seur = NucleosomeSignal(seur)
  seur = TSSEnrichment(seur)
  seur_list[[i]] <- seur
}

mergedSeurObj = merge(x=seur_list[['202402854b_count']],y=c(
                                             seur_list[['202402849_p1_count']],
                                             seur_list[['202402849_p2_count']],
                                             seur_list[['202402849_p3_count']],
                                             seur_list[['202402849_p4_count']],
                                             seur_list[['202402849_p5_count']]
                                          ))

saveRDS(mergedSeurObj,file = "./merged_unprocessedObj.rds")
