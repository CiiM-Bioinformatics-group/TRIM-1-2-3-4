#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

options(stringsAsFactors = FALSE)


#deg,day6_4h
combined=readRDS('../combined_integrated_harmony_sampleNpool_annot.rds')
celltypes=c('B','CD4_T','CD8_T','Mono/Macro','NK')

all_de=NULL
all_de_sig=NULL
for(cell in celltypes){
    submono=subset(combined,celltype==cell & stim2 == 'LPS')
    de.uhD <- FindMarkers(submono, group.by='stim1', ident.1='IVIG', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhD$gene=rownames(de.uhD)
    de.uhD$celltype=cell
    de.uhD$day='Day6_4h'
    de.uhD$comp='Day6IVIGLPSvsDay6RPMILPS'
    de.uhD$group='IVIG'
    de.uhD.sig = de.uhD %>% filter(p_val_adj < 0.05)
    
     
    submono=subset(combined,celltype==cell & stim2=='Candida')
    de.uhC <- FindMarkers(submono, group.by='stim1', ident.1='Asp', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhC$gene=rownames(de.uhC)
    de.uhC$celltype=cell
    de.uhC$group='Asper'
    de.uhC$day='Day6_4h'
    de.uhC$comp='Day6AspCandidavsDay6RPMICandida'
    de.uhC.sig = de.uhC %>% filter(p_val_adj < 0.05)

    submono=subset(combined,celltype==cell & stim2=='Staph')
    de.uhB <- FindMarkers(submono, group.by='stim1', ident.1='BCG', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhB$gene=rownames(de.uhB)
    de.uhB$celltype=cell
    de.uhB$day='Day6_4h'
    de.uhB$group='BCG'
    de.uhB$comp='Day6BCGStaphvsDay6RPMIStaph'
    de.uhB.sig = de.uhB %>% filter(p_val_adj < 0.05)

    submono=subset(combined,celltype==cell & stim2=='PolyIC')
    de.uhA <- FindMarkers(submono, group.by='stim1', ident.1='R848', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhA$gene=rownames(de.uhA)
    de.uhA$celltype=cell
    de.uhA$day='Day6_4h'
    de.uhA$group='R848'
    de.uhA$comp='Day6R848PolyICvsDay6RPMIPolyIC'
    de.uhA.sig = de.uhA %>% filter(p_val_adj < 0.05)
    
    all_de=rbind(all_de,de.uhA,de.uhB,de.uhC,de.uhD)
    all_de_sig=rbind(all_de_sig,de.uhA.sig,de.uhB.sig,de.uhC.sig,de.uhD.sig)
}
write.csv(all_de,file='./all_de_day6_4h.csv')
write.csv(all_de_sig,file='./all_de_sig_day6_4h.csv')

#day6
Day6=subset(rna,day_group=='Day6')
all_de=NULL
all_de_sig=NULL
for(cell in celltypes){
    sub=subset(Day6,celltype==cell)
    de.uhD <- FindMarkers(sub, group.by='stim1', ident.1='IVIG', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhD$gene=rownames(de.uhD)
    de.uhD$celltype=cell
    de.uhD$day='Day6'
    de.uhD$comp='IVIGvsRPMI'
    de.uhD$group='IVIG'
    de.uhD.sig = de.uhD %>% dplyr::filter(p_val_adj < 0.05)

    de.uhC <- FindMarkers(sub, group.by='stim1', ident.1='Asp', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhC$gene=rownames(de.uhC)
    de.uhC$celltype=cell
    de.uhC$group='Asper'
    de.uhC$day='Day6'
    de.uhC$comp='AspvsRPMI'
    de.uhC.sig = de.uhC %>% dplyr::filter(p_val_adj < 0.05)

    de.uhB <- FindMarkers(sub, group.by='stim1', ident.1='BCG', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhB$gene=rownames(de.uhB)
    de.uhB$celltype=cell
    de.uhB$day='Day6'
    de.uhB$group='BCG'
    de.uhB$comp='BCGvsRPMI'
    de.uhB.sig = de.uhB %>% dplyr::filter(p_val_adj < 0.05)

    de.uhA <- FindMarkers(sub, group.by='stim1', ident.1='R848', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhA$gene=rownames(de.uhA)
    de.uhA$celltype=cell
    de.uhA$day='Day6'
    de.uhA$group='R848'
    de.uhA$comp='R848vsRPMI'
    de.uhA.sig = de.uhA %>% dplyr::filter(p_val_adj < 0.05)
    print(paste0('finished',cell))
    all_de=rbind(all_de,de.uhA,de.uhB,de.uhC,de.uhD)
    all_de_sig=rbind(all_de_sig,de.uhA.sig,de.uhB.sig,de.uhC.sig,de.uhD.sig)
}
write.csv(all_de,file='./all_de_day6.csv')
write.csv(all_de_sig,file='./all_de_sig_day6.csv')

#Day0_4h,deg
D0h4=subset(rna,day_group=='Day0_4h')
all_de=NULL
all_de_sig=NULL
for(cell in celltypes){
    sub=subset(D0h4,celltype==cell)
    de.uhD <- FindMarkers(sub, group.by='stim1', ident.1='IVIG', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhD$gene=rownames(de.uhD)
    de.uhD$celltype=cell
    de.uhD$day='Day0_4h'
    de.uhD$comp='IVIGvsRPMI'
    de.uhD$group='IVIG'
    de.uhD.sig = de.uhD %>% dplyr::filter(p_val_adj < 0.05)
    
    de.uhC <- FindMarkers(sub, group.by='stim1', ident.1='Asp', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhC$gene=rownames(de.uhC)
    de.uhC$celltype=cell
    de.uhC$group='Asper'
    de.uhC$day='Day0_4h'
    de.uhC$comp='AspvsRPMI'
    de.uhC.sig = de.uhC %>% dplyr::filter(p_val_adj < 0.05)

    de.uhB <- FindMarkers(sub, group.by='stim1', ident.1='BCG', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhB$gene=rownames(de.uhB)
    de.uhB$celltype=cell
    de.uhB$day='Day0_4h'
    de.uhB$group='BCG'
    de.uhB$comp='BCGvsRPMI'
    de.uhB.sig = de.uhB %>% dplyr::filter(p_val_adj < 0.05)

    de.uhA <- FindMarkers(sub, group.by='stim1', ident.1='R848', ident.2='RPMI', min.pct=0.1, logfc.threshold=0)
    de.uhA$gene=rownames(de.uhA)
    de.uhA$celltype=cell
    de.uhA$day='Day0_4h'
    de.uhA$group='R848'
    de.uhA$comp='R848vsRPMI'
    de.uhA.sig = de.uhA %>% dplyr::filter(p_val_adj < 0.05)
    print(paste0('finished',cell))
    all_de=rbind(all_de,de.uhA,de.uhB,de.uhC,de.uhD)
    all_de_sig=rbind(all_de_sig,de.uhA.sig,de.uhB.sig,de.uhC.sig,de.uhD.sig)
}
write.csv(all_de,file='./all_de_day0h4.csv')
write.csv(all_de_sig,file='./all_de_sig_day0h4.csv')

##use tig perform pca
TImarker=D6h4[D6h4$celltype=='Mono/Macro' & D6h4$p_val_adj<0.05,]$gene
TI=subset(combined,day_group=='Day6_4h' & celltype=='Mono/Macro' ) 
TI$cond <- paste(TI$stim1, TI$stim2, sep = "_")
pseudo <- AggregateExpression(
  TI,
  assays = "RNA",
  slot = "counts",   # counts
  group.by = c("donor", "cond"),
  return.seurat = TRUE
)
expr <- GetAssayData(pseudo, slot = "counts")
expr <- expr[intersect(TImarker, rownames(expr)), ]

logcounts <- log1p(t(t(expr) / colSums(expr)) * 1e6)
mat <- t(logcounts)
vars <- apply(mat, 2, var)
mat2 <- mat[, vars > 0]
pca <- prcomp(mat2, scale. = TRUE)
df_pca <- as.data.frame(pca$x) 
meta <- pseudo@meta.data
df_pca$sample <- rownames(df_pca)
df_pca$cond <- sub("^[^_]+_", "", df_pca$sample)
type_colors <- c('R848_PolyIC'='#623A8A','BCG_Staph'='#99212C','Asp_Candida'='#428331',"IVIG_LPS"="#4774BA")
df_pca_fil=df_pca[df_pca$cond %in% c('R848_PolyIC','BCG_Staph','Asp_Candida',"IVIG_LPS"),]
pct_var <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
pc1_label <- paste0("PC1 (", pct_var[1], "%)")
pc2_label <- paste0("PC2 (", pct_var[2], "%)")
p <- ggplot(df_pca_fil, aes(PC1, PC2, color = cond)) +
  geom_point(size = 3) +stat_ellipse() +
  theme_classic()+scale_color_manual(values = type_colors)+labs(x = pc1_label, y = pc2_label)
ggsave(p,file='TIG_pca_donor_PC1PC2.pdf',height=3,width=4)    

#top loading
loadings_pc <- pca$rotation[, "PC2"]

top_pos <- sort(loadings_pc[loadings_pc > 0], decreasing = TRUE)[1:10]
top_neg <- sort(loadings_pc[loadings_pc < 0])[1:10]

df_loadings <- data.frame(
  gene = c(names(top_pos), names(top_neg)),
  loading = c(top_pos, top_neg)
) %>%
  arrange(loading) %>%
  mutate(gene = factor(gene, levels = gene))

p_loading <- ggplot(df_loadings, aes(x = loading, y = gene, fill = loading > 0)) +
  geom_col() +
  scale_fill_manual(values = c("TRUE" = "grey70", "FALSE" = "grey70"), guide = "none") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_classic() +
  labs(x = "PC2 Loading", y = NULL)

ggsave(p_loading, file = 'TIG_pca_PC2_loadings.pdf', height = 3, width = 3)
