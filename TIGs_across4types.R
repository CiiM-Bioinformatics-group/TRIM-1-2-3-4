#heatmap, use bcg
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(colorRamp2)
D6h4_fil=D6h4 %>% filter(celltype=='Mono/Macro' )
df_heat <- D6h4_fil %>%
  mutate(fc_signif = ifelse(p_val < 0.05, avg_log2FC, 0)) %>%
  group_by(gene) %>%
  mutate(any_sig = any(p_val_adj < 0.05)) %>%
  ungroup() %>%
  filter(any_sig) %>%
  select(gene, comp, fc_signif) %>% as.data.frame()
#head(df_heat)
df_heat_wide <- df_heat %>%
  pivot_wider(
    names_from = comp,       
    values_from = fc_signif, 
    values_fill = 0         
  )
df_heat_wide=data.frame(df_heat_wide)

mat <- as.matrix(df_heat_wide[, -1])         
rownames(mat) <- df_heat_wide$gene 
colnames(mat) <- c("R848", "BCG","Asp","IVIG")
mat <- mat[, c("BCG", "Asp", "R848", "IVIG")]
col<- list(
  Type = c('R848'='#623A8A','BCG'='#99212C','Asp'='#428331',"IVIG"="#4774BA")
)
la <- HeatmapAnnotation(
  Type = c("BCG", "Asp","R848","IVIG"),
  col=col,simple_anno_size=unit(2,"mm"),
    annotation_legend_param = list(
    title_gp  = gpar(fontsize = 8, fontface = "bold"),
    labels_gp = gpar(fontsize =8)
))
set.seed(123)
optimal_k=9
pa=cluster::pam(mat, k = optimal_k)
clustering=data.frame(pa$clustering)
clustering$gene=rownames(clustering)
module_counts <- table(clustering$pa.clustering)
module_info <- paste0("module ", names(module_counts), " (n=", module_counts, ")")
clustering$pa.clustering <- module_info[match(clustering$pa.clustering, names(module_counts))]
row_indices <- match(clustering$gene, rownames(mat))
clustering$index <- row_indices
module_gene_counts <- table(pa$clustering)
genelist <- c('IFI30','CLEC7A','CD274','CCL8','CXCL16','IDO1','STAT1','GBP1',
              'ISG15','OAS1','IFIT1','IFIT2','IFI44','IFIH1','CXCL10','CXCL11','CD163','GCH1',
              'IL1R1','IRAK1','CD1D','CEBPD','NFE2L2','SOD2','IL10RA','P2RX7','IL6','IL1B','TNF',
             'NLRP3','PSTPIP1','C5AR1','MMP7','SPP1','CCL22','CHI3L1','FN1','IL18','S100A4','NFKBIA','NFKBID','MAP3K11','CARD9','TREM2','SIGLEC9','CEBPA'
             )
index <- which(rownames(mat) %in% genelist)
col_fun <- colorRamp2(
  breaks = c(-2, 0, 2),        
  colors = c("#3F1B5A", "white", "#D55E00")  
)
labs <- rownames(mat)[index]
lab2=rowAnnotation(foo=anno_mark(at=index,labels=labs,labels_gp=gpar(fontsize=5),lines_gp = gpar()))

heatmap_legend <- Legend(
   title = "Effect size",
   col_fun = col_fun,
   at = c(-2, 0, 2),
   labels = c("-2", "0", "2"),title_gp  = gpar(fontsize = 8, fontface = "bold"),
  labels_gp = gpar(fontsize = 8)
 )
pdf('heatmap_TIGs_Orig_4Type.pdf',width=4,height=5)
ht=Heatmap(mat,
        name=NULL,
        cluster_columns = F,show_row_dend = T,col = col_fun,
        right_annotation = lab2,
        cluster_rows = T,top_annotation = la, 
        border = TRUE,clustering_method_columns='single',
        show_row_names = F,show_column_names = FALSE,
        row_split = clustering$pa.clustering,row_title_rot = 0,
        row_title_gp = gpar(col=rep('black',optimal_k),fontsize=rep(8,optimal_k)),
        row_names_gp = gpar(fontsize = 8),column_dend_height = unit(0.1, "cm"),show_heatmap_legend = FALSE)
draw(ht, annotation_legend_list = list(heatmap_legend))
dev.off()
write.csv(clustering,file='cluster_Orig_4Type.csv')


###ORA
m_df <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)


run_enrichments <- function(genes, module_name) {
  # ID
  m.df <- bitr(genes, fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Hs.eg.db)
  
  # Reactome
  reactome_res <- enrichPathway(gene = m.df$ENTREZID,
                                pvalueCutoff = 0.05,
                                readable = TRUE)
  reactome_df <- as.data.frame(reactome_res)
  if (nrow(reactome_df) > 0) {
    reactome_df$Source <- "Reactome"
  }
  
  # Hallmark
  hallmark_res <- enricher(gene = genes,
                           TERM2GENE = m_df,
                           pvalueCutoff = 0.05)
  hallmark_df <- as.data.frame(hallmark_res)
  if (nrow(hallmark_df) > 0) {
    hallmark_df$Source <- "Hallmark"
  }
  
  # 
  enrich_all <- bind_rows(reactome_df, hallmark_df)
  if (nrow(enrich_all) > 0) {
    enrich_all <- enrich_all %>%
      filter(p.adjust < 0.05) %>%
      mutate(qscore = -log10(p.adjust),
             Module = module_name)
  }
  return(enrich_all)
}

clustering=read.csv('cluster_Orig_4Type.csv')
modules <- unique(clustering$pa.clustering)
results_list <- list()

#modules <- c("module 3 (n=84)")


for (m in modules) {
  cat("Processing:", m, "\n")
  genes_m <- subset(clustering, pa.clustering == m)$gene
  enrich_m <- run_enrichments(genes_m, m)
  if (nrow(enrich_m) > 0) {
    results_list[[m]] <- enrich_m
  }
}

all_results <- bind_rows(results_list)

top_df <- all_results %>%
  group_by(Module, Source) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  ungroup()

top_df[top_df$Module=='module 3 (n=84)',]$Module='3'
top_df[top_df$Module=='module 2 (n=96)',]$Module='2'
top_df[top_df$Module=='module 6 (n=133)',]$Module='6'
top_df[top_df$Module=='module 5 (n=162)',]$Module='5'
top_df[top_df$Module=='module 8 (n=82)',]$Module='8'
top_df[top_df$Module=='module 1 (n=41)',]$Module='1'
top_df[top_df$Module=='module 4 (n=27)',]$Module='4'

top_df=top_df[top_df$Module %in% c('3','2','6','5','8','1','4'),]


pdf("Modules_ORA_Top3.pdf", width = 7, height = 8)

ggplot(top_df, aes(
  x = qscore,
  y = fct_reorder(Description, qscore),
  fill = qscore          # 
)) +
  geom_col() +
  facet_grid(
    rows = vars(Module),
    scales = "free_y",
    space = "free_y"
  ) +
  scale_fill_gradient(
    low = "#FCEEE8",      # 0 
    high = '#FAAE9E'         ###'#d3968c'    
  ) +
  scale_x_continuous(position = "bottom") +
  theme_bw(base_size = 13) +
  labs(
    x = expression(-log[10]("p.adjust")),
    y = NULL,
    fill = expression(-log[10]("p.adjust"))
  ) +
  theme(
    axis.text.y = element_text(size = 10, hjust = 1),
    axis.text.x.bottom = element_text(size = 10),
    axis.title.x.bottom = element_text(size = 12, vjust = -1),
    legend.position = "top",
    strip.placement = "outside",strip.background = element_rect(fill="white"),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 10)
  )

dev.off()
