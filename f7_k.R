library(tidyverse)
library(data.table)
library(pheatmap)

# fread('/Users/masahiro/Desktop/collaboration/Nosaka/Nosaka_data/SC3_RNA_expression_all.tsv') %>%
#   dplyr::rename(gene_name = 'V1') -> expression_all
# 
# 
# 
# # For counting export defined regions
# fread("/Volumes/MasaHD1/local/NGS_data_analysis/reference_genome/GRCm38p3_wo_pseudo_v1b_chrx_gene.bed", select = c(1:4, 5), col.names = c('chr','start','end', 'strand', 'name')) -> genes
# 
# 
# marks <- c('K4me3', 'K27me3', 'K36me3', 'Input')
# 
# 
# list.files('/Users/masahiro/Desktop/collaboration/Nosaka/counts/genes' ,pattern = '_pool_mm10.genes_fpm.bedGraph') %>%
#   setNames(. ,sub('_pool_mm10.genes_fpm.bedGraph', '', .)) -> paths
# 
# paths %>%
#   mapply(function(x, n){
#     fread(sprintf('/Users/masahiro/Desktop/collaboration/Nosaka/counts/genes/%s', x), col.names = c('chr', 'start', 'end', 'count')) %>%
#       mutate(count = log2((count+1)/(end - start)*1000) +1) %>%
#       `colnames<-`(c('chr', 'start', 'end', n))
#   }, ., names(.), SIMPLIFY = F) %>%
#   Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>%
#   left_join(genes, .) -> epigenome.data
# 
# list.files('/Volumes/MasaHD1/local/NGS_data_analysis/Methylome/emseq/bgs' ,pattern = '_methyl_Gene.bedGraph', full.names = T) %>%
#   setNames(. ,sub('.*\\/(.*)_Gene.bedGraph', '\\1', .))  -> methyl.paths
# 
# methyl.paths %>%
#   mapply(function(x, n){
#     fread(x) %>%
#       `colnames<-`(c('chr', 'start', 'end', n, 'Name'))
#   }, ., names(.), SIMPLIFY = F) %>%
#   Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end", Name = 'Name'), all.x = T), .) %>%
#   left_join(genes, ., by = c(chr = "chr", start = "start", end = "end", name = 'Name')) -> methyl.data
# 
# ########################################################################################################################################################
# ########################################################################################################################################################
# # all gene with RNA
# expression_all %>%
#   # left_join(., methyl.data, by = (c(gene_name ='name'))) %>%
#   left_join(., epigenome.data, by = (c(gene_name ='name'))) %>%
#   na.omit() -> all.epi.rna.dat
# 
# save(file = '/Volumes/MasaHD1/local/NGS_data_analysis2/github/Nosaka_2025_DevelopmentalCell/source/all.epi.rna.dat.rda', all.epi.rna.dat)

load('./source/all.epi.rna.dat.rda') 
# without gene expression
all.epi.rna.dat %>% 
  dplyr::select(c("GV_house_K36me3", "GV_house_K27me3", "GV_house_K4me3")) %>%
  na.omit() %>%
  # Z-score
  lapply(., function(x){
    mean(x)-> m
    sd(x) -> s
    ((x-m)/s)
  }) %>%  
  data.frame()%>%
  pheatmap::pheatmap( annotation_legend = TRUE,
                      annotation_names_row = F, annotation_names_col = TRUE,
                      clustering_distance_rows = "correlation",
                      drop_levels = TRUE, show_rownames = F, show_colnames = T,
                      cluster_cols = T, cluster_rows = T, color=colorRampPalette(c("blue", "white", "red"))(n=299)) -> p

p$tree_row$order -> heat.order

# in each 
all.epi.rna.dat %>% 
  dplyr::select(c("Antral_FGO_vivo_RNA","GV_house_K36me3", "GV_house_K27me3", "GV_house_K4me3")) %>%
  .[heat.order, ] -> resorted

lim <- 2
resorted %>%
  # Z-score
  lapply(., function(x){
    mean(x)-> m
    sd(x) -> s
    ((x-m)/s)
  }) %>%  
  data.frame() %>%
  # make more than 6 are all red
  mutate_all(function(x){
    x= case_when(x > lim ~ lim,
                 x < -lim ~ -lim,
                 TRUE ~ x)
  }) %>%
  pheatmap::pheatmap(annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = TRUE,
                     drop_levels = TRUE,
                     show_rownames = FALSE,
                     show_colnames = T,
                     cluster_cols = FALSE,
                     cluster_rows = F,
                     color = colorRampPalette(c("blue", "white", "red"))(n = 299)) -> p

pdf(file = 'f7_k1.pdf', width =4, height = 6.8)
p
dev.off()


# each cell type (c34)
all.epi.rna.dat %>% 
  dplyr::select(c("c34_H18_RNA", "c34_H18_K36me3", "c34_H18_K27me3","c34_H18_K4me3")) %>%
  .[heat.order, ] -> resorted

lim <- 2
resorted %>%
  # Z-score
  lapply(., function(x){
    mean(x)-> m
    sd(x) -> s
    ((x-m)/s)
  }) %>%  
  data.frame() %>%
  # make more than 6 are all red
  mutate_all(function(x){
    x= case_when(x > lim ~ lim,
                 x < -lim ~ -lim,
                 TRUE ~ x)
  }) %>%
  pheatmap::pheatmap(annotation_legend = TRUE,
                     annotation_names_row = FALSE,
                     annotation_names_col = TRUE,
                     drop_levels = TRUE,
                     show_rownames = FALSE,
                     show_colnames = T,
                     cluster_cols = FALSE,
                     cluster_rows = F,
                     color = colorRampPalette(c("blue", "white", "red"))(n = 299)) -> p

pdf(file = 'f7_k2.pdf', width =4, height = 6.8)
p
dev.off()

all.epi.rna.dat %>% dplyr::transmute(gv = Antral_FGO_vivo_RNA, c34 = c34_H18_RNA) -> rna

all.epi.rna.dat %>% dplyr::transmute(gv = GV_house_K4me3, c34 = c34_H18_K4me3) -> k4me3

all.epi.rna.dat %>% dplyr::transmute(gv = GV_house_K27me3, c34 = c34_H18_K27me3) -> k27me3

all.epi.rna.dat %>% dplyr::transmute(gv = GV_house_K36me3, c34 = c34_H18_K36me3) -> k36me3

bind_rows(rna, k4me3, k27me3, k36me3) %>%
  na.omit() %>% cor() -> c
