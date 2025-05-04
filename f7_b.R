library(data.table)
library(tidyverse)
library(corrplot)

#K36me3
mark <- 'K36me3'

# check correlation's UHC
list.files(path = './source', pattern = '_K36me3_pool_mm10.10kb_fpm.bedGraph') -> nms
nms %>% setNames(., sub('_K36me3_pool_mm10.10kb_fpm.bedGraph', '', .)) %>%
  lapply(., function(n){
    fread(sprintf('./source/%s', n)) %>% .$V4 %>% {log2(. + .01)}
  }) %>% data.frame() %>%
  dplyr::select(-c(d4c7CsA_BDF121)) %>%
  mutate(m = apply(., 1, min)) %>%
  filter(!m == min(m)) %>%
  dplyr::select(-m) -> mat

cor(mat[, c("c13_H18", "E18_house",  "c34_H18", "GV_house")]) -> c

col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))

# K36me3
pdf('f7_b1.pdf', height = 4.8, width = 7.8)
corrplot(c, method="color", col=col(200),  
         diag=FALSE, tl.pos="d", tl.col = 'black',
         type="full", order='original', 
         addCoef.col = "black", # Add coefficient of correlation
         mar=c(0,0,1,0)
)
dev.off()


#K4me3
mark <- 'K4me3'

# check correlation's UHC
list.files(path = './source', pattern = '_K4me3_pool_mm10.10kb_fpm.bedGraph') -> nms
nms %>% setNames(., sub('_K4me3_pool_mm10.10kb_fpm.bedGraph', '', .)) %>%
  lapply(., function(n){
    fread(sprintf('./source/%s', n)) %>% .$V4 %>% {log2(. + .01)}
  }) %>% data.frame() %>%
  dplyr::select(-c(d4c7CsA_BDF121, d4c7RAB2_BDF121, GO_P7)) %>%
  mutate(m = apply(., 1, min)) %>%
  filter(!m == min(m)) %>%
  dplyr::select(-m) -> mat

cor(mat[, c("c13_H18", "E18_house",  "c34_H18", "GV_house")]) -> c

col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))

pdf('f7_b2.pdf', height = 5, width = 5)
corrplot(c, method="color", col=col(200),  
         diag=FALSE, tl.pos="d", tl.col = 'black',
         type="full", order='original', 
         # title=title, 
         addCoef.col = "black", # Add coefficient of correlation
         mar=c(0,0,1,0) # http://stackoverflow.com/a/14754408/54964
)
dev.off()


#K27me3
mark <- 'K27me3'

# check correlation's UHC
list.files(path = './source', pattern = '_K27me3_pool_mm10.10kb_fpm.bedGraph') -> nms
nms %>% setNames(., sub('_K27me3_pool_mm10.10kb_fpm.bedGraph', '', .)) %>%
  lapply(., function(n){
    fread(sprintf('./source/%s', n)) %>% .$V4 %>% {log2(. + .01)}
  }) %>% data.frame() %>%
  dplyr::select(-c(d4c7CsA_BDF121, GO_P7)) %>%
  mutate(m = apply(., 1, min)) %>%
  filter(!m == min(m)) %>%
  dplyr::select(-m) -> mat

cor(mat[, c("c13_H18", "E18_house",  "c34_H18", "GV_house")]) -> c

col <- colorRampPalette(rev(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA")))

pdf('f7_b3.pdf', height = 5, width = 5)
corrplot(c, method="color", col=col(200),  
         diag=FALSE, tl.pos="d", tl.col = 'black',
         type="full", order='original', 
         # title=title, 
         addCoef.col = "black", # Add coefficient of correlation
         mar=c(0,0,1,0) # http://stackoverflow.com/a/14754408/54964
)
dev.off()
