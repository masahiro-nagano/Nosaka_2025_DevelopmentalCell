suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(rtracklayer))
suppressMessages(library(R.utils))
library(RColorBrewer)
library(ggh4x)
library(pals)
library(viridisLite)
library(ggrepel)

colorscale <- scale_color_gradientn(
  colors = rev(brewer.pal(11, "RdYlBu")), 
  values = rescale(c(-1, 0, 1), to = c(0, 1)) 
)

marks <- c('K4me3', 'K27me3', 'K36me3')

marks %>%
  lapply(., function(m){
    list.files('./source/' ,pattern = '_pool_mm10.10kb_fpm.bedGraph') %>%
      setNames(. ,sub('_pool_mm10.10kb_fpm.bedGraph', '', .)) -> paths
    
    paths %>%
      mapply(function(x, n){
        fread(sprintf('./source/%s', x), col.names = c('chr', 'start', 'end', 'count')) %>%
          mutate(count = log2((count+1)/(end - start)*1000) +1) %>%
          `colnames<-`(c('chr', 'start', 'end', n))
      }, ., names(.), SIMPLIFY = F) %>%
      Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .)
  }) %>% Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) -> epigenome.dat


list.files('./source' ,pattern = '_methyl_10000b.bedGraph', full.names = T) %>%
  setNames(. ,sub('.*\\/(.*)_methyl_10000b.bedGraph', '\\1', .))  -> methyl.paths

methyl.paths %>%
  mapply(function(x, n){
    fread(x) %>%
      `colnames<-`(c('chr', 'start', 'end', n))
  }, ., names(.), SIMPLIFY = F) %>%
  Reduce(function(...) left_join(..., by =c(chr="chr", start="start", end = "end"), all.x = T), .) %>%
  mutate(start = start -1) -> methyl.data

colorscale <- scale_color_gradientn(
  colors = rev(brewer.pal(11, "RdYlBu")), 
  limits = c(-4, 4), # actual limits
  oob = scales::squish, # values beyond limit/out-of-bound are "squished" to the limit
  breaks = c(-4, -2, 0, 2, 4)
)


epigenome.dat %>%
  left_join(., methyl.data, by =c(chr="chr", start="start", end = "end")) %>%
  na.omit() %>%
  mutate(d_K4me3 = c34_H18_K4me3 - GV_house_K4me3, d_K36me3 = c34_H18_K36me3 - GV_house_K36me3, d_K27me3 = c34_H18_K27me3 - GV_house_K27me3, d_methyl = d4c34 - FGO) %>%
  ggplot(aes(x = d_methyl, y = d_K36me3, col = d_K4me3)) +
  geom_point(alpha=.4) +
  geom_density2d(col = 'red') +
  colorscale +
  theme_bw()-> p
p

ggsave(p, filename = "sf5_i.pdf", width = 8.11 , height = 4.61)
