library(tidyverse)
library(data.table)
library(rtracklayer)
library(R.utils)
library(RColorBrewer)
library(ggh4x)
library(pals)
library(viridisLite)
library(ggrepel)

colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

dss <- c(
  d4c13 = "./source/d4c13_methyl_10000b.bedGraph",
  d4c34= "./source/d4c34_methyl_10000b.bedGraph",
  NG = "./source/NG_methyl_10000b.bedGraph",
  FGO = "./source/FGO_methyl_10000b.bedGraph")

samps <- c( "d4c13", "d4c34", "NG",  "FGO")

dss %>% lapply(., function(x){
  fread(x, col.names = c('chr', 'start', 'end', 'mC')) 
}) %>% 
  Reduce(function(...) left_join(..., by=c(chr= 'chr', start = 'start', end = 'end'), all.x=TRUE), .) %>%
  `colnames<-`(c('chr', 'start', 'end', names(dss))) %>%
  mutate(start = start -1) %>% na.omit() -> mc.dat

round(cor(mc.dat$d4c13, mc.dat$NG), 2) -> d4c13_NG_cor

mc.dat %>%
  ggplot(aes(x= d4c13*100, y = NG*100))+
  geom_hex(binwidth = c(1, 1)) + 
  colorscale +
  coord_fixed() +
  geom_abline(slope = 1, col="red") +
  annotate(geom="text", label= sprintf("R=%s", d4c13_NG_cor), x=80, y=95) +
  theme_bw() +
  labs(x = 'd4c13', y = 'NG')-> p
p
ggsave('f6_c1.pdf', p, height = 4, width = 5.2, device = cairo_pdf, bg = 'transparent')


round(cor(mc.dat$d4c34, mc.dat$FGO), 2) -> d4c34_FGO_cor

mc.dat %>%
  ggplot(aes(x= d4c34*100, y = FGO*100))+
  geom_hex(binwidth = c(1, 1)) + 
  colorscale +
  coord_fixed() +
  geom_abline(slope = 1, col="red") +
  annotate(geom="text", label= sprintf("R=%s", d4c34_FGO_cor), x=80, y=95) +
  theme_bw() +
  labs(x = 'd4c34', y = 'FGO')-> p
p
ggsave('f6_c2.pdf', p, height = 4, width = 5.2, device = cairo_pdf, bg = 'transparent')

