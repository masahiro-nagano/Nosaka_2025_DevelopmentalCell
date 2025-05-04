suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(rtracklayer))
suppressMessages(library(R.utils))
library(RColorBrewer)
library(ggh4x)
library(pals)
library(viridisLite)
library(ggrepel)

colorscale = scale_fill_gradientn(
  colors = rev(brewer.pal(9, "YlGnBu")),
  values = c(0, exp(seq(-5, 0, length.out = 100))))

dss <- c(
  d4c7 = "./source/d4c7_methyl_10000b.bedGraph",
  d4c13 = "./source/d4c13_methyl_10000b.bedGraph",
  d4c24 = "./source/d4c24_methyl_10000b.bedGraph",
  d4c34= "./source/d4c34_methyl_10000b.bedGraph",
  NG = "./source/NG_methyl_10000b.bedGraph",
  FGO = "./source/FGO_methyl_10000b.bedGraph",
  P12 = "./source/P12_methyl_10000b.bedGraph",
  P15 ="./source/P15_methyl_10000b.bedGraph"
)

samps <- c( "d4c7", "d4c13", "d4c24", "d4c34", "NG", "P12", "P15", "FGO")

dss %>% lapply(., function(x){
  fread(x, col.names = c('chr', 'start', 'end', 'mC')) 
}) %>% 
  Reduce(function(...) left_join(..., by=c(chr= 'chr', start = 'start', end = 'end'), all.x=TRUE), .) %>%
  `colnames<-`(c('chr', 'start', 'end', names(dss))) %>%
  mutate(start = start -1) -> mc.dat

list.files('./source', pattern = 'K36me3_pool_mm10.10kb_fpm.bedGraph') -> K36me3.path

K36me3.path %>% setNames(., sub('_K36me3_pool_mm10.10kb_fpm.bedGraph', '', .)) %>%
  lapply(function(x){
  fread(sprintf('./source/%s', x)) %>% .$V4 %>% {log2(. + .1)}
}) %>% data.frame() %>%
  bind_cols(fread("./source/c13_H18_K36me3_pool_mm10.10kb_fpm.bedGraph", select = c(1:3), col.names = c('chr', 'start', 'end')), .)-> K36me3.dat 

# include K9 data
fread('./source/d4c7PGCLC_genome_10kb_fpm_matrix.bedGraph') -> CsA.dat


# most of the regions enriched K9me3 are highly repetative
left_join(K36me3.dat, mc.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  left_join(., CsA.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  na.omit() %>%
  filter(Input != 0) %>%
  mutate(K9me3 = {log2(K9me3/Input + .1)}) %>%
  mutate(K9me2 = {log2(K9me2/Input + .1)}) %>%
  ggplot(aes(x = c13_H18, y = d4c13, col = K9me3)) +
  geom_point() +
  geom_density2d() + 
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "RdYlBu")), 
    limits = c(-5, 5), # actual limits
    oob = scales::squish, # values beyond limit/out-of-bound are "squished" to the limit
    breaks = c(-5, -2.5, 0, 2.5, 5)
  ) +
  labs(x = 'K36me3 at c13', y = 'DNAme at c13')+
  theme_bw() -> p
p
ggsave(p, filename = "sf5_h1.pdf", width = 4, height = 4)

# NG ver.
left_join(K36me3.dat, mc.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  left_join(., CsA.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  na.omit() %>%
  filter(Input != 0) %>%
  mutate(K9me3 = {log2(K9me3/Input + .1)}) %>%
  mutate(K9me2 = {log2(K9me2/Input + .1)}) %>%
  ggplot(aes(x = GC_E18, y = NG, col = K9me3)) +
  geom_point() +
  geom_density2d() + 
  scale_color_gradientn(
    colors = rev(brewer.pal(11, "RdYlBu")), 
    limits = c(-5, 5), # actual limits
    oob = scales::squish, # values beyond limit/out-of-bound are "squished" to the limit
    breaks = c(-5, -2.5, 0, 2.5, 5)
  ) +
  labs(x = 'K36me3 at E18', y = 'DNAme at NG')+
  theme_bw() -> p
p
ggsave(p, filename = "sf5_h2.pdf", width = 4, height = 4)


# remove K9me3 enriched and K4me3 enriched regions
left_join(K36me3.dat, mc.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  left_join(., CsA.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  na.omit() %>%
  filter(Input != 0) %>%
  mutate(K9me3 = {log2(K9me3/Input + .1)}) %>%
  mutate(K9me2 = {log2(K9me2/Input + .1)}) %>%
  mutate(K4me3 = {log2(K4me3/Input + .1)}) %>%
  filter(K9me3 < -0.5 ) %>%
  filter(K4me3 < 0 ) -> dat

dat %>% {cor(.$d4c13, .$c13_H18)}

lm(d4c13 ~ c13_H18, data = dat) -> res
summary(res)

dat %>% 
  ggplot(aes(x = c13_H18, y = d4c13))+
  geom_hex(binwidth = c(0.1, 0.004)) + 
  colorscale +
  coord_cartesian(ylim = c(0, 0.3), xlim = c(-2, 5)) +
  geom_abline(slope = 0.0261828 , intercept = 0.0722378, col="red", alpha =.8) +
  theme_bw() +
  labs(x = 'K36me3 at c13', y = 'DNAme at c13') -> p

p
ggsave(p, filename = "sf5_h3.pdf", width = 4, height = 4)

left_join(K36me3.dat, mc.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  left_join(., CsA.dat, by=c(chr= 'chr', start = 'start', end = 'end')) %>%
  na.omit() %>%
  filter(Input != 0) %>%
  mutate(K9me3 = {log2(K9me3/Input + .1)}) %>%
  mutate(K9me2 = {log2(K9me2/Input + .1)}) %>%
  mutate(K4me3 = {log2(K4me3/Input + .1)}) %>%
  filter(K9me3 < -0.5 ) %>%
  filter(K4me3 < 0 ) -> dat

dat %>% {cor(.$GC_E18, .$NG)}


lm(NG ~ GC_E18, data = dat) -> res
summary(res)

dat %>%
  ggplot(aes(x = GC_E18, y = NG))+
  geom_hex(binwidth = c(0.1, 0.004)) + 
  colorscale +
  coord_cartesian(ylim = c(0, 0.3), xlim = c(-2, 5)) +
  geom_abline(slope = 1.186e-03 , intercept = 8.792e-03, col="red", alpha =.8) +
  theme_bw() +
  labs(x = 'K36me3 at NG', y = 'DNAme at NG') -> p

p
ggsave(p, filename = "sf5_h4.pdf", width = 4, height = 4)
