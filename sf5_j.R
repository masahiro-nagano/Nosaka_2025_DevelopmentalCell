library(tidyverse)
library(data.table)
library(rtracklayer)
library(eulerr)
library(UpSetR2)
library(ggsignif)

vivocolor=c(E9.5="#D9D9D9", E12.5="#BDBDBD" ,E14.5="#969696", E16.5="#737373", E18.5="#525252",
            Primordial="#00FFFF",Primary="#00CCD6",Secondary="#0099AD",Preantral="#006685" ,Antral="#00335C","MII" ="#000033")
vitrocolor=c(c0="#CCFFCC",c3="#CCBFD9",c5="#CC80E6",c7="#CC40F2",c9="#CC00FF",c11="#E64D80",
             c13= "#FF9900",c24="#CC4D00",c34="#990000")


# How is the DNA methylation in differential regions of oocyte specific peaks
fread('./source/gv_c34_K4me3_peaks.bed', col.names = c('chr', 'start', 'end')) %>% makeGRangesFromDataFrame() -> gv_c34
fread('./source/gv_specific_K4me3_peaks.bed', col.names = c('chr', 'start', 'end')) %>% makeGRangesFromDataFrame() -> gv_specific
fread('./source/c13_specific_K4me3_peaks.bed', col.names = c('chr', 'start', 'end')) %>% makeGRangesFromDataFrame() -> c13_specific
fread('./source/all_shared_K4me3_peaks.bed', col.names = c('chr', 'start', 'end')) %>% makeGRangesFromDataFrame() -> all_shared

# 2kb resolution mC info
dss <- c(
  d4c13 = "./source/d4c13_methyl_2000b.bedGraph",
  d4c24 = "./source/d4c24_methyl_2000b.bedGraph",
  d4c34= "./source/d4c34_methyl_2000b.bedGraph",
  NG = "./source/NG_methyl_2000b.bedGraph",
  FGO = "./source/FGO_methyl_2000b.bedGraph",
  P12 = "./source/P12_methyl_2000b.bedGraph",
  P15 ="./source/P15_methyl_2000b.bedGraph"
)

samps <- c( "d4c13", "d4c24", "d4c34", "NG", "P12", "P15", "FGO")

dss %>% lapply(., function(x){
  fread(x, col.names = c('chr', 'start', 'end', 'mC')) 
}) %>% 
  Reduce(function(...) left_join(..., by=c(chr= 'chr', start = 'start', end = 'end'), all.x=TRUE), .) %>%
  `colnames<-`(c('chr', 'start', 'end', names(dss))) %>%
  mutate(start = start -1) -> mc.dat


mc.dat %>% makeGRangesFromDataFrame() -> mc.dat.gr

# More precisely compare between 
gv %>% makeGRangesFromDataFrame() -> gv.gr
c34 %>% makeGRangesFromDataFrame() -> c34.gr

all.peaks %>%
  data.frame() %>%
  bind_cols(ol.mat) %>%
  filter(c34 == T & gv == F & e18 == F & c13 == F) %>%
  makeGRangesFromDataFrame() -> c34_only

setdiff(gv.gr, c34.gr) -> gv_c34
setdiff(c34.gr, gv.gr) -> c34_gv
intersect(gv.gr, c34.gr) -> gv_c34_intersect

mc.dat %>%
  mutate(cat = case_when(overlapsAny(mc.dat.gr, c34_only) ~ 'c34_only',
                         overlapsAny(mc.dat.gr, gv_c34) ~ 'GV_specific',
                         overlapsAny(mc.dat.gr, c34_gv) ~ 'c34_specific',
                         overlapsAny(mc.dat.gr, gv_c34_intersect) ~ 'Shared_K4me3',
                         TRUE ~ 'Other')) %>%
  filter(cat %in% c('c34_specific', 'Shared_K4me3')) %>%
  mutate(d = 100*(d4c34 - FGO)) %>%
  ggplot(aes(x = cat, y = d)) +
  geom_violin(aes(fill = cat)) +
  geom_boxplot() +
  geom_signif(comparisons =  list(c("c34_specific", "Shared_K4me3")),
              test = wilcox.test, 
              map_signif_level = T, step_increase = 0.1) +
  coord_cartesian(ylim = c(-100, 100)) +
  labs(x = '', y = 'Diff Methyl (%) (c34 - FGO)', fill = 'K4me3 peak category') +
  theme_bw() -> p

ggsave(file ='sf5_j.pdf', width = 5, height = 4)
