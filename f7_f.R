library(tidyverse)
library(data.table)
library(rtracklayer)
library(UpSetR2)

vivocolor=c(E9.5="#D9D9D9", E12.5="#BDBDBD" ,E14.5="#969696", E16.5="#737373", E18.5="#525252",
            Primordial="#00FFFF",Primary="#00CCD6",Secondary="#0099AD",Preantral="#006685" ,Antral="#00335C","MII" ="#000033")
vitrocolor=c(c0="#CCFFCC",c3="#CCBFD9",c5="#CC80E6",c7="#CC40F2",c9="#CC00FF",c11="#E64D80",
             c13= "#FF9900",c24="#CC4D00",c34="#990000")

fread('./source/c13_H18_K4me3_peaks_1000.bed', col.names = c('chr', 'start', 'end')) -> c13

fread('./source/c34_H18_K4me3_peaks_1000.bed', col.names = c('chr', 'start', 'end')) -> c34

fread('./source/E18_house_K4me3_peaks_1000.bed', col.names = c('chr', 'start', 'end')) -> e18

fread('./source/GV_house_K4me3_peaks_1000.bed', col.names = c('chr', 'start', 'end')) -> gv

bind_rows(c13, c34, e18, gv) %>%
  makeGRangesFromDataFrame() %>%
  reduce() -> all.peaks

list(e18 = e18, c13 = c13,  gv = gv, c34 = c34) %>%
  lapply(., function(x){
    x %>% makeGRangesFromDataFrame() %>%
    {overlapsAny(all.peaks, .)}
  }) %>% data.frame() -> ol.mat

# upset
ol.mat %>%
  mutate_all(function(x){x=case_when(x== TRUE ~ 1,
                                     x == FALSE ~ 0)}) -> ol.mat2

pdf('f7_f.pdf',width=4.73,height=5.5)
print(upset(ol.mat2,
            order.by = "freq",
            show.numbers = F,
            keep.order = T,
            sets = colnames(ol.mat2),
            sets.bar.color = c("#525252", "#FF9900", "#00335C", "#990000"),
            set_size.angles = 60,
            ncut = 5))
dev.off()