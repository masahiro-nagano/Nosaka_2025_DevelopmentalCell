library(data.table)
library(tidyverse)
library(rtracklayer)
library(ggsignif)
library(patchwork)

# K27me3 imprint gene analysis

# load TSS info
fread("./source/GRCm38p3_wo_pseudo_v1b_chrx_gene.bed", select = c(1:4, 5), col.names = c('chr','start','end', 'strand', 'name')) %>% 
  mutate(TSS = case_when(strand == ' +' ~ start,
                         TRUE ~ end)) %>%
  dplyr::transmute(chr, start = TSS - 1000, end = TSS + 1000, name) %>%
  arrange(chr, start, end) -> tss.info

list.files('./source/', pattern = '_pool_mm10.TSS1kb_fpm.bedGraph', full.names = T) %>%
  setNames(., sub('.*\\/(.*)_pool_mm10.TSS1kb_fpm.bedGraph', '\\1', .)) -> all
all %>%  {.[grepl(., pattern = 'K27me3')]} -> K27me3
all %>%  {.[grepl(., pattern = 'Input')]} -> Input


k27.imprint.genes <- c("Gab1", "Sfmbt2", "Bbx", "Mbnl2", "Epas1", "C430002E04Rik", "Phf17", "Etv6", "Slc38a4", "Gramd1b", "Otx2", "Tle3", "E2f3", "Smoc1", "Sox21", "Adamts2")


K27me3 %>%
  mapply(function(x, n){
    fread(x, col.names = c('chr', 'start', 'end', 'count')) %>%
      mutate(norm_cts = count/2) %>% # FPKM
      left_join(., tss.info) %>%
      filter(name %in% k27.imprint.genes) %>%
      mutate(samp = n)
  }, ., names(.), SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  remove_rownames() %>% separate(samp, into = c('cell', 'line', 'mark'), remove = F) -> imprint.k27


Input %>%
  mapply(function(x, n){
    fread(x, col.names = c('chr', 'start', 'end', 'count')) %>%
      mutate(input = count/2) %>% # FPKM
      left_join(., tss.info) %>%
      filter(name %in% k27.imprint.genes) %>%
      mutate(samp = n)
  }, ., names(.), SIMPLIFY = F) %>%
  do.call(rbind, .) %>%
  remove_rownames() %>% separate(samp, into = c('cell', 'line', 'mark')) -> imprint.input

imprint.k27 %>%
  left_join(., imprint.input, by = c(chr = 'chr', start = 'start', end = 'end', name = 'name', cell = 'cell', line = 'line')) %>%
  mutate(val = case_when(is.na(input) ~ log2(norm_cts + 0.1),
                         TRUE ~ log2((norm_cts +0.1)/(input +0.1)))) -> dat

vitro.samps <- c("c7_H18_K27me3", "c13_H18_K27me3", "c34_H18_K27me3")

dat %>%
  filter(samp %in% vitro.samps) %>%
  dplyr::select(samp, val, name) %>%
  spread(key = name, value =val) %>%
  column_to_rownames('samp') %>%
  mutate_all(function(x){
    mean(x)-> m
    sd(x) -> s
    ((x-m)/s)
  })  %>%
  rownames_to_column('samp')%>%
  gather(-samp, key = name, value = val) %>%
  mutate(samp = factor(samp, levels = vitro.samps)) -> vitro.dat

vivo.samps <- c("PGC_E13_K27me3", "E18_house_K27me3", "GV_house_K27me3")


dat %>%
  filter(samp %in% vivo.samps) %>%
  dplyr::select(samp, val, name) %>%
  spread(key = name, value =val) %>%
  column_to_rownames('samp') %>%
  mutate_all(function(x){
    mean(x)-> m
    sd(x) -> s
    ((x-m)/s)
  })  %>%
  rownames_to_column('samp')%>%
  gather(-samp, key = name, value = val) %>%
  mutate(samp = factor(samp, levels = vivo.samps)) -> vivo.dat

clrs <- c("#CC40F2", "#FF9900", "#990000", "#BDBDBD", "#525252","#00335C")

bind_rows(vitro.dat) %>%
  mutate(col = case_when(samp %in% vitro.samps ~ 'vitro',
                         samp %in% vivo.samps ~ 'vivo')) %>%
  mutate(samp = factor(samp, levels = c(vitro.samps, vivo.samps))) %>%
  ggplot(aes(x= samp, y = val)) +
  geom_violin(aes(fill = samp)) +
  geom_boxplot(width = .2) +
  geom_signif(comparisons =  list(c("c7_H18_K27me3", "c13_H18_K27me3"), c("c13_H18_K27me3", "c34_H18_K27me3"), c("c7_H18_K27me3", "c34_H18_K27me3")),
              test = wilcox.test, 
              # test.args=list(alternative = "two.sided", paired=T),
              map_signif_level = TRUE, step_increase = 0.1) +
  theme_bw() +
  facet_wrap(.~col, scales = 'free') +
  labs(x='', y = 'Z-score') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = clrs) -> p1

p1

bind_rows(vivo.dat) %>%
  mutate(col = case_when(samp %in% vitro.samps ~ 'vitro',
                         samp %in% vivo.samps ~ 'vivo')) %>%
  mutate(samp = factor(samp, levels = c(vitro.samps, vivo.samps))) %>%
  ggplot(aes(x= samp, y = val)) +
  geom_violin(aes(fill = samp)) +
  geom_boxplot(width = .2) +
  geom_signif(comparisons =  list(c("PGC_E13_K27me3", "E18_house_K27me3"), c("E18_house_K27me3", "GV_house_K27me3"), c("PGC_E13_K27me3", "GV_house_K27me3")),
              test = wilcox.test, 
              # test.args=list(alternative = "two.sided", paired=T),
              map_signif_level = TRUE, step_increase = 0.1) +
  theme_bw() +
  facet_wrap(.~col, scales = 'free') +
  labs(x='', y = 'Z-score') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90))+
  scale_fill_manual(values = clrs) -> p2
p2

p2 | p1 -> p

p
ggsave('f7_i.pdf', p, height = 6, width = 4, device = cairo_pdf, bg = 'transparent')
