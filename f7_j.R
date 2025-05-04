library(data.table)
library(tidyverse)
library(rtracklayer)
library(ggsignif)
library(patchwork)


# K27me3 imprint gene analysis
k27.imprint.genes <- c("Gab1", "Sfmbt2", "Bbx", "Mbnl2", "Epas1", "C430002E04Rik", "Phf17", "Etv6", "Slc38a4", "Gramd1b", "Otx2", "Tle3", "E2f3", "Smoc1", "Sox21", "Adamts2")

fread('./source/SC3_RNA_expression_all.tsv') %>%
  dplyr::rename(gene_name = 'V1') -> expression_all

expression_all %>%
  filter(gene_name %in% k27.imprint.genes) %>%
  dplyr::select(gene_name, E12.5_vivo_RNA,  Primordial_NGO_vivo_RNA, Preantral_P14_vivo_RNA, Antral_FGO_vivo_RNA, c7_H18_RNA, c13_H18_RNA, c34_H18_RNA) -> rna.dat

clrs <- c("#CC40F2", "#FF9900", "#990000", "#BDBDBD", "#525252","#00335C")

# zscore
rna.dat %>%
  gather(-gene_name, key = samp, value = val) %>%
  spread(key = gene_name, value =val) %>%
  column_to_rownames('samp') %>%
  mutate_all(function(x){
    mean(x)-> m
    sd(x) -> s
    ((x-m)/s)
  })  %>%
  rownames_to_column('samp')%>%
  gather(-samp, key = name, value = val) %>%
  mutate(samp = factor(samp, levels = c('c7_H18_RNA', 'c13_H18_RNA', 'c34_H18_RNA', 'E12.5_vivo_RNA',  'Primordial_NGO_vivo_RNA', 'Preantral_P14_vivo_RNA', 'Antral_FGO_vivo_RNA'))) %>%
  mutate(col = case_when(samp %in% c('c7_H18_RNA', 'c13_H18_RNA', 'c34_H18_RNA') ~ 'vitro',
                         samp %in% c('E12.5_vivo_RNA',  'Primordial_NGO_vivo_RNA','Antral_FGO_vivo_RNA') ~ 'vivo')) %>%
  na.omit() %>%
  mutate(samp = factor(samp, levels = c('c7_H18_RNA', 'c13_H18_RNA', 'c34_H18_RNA', 'E12.5_vivo_RNA',  'Primordial_NGO_vivo_RNA','Antral_FGO_vivo_RNA'))) -> dat

dat %>%
  filter(col == 'vitro') %>%
  ggplot(aes(x= samp, y = val)) +
  geom_violin(aes(fill = samp)) +
  geom_boxplot(width = .3) +
  geom_signif(comparisons =  list(c("c7_H18_RNA", "c13_H18_RNA"), c("c13_H18_RNA", "c34_H18_RNA"), c("c7_H18_RNA", "c34_H18_RNA")),
                test = wilcox.test, 
                map_signif_level = TRUE, step_increase = 0.1) +
  theme_bw() +
  labs(x='', y = 'Z-score') +
  facet_wrap(.~col, scales = 'free') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = clrs)-> p1
p1

dat %>%
  filter(col == 'vivo') %>%
  ggplot(aes(x= samp, y = val)) +
  geom_violin(aes(fill = samp)) +
  geom_boxplot(width = .3) +
  geom_signif(comparisons =  list(c("E12.5_vivo_RNA", "Primordial_NGO_vivo_RNA"), c("Primordial_NGO_vivo_RNA", "Antral_FGO_vivo_RNA"), c("E12.5_vivo_RNA", "Antral_FGO_vivo_RNA")),
              test = wilcox.test, 
              map_signif_level = TRUE, step_increase = 0.1) +
  theme_bw() +
  labs(x='', y = 'Z-score') +
  facet_wrap(.~col, scales = 'free') +
  theme(legend.position = 'none',
        axis.text.x = element_text(angle = 90)) +
  scale_fill_manual(values = clrs)-> p2

p2

p2 | p1 -> p

p
ggsave('f7_j.pdf', p, height = 6, width = 4, device = cairo_pdf, bg = 'transparent')