library(tidyverse)
library(data.table)
library(plotly)

# K36me3
# visualize deeptools output
dat <- read.table('./source/K36me3_epic_domain.out.tab', skip = 2, stringsAsFactors = F) %>%
  reshape2::melt(id.vars = c("V1", "V2"),
                 value.name = 'y', variable.name = 'x') %>%
  mutate(x = as.numeric(sub('V', '', x)) - 2) %>%
  dplyr::rename(reg = V2, bws = V1) %>% mutate(reg=str_replace(reg, "_peaks.bed", "")) %>% 
  separate(., bws, into=c("cell", "line", "mark", "rep")) %>%
  mutate(samp = sprintf('%s_%s', cell, line)) 
# unnorm
dat %>%  
  na.omit() %>%
  mutate(class=case_when(samp%in% c('GV_house', 'E18_house') ~ "vivo",
                         cell %in% c("PGC", "GO", "GV", "MII") ~ "vivo",
                         TRUE~"vitro")) %>%
  filter(samp %in% c("GV_house", "c13_H18", "c34_H18", 'E18_house')) %>%
  mutate(samp = factor(samp, levels = c('E18_house', "c13_H18", "GV_house", "c34_H18"))) %>%
  ggplot(., aes(x, y)) +
  geom_vline(xintercept = c(20, 70), color = "grey80",
             size = 0.5) +
  geom_line(aes(color = samp)) +
  facet_wrap(. ~ class, ncol=4) +
  ylab("Enrichment") +
  scale_x_continuous(breaks = c(1, 20, 70, 90),
                     labels = c("-200kb", "start", "end", "+200kb")) +
  scale_color_manual(values = c("#525252", "#FF9900", "#00335C", "#990000")) +
  theme(legend.position = "bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.line.x = element_line(size = 0.5, color = "black"),
        axis.text = element_text(size = 12, color = "black"),
        axis.ticks.y = element_blank(),
        axis.title = element_text(size = 12, color = "black"),
        panel.grid.major = element_line(color = "grey", linetype = "dashed"),
        strip.background = element_rect(fill = "black"),
        strip.text = element_text(size = 12, color = "white"),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 12, color = 'black'),
        plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")) ->p
p
ggsave('f7_c.pdf', p, height = 3.78, width = 5.68, device = cairo_pdf, bg = 'transparent')
