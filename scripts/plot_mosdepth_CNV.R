library(tidyverse)
library(purrr)
library(readr)
library(fs)
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(scales)
library(viridis)

keydiversitystrains = tibble(infiles=c("ATCC_42720","CBS_6936", "AR_0398",
"UCDFST_16-4994.55","UCDFST_61-4","UCDFST_71-129",
"UCDFST_76-31","UCDFST_79-1","UCDFST_80-11",
"UCDFST_80-12","UCDFST_82-606.2","CL_1A","A_L17",
"C_L01","B_L01","FDAARGOS_655","DSG_P1","N2_070_000G1"))

Prefix = "Supercontig_1." # fixme
mosdepthdir = "coverage/mosdepth"
patternstr=".10000bp.regions.bed.gz$"
infiles <- list.files(path=mosdepthdir,pattern = patternstr, recursive = FALSE,full.names=TRUE)
mosdepthdir <- paste(mosdepthdir,"/",sep = "")
pattern <- str_replace(patternstr,"\\$","")
mdfiles <- tibble::rowid_to_column(data.frame(infiles),"source") %>% mutate_at("infiles",str_replace,mosdepthdir,"") %>% mutate_at("infiles",str_replace,pattern,"")
mdfiles <- tibble(mdfiles)  %>% type_convert(col_types="ic")
mosdepth <- infiles %>% map_dfr(read_tsv, .id="source",
                                col_names = c("Scaffold","Start","End","Coverage"), 
                                col_types=c("ciid")) %>% filter(Scaffold != "MT_CBS_6936") %>%
  mutate_at("Scaffold",str_replace,Prefix,"") %>% type_convert(col_types="iiiid")
mosdepth <- mosdepth %>% left_join(mdfiles,by="source") %>% select(infiles,Scaffold,Start,End,Coverage) 
mosdepthfilter <- mosdepth %>% semi_join(keydiversitystrains,by="infiles")
meancoverage = mosdepthfilter %>% group_by(infiles) %>% summarize(medianDepth = median(Coverage))
mosdepthfilter <- mosdepthfilter %>% left_join(meancoverage,by="infiles") %>% mutate(normcoverage = Coverage/medianDepth)
chrlist = 1:8

lenplot <- mosdepthfilter%>% filter(Scaffold != 'NA') %>% ggplot(aes(x = Start, y = normcoverage,color = factor(infiles))) +
  geom_point() + scale_color_viridis(discrete=TRUE) + theme_bw() +
  labs(x="Position",y="Normalized coverage",
       colour="Strains",
       title = "Chrom Coverage") + facet_grid(Scaffold ~ .,scales = "free",space = "free")

lenplot
ggsave("plots/CNV_depth_plot_openscales.pdf",lenplot,width=12,height=16)
lenplot <- mosdepthfilter%>% filter(Scaffold != 'NA') %>% ggplot(aes(x = Start, y = normcoverage,color = factor(infiles))) +
  geom_point() + scale_color_viridis(discrete=TRUE) + theme_bw() +
  labs(x="Position",y="Normalized coverage",
       colour="Strains",
       title = "Chrom Coverage") + facet_grid(Scaffold ~ .)
ggsave("plots/CNV_depth_plot_fixedscales.pdf",lenplot,width=12,height=16)
heatplot <- mosdepthfilter %>% filter(Scaffold != 'NA') %>%
  ggplot(aes(x = Start, fill = normcoverage,y = infiles)) + 
  geom_tile() + theme_bw() + scale_fill_viridis(option="magma") +
  labs(title = "Coverage",x="Chrom position",y="Strain",colour ="Normalized Coverage") + 
  facet_grid(Scaffold ~ .)
heatplot
ggsave("plots/CNV_heatmap_facet1.pdf",heatplot,width=12,height=16)

mosdepthfilter <- mosdepth %>% semi_join(keydiversitystrains,by="infiles")
meancoverage = mosdepthfilter %>% group_by(infiles) %>% summarize(medianDepth = median(Coverage))
mosdepthfilter <- mosdepthfilter %>% left_join(meancoverage,by="infiles") %>% 
  mutate(normcoverage = pmin(4,Coverage/medianDepth))

lenplot <- mosdepthfilter %>% filter(Scaffold != 'NA') %>% 
  ggplot(aes(x = Start, y = normcoverage,color = factor(infiles))) +
  geom_point() + scale_color_viridis(discrete=TRUE) + theme_bw() +
  labs(x="Position",y="Normalized coverage",
       colour="Strains",
       title = "Chrom Coverage") + facet_grid(Scaffold ~ .,scales = "free",space = "free")

lenplot
ggsave("plots/CNV_depth_trunc_openscales.pdf",lenplot,width=12,height=16)
lenplot <- mosdepthfilter%>% filter(Scaffold != 'NA') %>% 
  ggplot(aes(x = Start, y = normcoverage,color = factor(infiles))) +
  geom_point() + scale_color_viridis(discrete=TRUE) + theme_bw() +
  labs(x="Position",y="Normalized coverage",
       colour="Strains",
       title = "Chrom Coverage") + facet_grid(Scaffold ~ .)
ggsave("plots/CNV_depth_trunc_fixedscales.pdf",lenplot,width=12,height=16)
heatplot <- mosdepthfilter %>% filter(Scaffold != 'NA') %>%
  ggplot(aes(x = Start, fill = normcoverage,y = infiles)) + 
  geom_tile() + theme_bw() + scale_fill_viridis(option="magma") +
  labs(title = "Coverage",x="Chrom position",y="Strain",colour ="Normalized Coverage") + 
  facet_grid(Scaffold ~ .)
heatplot
ggsave("plots/CNV_heatmap_trunc_facet1.pdf",heatplot,width=12,height=16)

