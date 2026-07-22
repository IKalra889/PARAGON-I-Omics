### Fig 5 ###
### Phototrophy:Heterotrophy ratios ###
### Deseq of KOs ###
### Daily PIT, NetTraps and WaterColumn ###
### By: Isha Kalra ###
### Last Updated: 07/21/2026 ###

library(tidyverse)
library(reshape2)
library(ggpubr)

###Daily-PITs####
#read the full normalized data (counts+tax+function) -no metazoa
df <- read.csv("metaT_dpit_wc_normalized_filtered.csv") #277486 rows

#remove fungiand late water column sample
df2 <- df %>% filter(!grepl("Fungi", taxonomy)) %>%
  mutate(water_01 = rowMeans(cbind(wc1_1, wc1_2), na.rm = TRUE)) %>%
  mutate(water_05 = rowMeans(cbind(wc2_1, wc2_2), na.rm = TRUE)) %>%
  select(-starts_with("wc"))

#read custom ko list
ko <- read.csv("ko_list_full.csv")
ko$Trophic_status <- ifelse(grepl("Phagotrophy",ko$Metabolism),"Heterotrophy",NA)
ko$Trophic_status <- ifelse(grepl("Photosynthesis|Calvin cycle",ko$Metabolism),"Phototrophy",ko$Trophic_status)
ko$Trophic_status <- ifelse(is.na(ko$Trophic_status),"Other",ko$Trophic_status)

#join the normalized data with the custom KO categories
CPM_wKO <- left_join(df2, ko, by="KO")

#no of transcripts assigned trophic status
sum(!is.na(CPM_wKO$Trophic_status)) #47033
#total
length(CPM_wKO$Trophic_status) #273568
# % of transcripts assigned 
perc_trophic_dpit <- 47033*100/273568 #17.19 %

##dpit p:h ratio dataframe
dpit_ph <- CPM_wKO %>%
  select(-geneID,-taxonomy,-description,-GeneID, -Source) %>%
  melt() %>%
  group_by(variable, Trophic_status) %>%
  summarise(cpm=sum(value)) %>%
  drop_na() %>%
  pivot_wider(names_from = Trophic_status, values_from = cpm) %>%
  mutate(
    ratio = Phototrophy / Heterotrophy,
    log2_ratio = log2(ratio)
  )

#add metadata colz
colz <- colsplit(dpit_ph$variable, "_", c("collection", "day"))
dpit_ph <- cbind(dpit_ph, colz)

#save
write.csv(dpit_ph, file = "dpit_ph_ratio.csv")

###Net-traps###

#read the full normalized data (counts+tax+function) for net traps only
df <- read.csv("metaT_net_only_normalized.csv") 

#remove fungiand late water column sample
df2 <- df %>% filter(!grepl("Fungi|Metazoa", taxonomy))

#join the normalized data with the custom KO categories
CPM_wKO_net <- left_join(df2, ko, by="KO")

## PH dataframe for NetTraps
net_ph <- CPM_wKO_net %>%
  select(-geneID,-taxonomy,-description,-GeneID, -Source) %>%
  melt() %>%
  group_by(variable, Trophic_status) %>%
  summarise(cpm=sum(value)) %>%
  drop_na() %>%
  pivot_wider(names_from = Trophic_status, values_from = cpm) %>%
  mutate(
    ratio = Phototrophy / Heterotrophy,
    log2_ratio = log2(ratio)
  )

net_ph$sample <- factor(net_ph$variable, levels=c("NetTrap_150_1", "NetTrap_150_2", "NetTrap_150_3", 
                                                  "NetTrap_175_1", "NetTrap_175_2", 
                                                  "NetTrap_200_1", "NetTrap_200_2", "NetTrap_200_3", 
                                                  "NetTrap_300_1", "NetTrap_300_2", "NetTrap_300_3"),
                        labels=c("net150_1","net150_2","net150_3",
                                 "net175_1","net175_2",
                                 "net200_1", "net200_2", "net200_3",
                                 "net300_1", "net300_2", "net300_3"))

#add collection metadata
colz <- colsplit(net_ph$sample, "_", c("collection", "rep"))
net_ph <- cbind(net_ph, colz)

##add both dpit and net trap
all_ph <- rbind(dpit_ph, net_ph)
all_ph <- all_ph %>% group_by(collection) %>% 
  mutate(mean = mean(log2_ratio))

##plot

#at 150 m particle vs water
fig_5a <- all_ph %>%
  filter(collection %in% c("dpit", "water", "net150")) %>%
  ggplot(aes(x=collection, y = log2_ratio, fill = collection)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8)+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))+
  labs(
    title = "",
    y = expression(log[2](Phototrophy:Heterotrophy)),
    x = ""
  )+
  scale_fill_manual(values = c("salmon", "cyan", "#0099FF"))+
  scale_x_discrete(labels = c("dpit" = "Daily PIT 150m", "water" = "WaterColumn 150m", "net150" = "NetTrap 150m"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust =1))+
  scale_y_continuous(limits = c(-5, 0))+
  theme(legend.position = "none")

# net traps at different depths
fig_5c <- all_ph %>%
  filter(grepl("net", collection)) %>%
  ggplot(aes(x=collection, y = log2_ratio, fill=collection)) +
  geom_boxplot()+
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))+
  labs(
    title = "",
    y = expression(log[2](Phototrophy:Heterotrophy)),
    x = ""
  )+
  scale_x_discrete(labels = c("net150" = "NetTrap 150m", "net175" = "NetTrap 175m", 
                              "net200" = "NetTrap 200m", "net300" = "NetTrap 300m"))+
  scale_fill_manual(values = c("cyan", "pink", "lightgreen", "purple"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust =1))+
  scale_y_continuous(limits = c(-5, 0))+
  theme(legend.position = "none")

##statistics

#anova test on all samples
anova_all_ph <- aov(log2_ratio ~ collection, data = all_ph)
summary(anova_all_ph)

#Tukey's post hoc test
TukeyHSD(anova_all_ph)

#### Taxa specific P:H ratios ####

## DPIT + WC
taxa_list <- c("Chlorophyta","Chlorarachniophyta", "Dinophyceae","Polycystinea" ,"Phaeodaria", "Foraminifera")

#function to generate P:H ratio
tax_ph <- function(df, taxon) {
  df %>%
    drop_na() %>%
    filter(grepl(taxon, taxonomy)) %>%
    select(-geneID, -taxonomy, -description, -GeneID, -Source) %>%
    melt() %>%
    group_by(variable, Trophic_status) %>%
    summarise(cpm = sum(value)) %>%
    pivot_wider(names_from = Trophic_status, values_from = cpm) %>%
    mutate(
      ratio = Phototrophy / Heterotrophy,
      log2_ratio = log2(ratio),
      tax = taxon
    )
}

## Taxa specific data frame with P:H ratio for DPIT and WC
Chloro_dpit_ph <- tax_ph(CPM_wKO, "Chlorophyta")
Chlorarch_dpit_ph <- tax_ph(CPM_wKO, "Chlorarachniophyta")
Dino_dpit_ph <- tax_ph(CPM_wKO, "Dinophyceae")
Polycys_dpit_ph <- tax_ph(CPM_wKO, "Polycystinea")
Phaeo_dpit_ph <- tax_ph(CPM_wKO, "Phaeodaria")
Foram_dpit_ph <- tax_ph(CPM_wKO, "Foraminifera")
dpit_tax_ph <- rbind(Chloro_dpit_ph, Chlorarch_dpit_ph, Dino_dpit_ph, Polycys_dpit_ph, Phaeo_dpit_ph, Foram_dpit_ph)
#add metadata colz
colz <- colsplit(dpit_tax_ph$variable, "_", c("collection", "day"))

## dataframe with all taxa for daily PITs and water column
dpit_tax_ph <- cbind(dpit_tax_ph, colz)

## Taxa specific data frame with P:H ratio for NetTraps

Chloro_net_ph <- tax_ph(CPM_wKO_net, "Chlorophyta")
Chlorarch_net_ph <- tax_ph(CPM_wKO_net, "Chlorarachniophyta")
Dino_net_ph <- tax_ph(CPM_wKO_net, "Dinophyceae")
Polycys_net_ph <- tax_ph(CPM_wKO_net, "Polycystinea")
Phaeo_net_ph <- tax_ph(CPM_wKO_net, "Phaeodaria")
Foram_net_ph <- tax_ph(CPM_wKO_net, "Foraminifera")
net_tax_ph <- rbind(Chloro_net_ph, Chlorarch_net_ph, Dino_net_ph, Polycys_net_ph, Phaeo_net_ph, Foram_net_ph)

#relabel the samples
net_tax_ph$sample <- factor(net_tax_ph$variable, levels=c("NetTrap_150_1", "NetTrap_150_2", "NetTrap_150_3", 
                                                  "NetTrap_175_1", "NetTrap_175_2", 
                                                  "NetTrap_200_1", "NetTrap_200_2", "NetTrap_200_3", 
                                                  "NetTrap_300_1", "NetTrap_300_2", "NetTrap_300_3"),
                        labels=c("net150_1","net150_2","net150_3",
                                 "net175_1","net175_2",
                                 "net200_1", "net200_2", "net200_3",
                                 "net300_1", "net300_2", "net300_3"))
#add metadata colz
colz <- colsplit(net_tax_ph$sample, "_", c("collection", "rep"))

## dataframe with all taxa for NetTraps
net_tax_ph <- cbind(net_tax_ph, colz)

#combine both df
all_tax_ph <- rbind(dpit_tax_ph, net_tax_ph)
all_tax_ph$taxa <- factor(all_tax_ph$tax, 
                          levels=c("Chlorophyta","Chlorarachniophyta", "Dinophyceae",
                                   "Polycystinea" ,"Phaeodaria", "Foraminifera"),
                          labels=c("Chlorophytes","Chlorarachniophytes", "Dinoflagellates",
                                   "Polycystines" ,"Phaeodarians", "Foraminiferans"))

labels <- c(
  "Chlorophytes" = "B Chlorophytes",
  "Chlorarachniophytes" = "C Chlorarachniophytes",
  "Dinoflagellates" = "D Dinoflagellates",
  "Polycystines" = "E Polycystines",
  "Phaeodarians" = "F Phaeodarians",
  "Foraminiferans" = "G Foraminiferans"
)
## plot the taxa specific P:H ratios

fig_5b <- all_tax_ph %>%
  filter(collection %in% c("dpit", "water", "net150")) %>%
  ggplot(aes(x=collection, y = log2_ratio, fill=collection)) +
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  facet_wrap(~taxa)+
  scale_x_discrete(labels = c("dpit" = "Daily PIT 150m", "water" = "WaterColumn 150m", "net150" = "NetTrap 150m"))+
  labs(
    title = "",
    y = expression(log[2](Phototrophy:Heterotrophy)),
    x = ""
  ) +
  scale_fill_manual(values = c("salmon", "cyan", "#0099FF"))+
  scale_x_discrete(labels = c("dpit" = "Daily PIT 150m", "water" = "WaterColumn 150m", "net150" = "NetTrap 150m"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust =1))+
  scale_y_continuous(limits = c(-15, 5))+
  theme(legend.position = "none")+
  theme(
    strip.text = element_text(size = 6, face = "bold")
  )

fig_5d <- all_tax_ph %>%
  filter(grepl("net", collection)) %>%
  ggplot(aes(x=collection, y = log2_ratio, fill = collection)) +
  geom_boxplot()+
  geom_hline(yintercept = 0, linetype = "dashed", color = "black")+
  facet_wrap(~taxa)+
  scale_x_discrete(labels = c("net150" = "NetTrap 150m", "net175" = "NetTrap 175m", 
                              "net200" = "NetTrap 200m", "net300" = "NetTrap 300m"))+
  labs(
    title = "",
    y = expression(log[2](Phototrophy:Heterotrophy)),
    x = ""
  ) +
  scale_fill_manual(values = c("cyan", "pink", "lightgreen", "purple"))+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust =1))+
  scale_y_continuous(limits = c(-15, 5))+
  theme(legend.position = "none")+
  theme(
    strip.text = element_text(size = 6, face = "bold")
  )

## Statistics on taxa

#function to calculate anova for each taxa
tax_anova <- function(df, taxon){
  df %>%
    filter(taxa == taxon) %>%
    filter(is.finite(log2_ratio)) %>%
    aov(log2_ratio ~ collection, data = .)
}

#Chlorophytes
chloro_anova <- tax_anova(all_tax_ph, "Chlorophytes")
summary(chloro_anova)
chloro_tukey <- TukeyHSD(chloro_anova)

#Chlorachniophytes
chlorarach_anova <- tax_anova(all_tax_ph, "Chlorarachniophytes")
summary(chlorarach_anova)
chlorarach_tukey <- TukeyHSD(chlorarach_anova)

#Dinoflagellates
dino_anova <- tax_anova(all_tax_ph, "Dinoflagellates")
summary(dino_anova)
dino_tukey <- TukeyHSD(dino_anova)

#Polycystines
poly_anova <- tax_anova(all_tax_ph, "Polycystines")
summary(poly_anova)
poly_tukey <- TukeyHSD(poly_anova)

#Phaeodarians
phaeo_anova <- tax_anova(all_tax_ph, "Phaeodarians")
summary(phaeo_anova)
phaeo_tukey <- TukeyHSD(phaeo_anova)

#Foraminiferans
foram_anova <- tax_anova(all_tax_ph, "Foraminiferans")
summary(foram_anova)
foram_tukey <- TukeyHSD(foram_anova)

#plot figures together
# Fig 5
(fig_5a + fig_5b) /
  (fig_5c + fig_5d) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave("fig5abcd.pdf", width = 4.0, height = 8.0, units = "in")

## -------------- Volcano plots - Figure 5E-G -----------------------###
## -------- Fig 5E- DPIT vs Water Column -----------###

df_vol <- read.csv("deseq_dpit_vs_water_kegg_custom.csv")

#set factor
df_vol$metabolism <- factor(df_vol$Metabolism_2,
                                       levels=c("Photosynthesis", "Calvin cycle",
                                                 "C metabolism", "Energy Acquisition",
                                                 "Nitrogen metabolism", "P metabolism",
                                                 "Phagotrophy", "Phagotrophy-other"),
                            labels=c("Photosynthesis", "Calvin cycle",
                                     "Carbon metabolism", "Energy acquisition",
                                     "Nitrogen metabolism", "Phosphorus metabolism",
                                     "Phagotrophy", "Phagotrophy"))
#set colors
colors_metabolism <- c("#66A61E", "#E7298A", "blue", 
                       "brown", "purple","lightblue", "#C4A445", "lightgrey")

name_metabolism <- c("Photosynthesis", "Calvin cycle",
                     "Carbon metabolism", "Energy acquisition",
                     "Nitrogen metabolism", "Phosphorus metabolism",
                     "Phagotrophy", "")

names(colors_metabolism) <- name_metabolism

#plot
fig_5e <- df_vol %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = metabolism), alpha = 0.5)+
  geom_text_repel(
    aes(label = gene_label),
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.6,
    segment.curvature = -0.1,
    segment.color = 'grey30')+
  geom_text_repel(aes(label = poc_label),
    size = 3.5,
    color='red',
    max.overlaps = 40,
    segment.curvature = -0.1, # Adds a nice curve to lines
    segment.color ='red')+
  scale_color_manual(values=colors_metabolism)+
  theme_bw() +
  labs(
    x = "log2 Fold Change (Daily PITs vs WaterColumn)",
    y = "-log10 Adjusted P-value",
    color="")


## ----------- Fig 5 F and G - Net trap volcanos ----------- ##
df_vol2 <- read.csv("deseq_net200_vs150_customkegg.csv")

#set factor
df_vol2$metabolism <- factor(df_vol2$Metabolism_2,
                            levels=c("Photosynthesis", "Calvin cycle",
                                     "C metabolism", "Energy Acquisition",
                                     "Nitrogen metabolism", "P metabolism",
                                     "Phagotrophy", "Phagotrophy-other"),
                            labels=c("Photosynthesis", "Calvin cycle",
                                     "Carbon metabolism", "Energy acquisition",
                                     "Nitrogen metabolism", "Phosphorus metabolism",
                                     "Phagotrophy", "Phagotrophy"))

#plot
fig_5f <- df_vol2 %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = metabolism), alpha = 0.5)+
  geom_text_repel(aes(label = label_gene),
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.6,
    segment.curvature = -0.1,
    segment.color = 'grey30')+
  scale_color_manual(values=colors_metabolism)+
  theme_bw() +
  labs(
    x = "log2 Fold Change (NetTrap 200m vs NetTrap 150m)",
    y = "-log10 Adjusted P-value",
    color="") + theme (legend.position = "none")

# fig 5g
df_vol3 <- read.csv("deseq_net300_vs150_customkegg.csv")

#set factor
df_vol3$metabolism <- factor(df_vol3$Metabolism_2,
                             levels=c("Photosynthesis", "Calvin cycle",
                                      "C metabolism", "Energy Acquisition",
                                      "Nitrogen metabolism", "P metabolism",
                                      "Phagotrophy", "Phagotrophy-other"),
                             labels=c("Photosynthesis", "Calvin cycle",
                                      "Carbon metabolism", "Energy acquisition",
                                      "Nitrogen metabolism", "Phosphorus metabolism",
                                      "Phagotrophy", "Phagotrophy"))

#plot
fig_5g <- df_vol3 %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = metabolism), alpha = 0.5)+
  geom_text_repel(aes(label = label_gene),
    size = 2.5,
    max.overlaps = 20,
    box.padding = 0.6,
    segment.curvature = -0.1,
    segment.color = 'grey30')+
  scale_color_manual(values=colors_metabolism)+
  theme_bw() +
  labs(x = "log2 Fold Change (NetTrap 300m vs NetTrap 150m)",
    y = "-log10 Adjusted P-value",
    color="") +
  theme(legend.position = "none")

# plot together
(fig_5f / fig_5g) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("f","g"))) & 
  theme(plot.tag = element_text(face = "bold", size = 14))

