### PARAGON-I manuscript Figure 2 ###
### Taxonomy plots for mRNA data ###
### Daily-PITS, Net-trap and Water-column ###
### By : Isha Kalra ###
### Last updated: 02/18/2026 ###

library(tidyr)
library(dplyr)
library(reshape2)
library(ggplot2)
library(patchwork)
library(randomcoloR)

setwd("~/Desktop/Caron_lab_research/SCOPE/PARAGON-I/manuscript/analyses/compiled_data/")

## ------------------------------------------------------------------------- ##
## ---------------------------- Fig 2A ------------------------------------- ##
## Daily-PITs & WaterColumn taxa plots

df <- read.csv("metaT_dpit_wc_normalized_filtered.csv") #read protist+fungi normalized metaT data

#remove columns not needed
df$geneID <- NULL
df$KO <- NULL
df$description <- NULL

#remove fungi, average water-column data for each day, drop na values
df <- df %>%
  filter(!grepl("Fungi", taxonomy)) %>%
  mutate(
    WaterColumn_day1 = rowMeans(cbind(wc1_1, wc1_2), na.rm = TRUE),
    WaterColumn_day5 = rowMeans(cbind(wc2_1, wc2_2), na.rm = TRUE)
  ) %>%
  select(-starts_with("wc")) %>%
  drop_na() #180213 entries

#select taxa to plot and rename
df$tax <- ifelse(grepl("Chlorophyta",df$taxonomy),"Archaeplastid-Chlorophytes",NA)
df$tax <- ifelse(grepl("Apicomplexa",df$taxonomy),"Alveolates-Apicomplex",df$tax)
df$tax <- ifelse(grepl("Dinophyceae",df$taxonomy),"Alveolates-Dinoflagellates",df$tax)
df$tax <- ifelse(grepl("Syndinians",df$taxonomy),"Alveolates-Syndiniales",df$tax)
df$tax <- ifelse(grepl("Ciliate",df$taxonomy),"Alveolates-Ciliates",df$tax) 
df$tax <- ifelse(grepl("Amoebozoa",df$taxonomy),"Amoebozoa",df$tax) 
df$tax <- ifelse(grepl("Cryptophyta",df$taxonomy),"Cryptophytes",df$tax)
df$tax <- ifelse(grepl("Discoba",df$taxonomy),"Excavates-Discobids",df$tax)
df$tax <- ifelse(grepl("Stramenopile",df$taxonomy),"Stramenopiles-Other",df$tax)
df$tax <- ifelse(grepl("Pelagophyceae",df$taxonomy),"Stramenopiles-Pelagophytes",df$tax)
df$tax <- ifelse(grepl("Bacillariophyceae",df$taxonomy),"Stramenopiles-Diatoms",df$tax)
df$tax <- ifelse(grepl("MAST",df$taxonomy),"Stramenopiles-MAST",df$tax)
df$tax <- ifelse(grepl("Haptophyta",df$taxonomy),"Haptophytes",df$tax)
df$tax <- ifelse(grepl("Rhizaria",df$taxonomy),"Rhizaria-Other",df$tax)
#df$tax <- ifelse(grepl("Radiolaria",df$taxonomy),"Rhizaria-Radiolaria",df$tax)
df$tax <- ifelse(grepl("Polycystinea",df$taxonomy),"Rhizaria-Polycystines",df$tax)
df$tax <- ifelse(grepl("Acantharia",df$taxonomy),"Rhizaria-Acantharians",df$tax)
df$tax <- ifelse(grepl("Foraminifera",df$taxonomy),"Rhizaria-Foraminiferans",df$tax)
df$tax <- ifelse(grepl("Cercozoa",df$taxonomy),"Rhizaria-Cercozoans",df$tax)
#df$tax <- ifelse(grepl("Labyrinthulomycota",df$taxonomy),"Labyrinthulomycota",df$tax)
df$tax <- ifelse(grepl("Fungi",df$taxonomy),"Fungi",df$tax)
df$tax <- ifelse(grepl("Choano", df$taxonomy), "Choanoflagellates", df$tax)
#df$tax <- ifelse(grepl("Alveolate",df$taxonomy) & is.na(df$tax),"Alveolates-Other",df$tax)
df$tax <- ifelse(grepl("Archaeplastida",df$taxonomy) & is.na(df$tax),"Archaeplastid-Other",df$tax)
df$tax <- ifelse(df$taxonomy=="Eukaryota","Unknown Eukaryote",df$tax)
df$tax <- ifelse(is.na(df$tax),"Other Eukaryotes",df$tax)

#remove unidentified opisthokont and taxonomy column
df2_dpit <- df %>% 
  filter(taxonomy != "Opisthokont;") %>%
  select(-taxonomy) #176972 entries

#melt the dataframe, calculate relative abundance and total relabund for selected taxa
dfMelt_dpit <- df2_dpit %>%
  melt() %>%
  group_by(variable) %>%
  mutate(rel = 100* value/sum(value)) %>%
  ungroup() %>%
  group_by(tax, variable) %>%
  summarise(total_rel = sum(rel)) #220 observations

#taxa color
tax_color <- c('firebrick4','indianred1','tomato3','red3',
               'forestgreen','yellowgreen','purple',
               'darkblue','blue','magenta','lightblue','pink','grey',
               'gold2','moccasin','gold1','yellow3','yellow',
               'tan1','#DDAB4B','tan3','tan4')

#list of taxa
taxa_list <- c("Alveolates-Ciliates","Alveolates-Dinoflagellates","Alveolates-Syndiniales","Alveolates-Apicomplex",
               "Archaeplastid-Chlorophytes","Archaeplastid-Other","Choanoflagellates",
               "Cryptophytes","Excavates-Discobids","Fungi","Haptophytes","Labyrinthulomycota","Other Eukaryotes",
               "Rhizaria-Acantharians","Rhizaria-Cercozoans","Rhizaria-Foraminiferans","Rhizaria-Other",
               "Rhizaria-Polycystines","Stramenopiles-Diatoms","Stramenopiles-MAST",
               "Stramenopiles-Other","Stramenopiles-Pelagophytes")

#assign colors to each taxa
names(tax_color) <- taxa_list

#set variable names
dfMelt_dpit$variable <- factor(dfMelt_dpit$variable, 
                               levels=c("dpit_01","dpit_02",
                                        "dpit_03","dpit_05","dpit_06","dpit_07",
                                        "dpit_08","dpit_09","dpit_10", "WaterColumn_day1","WaterColumn_day5"),
                               labels=c("Daily PIT day1","Daily PIT day2",
                                        "Daily PIT day3", "Daily PIT day5","Daily PIT day6","Daily PIT day7",
                                        "Daily PIT day8","Daily PIT day9","Daily PIT day10", "WaterColumn day1","WaterColumn day5"))

## taxa bar plot
dpit_tax <- ggplot(dfMelt_dpit, aes(x=variable, y=total_rel, fill=tax))+
  geom_bar(colour = "black", stat="identity", position="fill")+
  labs(title="", x="",y="% relative transcript abundance")+
  theme(legend.title=element_blank(),legend.position="top",legend.text.align=0, 
        axis.text = element_text(color="black"),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line())+
  scale_y_continuous(position = "left")+
  scale_fill_manual(values=tax_color)+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(axis.text.x = element_text(angle=45,size=12, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title.x = element_text(size=14))+
  scale_y_continuous(labels = scales::percent)+
  labs(title="") + 
  geom_vline(xintercept = 9.5, linetype = "dashed", color = "black")

## ------------------------------------------------------------------------- ##
## ---------------------------- Fig 2B ------------------------------------- ##

### NetTrap & WaterColumn taxa plot

df <- read.csv("metaT_net_wc_normalized.csv") #584137 obs

#remove columns not needed
df$geneID <- NULL
df$KO <- NULL
df$description <- NULL

#remove fungi, average water-column and net trap data for each day, drop na values
df <- df %>% 
  filter(!grepl("Fungi", taxonomy)) %>%
  mutate(
    WaterColumn_150m = rowMeans(cbind(WaterColumn_1_1, WaterColumn_1_2, WaterColumn_2_1,
                                      WaterColumn_2_2), na.rm = TRUE),
    NetTrap_150m = rowMeans(cbind(NetTrap_150_1, NetTrap_150_2, NetTrap_150_3), na.rm = TRUE),
    NetTrap_175m = rowMeans(cbind(NetTrap_175_1, NetTrap_175_2), na.rm = TRUE),
    NetTrap_200m = rowMeans(cbind(NetTrap_200_1, NetTrap_200_2, NetTrap_200_3), na.rm = TRUE),
    NetTrap_300m = rowMeans(cbind(NetTrap_300_1, NetTrap_300_2, NetTrap_300_3), na.rm = TRUE)
  ) %>%
  select(-ends_with(c("_1","_2","_3"))) %>%
  drop_na()

#select taxa to plot and rename
#select taxa to plot and rename
df$tax <- ifelse(grepl("Chlorophyta",df$taxonomy),"Archaeplastid-Chlorophytes",NA)
df$tax <- ifelse(grepl("Apicomplexa",df$taxonomy),"Alveolates-Apicomplex",df$tax)
df$tax <- ifelse(grepl("Dinophyceae",df$taxonomy),"Alveolates-Dinoflagellates",df$tax)
df$tax <- ifelse(grepl("Syndinians",df$taxonomy),"Alveolates-Syndiniales",df$tax)
df$tax <- ifelse(grepl("Ciliate",df$taxonomy),"Alveolates-Ciliates",df$tax) 
df$tax <- ifelse(grepl("Amoebozoa",df$taxonomy),"Amoebozoa",df$tax) 
df$tax <- ifelse(grepl("Cryptophyta",df$taxonomy),"Cryptophytes",df$tax)
df$tax <- ifelse(grepl("Discoba",df$taxonomy),"Excavates-Discobids",df$tax)
df$tax <- ifelse(grepl("Stramenopile",df$taxonomy),"Stramenopiles-Other",df$tax)
df$tax <- ifelse(grepl("Pelagophyceae",df$taxonomy),"Stramenopiles-Pelagophytes",df$tax)
df$tax <- ifelse(grepl("Bacillariophyceae",df$taxonomy),"Stramenopiles-Diatoms",df$tax)
df$tax <- ifelse(grepl("MAST",df$taxonomy),"Stramenopiles-MAST",df$tax)
df$tax <- ifelse(grepl("Haptophyta",df$taxonomy),"Haptophytes",df$tax)
df$tax <- ifelse(grepl("Rhizaria",df$taxonomy),"Rhizaria-Other",df$tax)
#df$tax <- ifelse(grepl("Radiolaria",df$taxonomy),"Rhizaria-Radiolaria",df$tax)
df$tax <- ifelse(grepl("Polycystinea",df$taxonomy),"Rhizaria-Polycystines",df$tax)
df$tax <- ifelse(grepl("Acantharia",df$taxonomy),"Rhizaria-Acantharians",df$tax)
df$tax <- ifelse(grepl("Foraminifera",df$taxonomy),"Rhizaria-Foraminiferans",df$tax)
df$tax <- ifelse(grepl("Cercozoa",df$taxonomy),"Rhizaria-Cercozoans",df$tax)
#df$tax <- ifelse(grepl("Labyrinthulomycota",df$taxonomy),"Labyrinthulomycota",df$tax)
df$tax <- ifelse(grepl("Fungi",df$taxonomy),"Fungi",df$tax)
df$tax <- ifelse(grepl("Choano", df$taxonomy), "Choanoflagellates", df$tax)
#df$tax <- ifelse(grepl("Alveolate",df$taxonomy) & is.na(df$tax),"Alveolates-Other",df$tax)
df$tax <- ifelse(grepl("Archaeplastida",df$taxonomy) & is.na(df$tax),"Archaeplastid-Other",df$tax)
df$tax <- ifelse(df$taxonomy=="Eukaryota","Unknown Eukaryote",df$tax)
df$tax <- ifelse(is.na(df$tax),"Other Eukaryotes",df$tax)

#remove unidentified opisthokont and taxonomy column
df2_net <- df %>% 
  filter(taxonomy != "Opisthokont;") %>%
  select(-taxonomy) #419649 entries

#melt the dataframe, calculate relative abundance and total relabund for selected taxa
dfMelt_net <- df2_net %>%
  melt() %>%
  group_by(variable) %>%
  mutate(rel = 100* value/sum(value)) %>%
  ungroup() %>%
  group_by(tax, variable) %>%
  summarise(total_rel = sum(rel)) #100 observations

#reverse the levels for the sample names
dfMelt_net$variable <- factor(dfMelt_net$variable, 
                              levels=c("NetTrap_150m","NetTrap_175m", "NetTrap_200m", "NetTrap_300m", "WaterColumn_150m"),
                              labels=c("NetTrap 150m", "NetTrap 175m", "NetTrap 200m", "NetTrap 300m", "WaterColumn 150m"))

# plot
net_tax <- ggplot(dfMelt_net, aes(x=variable, y=total_rel, fill=tax))+
  geom_bar(colour = "black", stat="identity", position="fill")+
  labs(title="", x="",y="% relative transcript abundance")+
  theme(legend.title=element_blank(),legend.position="top",legend.text.align=0, 
        axis.text = element_text(color="black"),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),panel.background = element_blank(),
        panel.border = element_blank(), axis.line = element_line())+
  scale_y_continuous(position = "left")+
  scale_fill_manual(values=tax_color)+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(axis.text.x = element_text(angle=45,size=12, hjust = 1), 
        axis.text.y = element_text(size=12), axis.title.x = element_text(size=14))+
  scale_y_continuous(labels = scales::percent)+
  labs(title="") + 
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")

## -------------------------------------------------------------------------------------------------- ##

##plot both taxa plots
(dpit_tax+net_tax)+
  plot_layout(guides = "collect")+
  plot_annotation(tag_levels = "A") & 
  theme(plot.tag = element_text(face = "bold", size = 14))

#save figure
ggsave("Fig_2.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)

