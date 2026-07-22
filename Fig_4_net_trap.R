### PARAGON-I NetTrap plots - Figure 4 ###
### Ordination, Heat map and bubble plot###
### By: Isha Kalra ###
### Last updated: 07/21/2026 ###

library(tidyverse)
library(reshape2)
library(vegan)
library(pheatmap)
library(edgeR)

#### ---------------------------------------- Fig. 4A - PCA plot ---------------------------------------- ####
## read the full raw count data (counts+tax+function)
df <- read.csv("metaT_net_wc_wide.csv")

## remove fungi , metazoa, transcripts with no KO id and daily PIT day 4 sample
tmp <- df %>%
  filter(!grepl("Fungi|Metazoa", taxonomy)) %>%
  filter(KO != "")

colnames(tmp) <- c("geneID", "taxonomy", "description", "KO", "NetTrap_150_1","NetTrap_150_2", "NetTrap_150_3",
                   "NetTrap_175_1", "NetTrap_175_2", "NetTrap_200_1", "NetTrap_200_2", "NetTrap_200_3",
                   "NetTrap_300_1", "NetTrap_300_2", "NetTrap_300_3", "WaterColumn_1_1", "WaterColumn_1_2",
                   "WaterColumn_2_1", "WaterColumn_2_2")
# remove columns we don't need
tmp$geneID <- NULL
tmp$taxonomy <- NULL
tmp$description <- NULL

# group by kegg annotation
tmpSum <- tmp %>% 
  group_by(KO) %>% 
  summarise_all(sum) %>% 
  as.data.frame()

# average water column data
tmpSum <- tmpSum %>% 
  mutate(
  WaterColumn_150_day1 = rowMeans(cbind(WaterColumn_1_1, WaterColumn_1_2), na.rm = TRUE),
  WaterColumn_150_day5 = rowMeans(cbind(WaterColumn_2_1, WaterColumn_2_2), na.rm = TRUE)) %>%
  select(-WaterColumn_1_1, -WaterColumn_1_2, -WaterColumn_2_1, -WaterColumn_2_2)

#make rownames as KO
rownames(tmpSum) <- tmpSum$KO
tmpSum$KO <- NULL

# remove low counts , transform to make rows samples, and columns as KO
tmpSum <- as.data.frame(t(tmpSum))
tmpSum <- tmpSum[, colSums(tmpSum) > 0]

library(compositions)
# calculate proportions for clr transcformations (a.k.a relative abundance)
prop <- sweep(tmpSum, 1, rowSums(tmpSum), "/")
pseudocount <- min(prop[prop > 0]) / 2

# clr matrix
clr_mat <- clr(prop + pseudocount)

# pca
pca <- prcomp(clr_mat)

# look at the pca
summary(pca)

#sample wide results
pca$x

# extract the results for the plots
pca.plot <- as.data.frame(pca$x)
colz <- colsplit(rownames(pca.plot), "_", c("collection", "depth", "day/rep"))
# factor data
colz$collection <- factor(colz$collection, levels=c("NetTrap", "WaterColumn"),
                          labels=c("NetTrap", "WaterColumn"))
colz$depth <- factor(colz$depth, levels=c("150","175","200","300"),
                     labels=c("150m","175m","200m","300m"))
# add to pca df
pca.plot$collection <- colz$collection
pca.plot$depth <- colz$depth

fig_4a <- ggplot(pca.plot,aes(PC1,PC2,shape = collection, fill=depth))+
  geom_point(size = 6, aes(color=depth))+
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  labs(title = "", x="PC1: 38.91% variance", y="PC2: 11.77% variance",
       fill=NULL, shape="Sample Type", color="Depth")+
  scale_fill_manual(values = c("cyan3","pink","lightgreen","purple"))+
  scale_color_manual(values = c("cyan3","pink","lightgreen","purple"))+
  theme_bw()+
  theme(axis.text = element_blank(), legend.position = "right")+
  guides(fill = "none")

#save figure
ggsave("Fig_4a.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)

#### ---------------------------------------- Fig. 4B - Heat map ---------------------------------------- ####
#read the full normalized data (counts+tax+function) -no metazoa, only NetTrap samples
df <- read.csv("metaT_net_only_normalized.csv") #

#remove fungi
df2 <- df %>% filter(!grepl("Fungi", taxonomy))

#average data
df2 <- df2 %>%
  mutate(
    NetTrap_150m = rowMeans(cbind(NetTrap_150_1, NetTrap_150_2, NetTrap_150_3), na.rm = TRUE),
    NetTrap_175m = rowMeans(cbind(NetTrap_175_1, NetTrap_175_2), na.rm = TRUE),
    NetTrap_200m = rowMeans(cbind(NetTrap_200_1, NetTrap_200_2, NetTrap_200_3), na.rm = TRUE),
    NetTrap_300m = rowMeans(cbind(NetTrap_300_1, NetTrap_300_2, NetTrap_300_3), na.rm = TRUE)
  ) %>%
  select(-ends_with(c("_1","_2","_3"))) %>%
  drop_na()

#read custom ko list
ko <- read.csv("ko_list_full.csv")

#join the normalized data with the custom KO categories
CPM_wKO <- inner_join(df2, ko, by="KO")

#select samples and pathway
CPM_wKO <- CPM_wKO %>% select(starts_with("NetTrap"), Pathway)

# summarize by pathway
sum_cpm_wko <- CPM_wKO %>%
  group_by(Pathway) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  as.data.frame()

colnames(sum_cpm_wko) <- c("Pathway", "NetTrap 150m", "NetTrap 175m", "NetTrap 200m", "NetTrap 300m")

rm <- c("else", "por", "else", "other", "additional breakdown" )
sum_cpm_wko <- sum_cpm_wko %>% filter(!Pathway %in% rm)
# make matrix
rownames(sum_cpm_wko) <- sum_cpm_wko$Pathway
sum_cpm_wko <- sum_cpm_wko %>% select(-Pathway)
df_m <- as.matrix(sum_cpm_wko)

#p-heatmap
library(colorspace)
colors <- diverging_hcl(100, palette = "Blue-Red 3")
pheatmap(df_m, scale = "row", cluster_cols = FALSE,cluster_rows = TRUE, cellwidth=20, 
         cellheight = 14, angle_col = 45, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

#save figure
ggsave("Fig_4b.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)

#### ---------------------------------------- Fig. 4C - Bubble plots ---------------------------------------- ####
#read the full data (counts+tax+function)
df <- read.csv("metaT_net_wc_wide.csv")
colnames(df) <- c("geneID", "taxonomy", "description", "KO", "NetTrap_150_1","NetTrap_150_2", "NetTrap_150_3",
                   "NetTrap_175_1", "NetTrap_175_2", "NetTrap_200_1", "NetTrap_200_2", "NetTrap_200_3",
                   "NetTrap_300_1", "NetTrap_300_2", "NetTrap_300_3", "WaterColumn_1_1", "WaterColumn_1_2",
                   "WaterColumn_2_1", "WaterColumn_2_2")
## NetTraps

#filter out metazoa, remove unknown KO, remove NA values
df2_net <- df %>% 
  select(-starts_with("Water")) %>%
  filter(!grepl("Metazoa", taxonomy)) %>% 
  filter(KO != "") %>%
  drop_na()

#filter rows with less than 10 sequeneces
df3_net <- df2_net %>% 
  filter(rowSums(df2_net[,5:15]) > 10)

#list of taxa
tax <- c("Chlorophyta","Chlorarachniophyta",
         "Dinophyceae",
         "Foraminifera","Polycystinea","Phaeodaria")

##loop for normalization
for(i in tax){
  tmp_counts <- subset(df3_net, grepl(i, df3_net$taxonomy))
  y <- dim(tmp_counts) [2]
  #using tmp_counts to perform edgeR normalization
  dge_obj <- DGEList(counts=tmp_counts[5:y],
                     genes=tmp_counts[1:4],
                     group=c(rep("NetTrap_150",3), rep("NetTrap_175",2), rep("NetTrap_200",3), 
                             rep("NetTrap_200",3)))
  dge_obj$samples
  data<-calcNormFactors(dge_obj, method="TMM") # TMM normalization
  data$samples # Normalized library values
  cpm_data<-cpm(data, normalized.lib.sizes=TRUE, log=FALSE) #obtain only CPM (not logged)
  cpm_data<-as.data.frame(cpm_data)
  data_CPM <- data.frame(data$genes,cpm_data)
  data_CPM$taxa <- i
  name<- paste("dfnorm", i, sep="_")
  assign(name, data_CPM)
  print("done with"); print(i)
}

all_tax_net <- rbind(dfnorm_Chlorophyta,dfnorm_Chlorarachniophyta,dfnorm_Dinophyceae,  
                 dfnorm_Polycystinea, dfnorm_Phaeodaria, dfnorm_Foraminifera)

# take average of NetTrap samples
all_tax_net <- all_tax_net %>%
  mutate(
    NetTrap_150m = rowMeans(cbind(NetTrap_150_1, NetTrap_150_2, NetTrap_150_3), na.rm = TRUE),
    NetTrap_175m = rowMeans(cbind(NetTrap_175_1, NetTrap_175_2), na.rm = TRUE),
    NetTrap_200m = rowMeans(cbind(NetTrap_200_1, NetTrap_200_2, NetTrap_200_3), na.rm = TRUE),
    NetTrap_300m = rowMeans(cbind(NetTrap_300_1, NetTrap_300_2, NetTrap_300_3), na.rm = TRUE)
  ) %>%
  select(-ends_with(c("_1","_2","_3")))

## water column

#filter out metazoa, remove unknown KO, remove NA values
df2_wc <- df %>% 
  select(-starts_with("NetTrap")) %>%
  filter(!grepl("Metazoa", taxonomy)) %>% 
  filter(KO != "") %>%
  drop_na()

#filter rows with less than 10 sequeneces
df3_wc <- df2_wc %>%
  filter(rowSums(df2_wc[,5:8]) > 10)

##loop for normalization
rm(tmp_counts, dge_obj,data,cpm_data,data_CPM,tmp2,cpm_tmp,long_avg)
for(i in tax){
  tmp_counts <- subset(df3_wc, grepl(i, df3_wc$taxonomy))
  y <- dim(tmp_counts) [2]
  #using tmp_counts to perform edgeR normalization
  dge_obj <- DGEList(counts=tmp_counts[5:y],
                     genes=tmp_counts[1:4],
                     group=c(rep("wc1",2), rep("wc2",2)))
  dge_obj$samples
  data<-calcNormFactors(dge_obj, method="TMM") # TMM normalization
  data$samples # Normalized library values
  cpm_data<-cpm(data, normalized.lib.sizes=TRUE, log=FALSE) #obtain only CPM (not logged)
  cpm_data<-as.data.frame(cpm_data)
  data_CPM <- data.frame(data$genes,cpm_data)
  data_CPM$taxa <- i
  name<- paste("dfnorm_wc", i, sep="_")
  assign(name, data_CPM)
  print("done with"); print(i)
}

# combine all watercolumn dfs
all_tax_wc <- rbind(dfnorm_wc_Chlorophyta,dfnorm_wc_Chlorarachniophyta,dfnorm_wc_Dinophyceae,  
                    dfnorm_wc_Polycystinea, dfnorm_wc_Phaeodaria, dfnorm_wc_Foraminifera)

# take average of water column samples
all_tax_wc <- all_tax_wc %>%
  mutate(
    WaterColumn_150m = rowMeans(cbind(WaterColumn_1_1, WaterColumn_1_2, WaterColumn_2_1, WaterColumn_2_2), na.rm = TRUE)
  ) %>%
  select(-ends_with(c("_1", "_2")))


## NetTrap and Water column ##
# melt the dataframes and combine
all_tax_net_melt <- melt(all_tax_net)
all_tax_wc_melt <- melt(all_tax_wc)
melt_all_tax <- rbind(all_tax_net_melt, all_tax_wc_melt)

# join the combined dataframes with ko information
all_tax_wko <- left_join(melt_all_tax, ko, by="KO")

# factor and relabel samples
all_tax_wko$Sample <- factor(all_tax_wko$variable, 
                             levels=c("NetTrap_150m","NetTrap_175m","NetTrap_200m","NetTrap_300m","WaterColumn_150m"),
                             labels=c("NetTrap 150m","NetTrap 175m","NetTrap 200m","NetTrap 300m","WaterColumn 150m"))

# factor and rename the taxa
all_tax_wko$Taxa <- factor(all_tax_wko$taxa, levels=c("Chlorophyta", "Chlorarachniophyta",
                                                      "Dinophyceae", "Polycystinea",
                                                      "Phaeodaria", "Foraminifera"),
                           labels=c("Chlorophytes", "Chlorarachniophytes",
                                    "Dinoflagellates", "Polycystines",
                                    "Phaeodarians", "Foraminiferans"))

# summarize based on Pathway and clean the dataframe for bubble plot
bubble <- all_tax_wko %>%
  select(Sample,value, Taxa, Metabolism_2) %>%
  drop_na() %>%
  group_by(Sample, Taxa, Metabolism_2) %>%
  summarise(total_cpm = sum(value)) %>%
  filter(!Metabolism_2 %in% c("Phagotrophy-other", "P metabolism"))

# factor the pathways
bubble$Metabolism_2 <- factor(bubble$Metabolism_2, levels=rev(c("Photosynthesis", "Calvin cycle", "C metabolism",
                                                                "Energy Acquisition", "Nitrogen metabolism", "Phagotrophy")),
                              labels=rev(c("Photosynthesis", "Carbon fixation", "C metabolism",
                                           "Energy Acquisition", "Nitrogen metabolism", "Phagotrophy")))

# pothway color
pathway_color <- c("lightgreen","yellow4","white","salmon","skyblue","#C4A445")
names(pathway_color) <- c("Photosynthesis", "Carbon fixation", "C metabolism",
                          "Energy Acquisition", "Nitrogen metabolism", "Phagotrophy")

# plot the combined bubble plot
fig_4c <- bubble %>%
  ggplot(aes(x=Sample, y=Metabolism_2, fill=Metabolism_2, size=total_cpm))+
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(1, 10)) +  # adjust point size range
  scale_fill_manual(values=pathway_color)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1, size=10),
        legend.position = "right")+
  guides(fill="none")+
  labs(title="", x="",y="", fill="Pathway", size="Counts per million") +
  annotate("rect", xmin = 4.5, xmax = 5.5, ymin = 0, ymax = Inf, 
           fill = 'white', alpha = 0.2) +
  annotate("rect", xmin = 0.5, xmax = 4.5, ymin = 0, ymax = Inf, 
           fill = "#C4A484", alpha = 0.2)+
  facet_wrap(~Taxa, nrow=6) 

#save figure
ggsave("Fig_4c.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)
