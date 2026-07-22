### PARAGON-I Daily PIT plot - Figure 3 ###
### POC, Ordination, Heat map and bubble plot###
### By: Isha Kalra ###
### Last updated: 02/24/2026 ###

library(tidyverse)
library(reshape2)
library(vegan)
library(pheatmap)
library(edgeR)

#### ------------------------------------- Fig 3A - POC plot ---------------------------------------------------- ####
cnp <- read.csv("dailyPIT_CNP.csv")

cnp$sample <- as.factor(c("Daily PIT day1","Daily PIT day2",
                "Daily PIT day3","Daily PIT day4","Daily PIT day5","Daily PIT day6","Daily PIT day7",
                "Daily PIT day8","Daily PIT day9","Daily PIT day10"))

cnp_df <- cnp %>% 
  filter(!Day == "Day-04") %>% 
  select(POC, sample)

cnp_df$sample <- factor(cnp_df$sample, levels = c("Daily PIT day1","Daily PIT day2",
                                           "Daily PIT day3","Daily PIT day5","Daily PIT day6","Daily PIT day7",
                                           "Daily PIT day8","Daily PIT day9","Daily PIT day10"),
                        labels = c("Daily PIT day1","Daily PIT day2",
                                   "Daily PIT day3","Daily PIT day5","Daily PIT day6","Daily PIT day7",
                                   "Daily PIT day8","Daily PIT day9","Daily PIT day10"))

cnp_plot <- ggplot(cnp_df, aes(x=sample, y=POC)) +
  geom_point(size = 4)+
  geom_line(group="sample", linewidth = 1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45,size=12, hjust = 1), 
        axis.text.y = element_text(size=11), axis.title.x = element_blank(), axis.title.y = element_text(size=12))+
  labs(y = expression(POC~(mg~m^{-2}~day^{-1})))

ggsave("Fig_3a.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)

##### ------------------------------------- Fig. 3B - PCA plot ------------------------------------------------- #####
## read the full raw count data (counts+tax+function)
df <- read.csv("metaT_dpit_wc_wide.csv")

## remove fungi , metazoa, transcripts with no KO id and daily PIT day 4 sample
tmp <- df %>%
  filter(!grepl("Fungi|Metazoa", taxonomy)) %>%
  filter(KO != "") %>%
  select(-X , -dpit_04)

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
  mutate(WaterColumn_day1 = rowMeans(cbind(wc1_1, wc1_2, na.rm = TRUE))) %>%
  mutate(WaterColumn_day5 = rowMeans(cbind(wc2_1, wc2_2, na.rm = TRUE))) %>%
  select(-starts_with("wc"))

#make rownames as KO
rownames(tmpSum) <- tmpSum$KO
tmpSum$KO <- NULL

# transcform to make rows samples, and columns as KO
tmpSum <- as.data.frame(t(tmpSum))
tmpSum <- tmpSum[, colSums(tmpSum) > 0]

library(compositions)
# calculate proportions for clr transformations (a.k.a relative abundance)
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
colz <- colsplit(rownames(pca.plot), "_", c("collection", "day"))
# factor data
colz$collection <- factor(colz$collection, levels=c("dpit", "WaterColumn"),
                          labels=c("Daily PIT", "WaterColumn"))
colz$day <- factor(colz$day, levels=c("01", "02", "03","05","06","07","08","09","10","day1","day5"),
                   labels=c("1", "2", "3","5","6","7","8","9","10","1","5"))

# add to pca df
pca.plot$collection <- colz$collection
pca.plot$day <- colz$day

pca_plot <- ggplot(pca.plot,aes(PC1,PC2,shape = collection, fill=collection,
                    label=day))+
  geom_point(size = 6, aes(fill=collection))+
  geom_text(color = "white", size = 3)+
  scale_shape_manual(values = c(21, 22, 23, 24, 25))+
  labs(title = "", x="PC1: 28.77% variance", y="PC2: 14.56% variance",
       fill=NULL, shape=NULL)+
  theme_bw()+
  theme(axis.text = element_blank(), legend.position = "right")

ggsave("Fig_3b.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)

#### ------------------------------------------------- Fig. 3C - Heatmap ------------------------------------------------- ####
#read the full normalized data (counts+tax+function) -no metazoa, only dpit
df <- read.csv("metaT_dpit_only_normalized_filtered.csv") #277486 rows

#remove fungi
df2 <- df %>% filter(!grepl("Fungi", taxonomy))

#read custom ko list
ko <- read.csv("ko_list_full.csv")

#join the normalized data with the custom KO categories
CPM_wKO <- inner_join(df2, ko, by="KO")

#select samples and pathway
CPM_wKO <- CPM_wKO %>% select(starts_with("dpit"), Pathway)

# summarize by pathway
sum_cpm_wko <- CPM_wKO %>%
  group_by(Pathway) %>%
  summarise_all(sum) %>%
  ungroup() %>%
  as.data.frame()

colnames(sum_cpm_wko) <- c("Pathway", "Daily PIT day1","Daily PIT day2",
                           "Daily PIT day3", "Daily PIT day5","Daily PIT day6","Daily PIT day7",
                           "Daily PIT day8","Daily PIT day9","Daily PIT day10")

rm <- c("else", "por", "else", "other", "additional breakdown" )
sum_cpm_wko <- sum_cpm_wko %>% filter(!Pathway %in% rm)

# make matrix
rownames(sum_cpm_wko) <- sum_cpm_wko$Pathway
sum_cpm_wko <- sum_cpm_wko %>% select(-Pathway)
df_m <- as.matrix(sum_cpm_wko)

#p-heatmap
library(colorspace)
colors <- diverging_hcl(100, palette = "Blue-Red 3")
pheatmap(df_m, scale = "row", cluster_cols = FALSE,
         cluster_rows = TRUE, cellwidth=14, 
         cellheight = 14, angle_col = 45, 
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

pheatmap(df_m, scale = "column", cluster_cols = FALSE,cluster_rows = TRUE, cellwidth=14, 
         cellheight = 14, angle_col = 45, color = colorRampPalette(c("navy", "white", "firebrick3"))(100))

ggsave("Fig_3c.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)


#### ------------------------------------------------- Fig. 3D - Bubble plot ------------------------------------------------- ####

#read the full data (counts+tax+function)
df <- read.csv("metaT_dpit_wc_wide.csv")

## Daily PITs

#filter out metazoa, remove unknown KO, remove NA values
df2_dpit <- df %>% 
  select(-starts_with("wc")) %>%
  filter(!grepl("Metazoa", taxonomy)) %>% 
  select(-X) %>%
  filter(KO != "") %>%
  drop_na()

#filter rows with less than 10 sequeneces
df3_dpit <- df2_dpit %>% 
  select(-dpit_04) %>%
  filter(rowSums(df2_dpit[,5:13]) > 10)

#list of taxa
tax <- c("Chlorophyta","Chlorarachniophyta",
         "Dinophyceae",
         "Foraminifera","Polycystinea","Phaeodaria")

##loop for normalization
for(i in tax){
  tmp_counts <- subset(df3_dpit, grepl(i, df3_dpit$taxonomy))
  y <- dim(tmp_counts) [2]
  #using tmp_counts to perform edgeR normalization
  dge_obj <- DGEList(counts=tmp_counts[5:y],
                     genes=tmp_counts[1:4],
                     group=c(rep("dpit_01",1), rep("dpit_02",1), rep("dpit_03",1), 
                            rep("dpit_05",1), rep("dpit_06",1),
                             rep("dpit_07",1), rep("dpit_08",1), rep("dpit_09",1),
                             rep("dpit_10",1)))
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

all_tax <- rbind(dfnorm_Chlorophyta,dfnorm_Chlorarachniophyta,dfnorm_Dinophyceae,  
                 dfnorm_Polycystinea, dfnorm_Phaeodaria, dfnorm_Foraminifera)


## water column

#filter out metazoa, remove unknown KO, remove NA values
df2_wc <- df %>% 
  select(-starts_with("dpit")) %>%
  filter(!grepl("Metazoa", taxonomy)) %>% 
  select(-X) %>%
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
    WaterColumn_day1 = rowMeans(cbind(wc1_1, wc1_2), na.rm = TRUE),
    WaterColumn_day5 = rowMeans(cbind(wc2_1, wc2_2), na.rm = TRUE)
  ) %>%
  select(-starts_with("wc"))

## DPIT and Water column ##
# melt the dataframes and combine
all_tax_melt <- melt(all_tax)
all_tax_wc_melt <- melt(all_tax_wc)
melt_all_tax <- rbind(all_tax_melt, all_tax_wc_melt)

# join the combined dataframes with ko information
all_tax_wko <- left_join(melt_all_tax, ko, by="KO")

# factor and relabel samples
all_tax_wko$Sample <- factor(all_tax_wko$variable, 
                             levels=c("dpit_01", "dpit_02", "dpit_03",
                                      "dpit_05","dpit_06","dpit_07","dpit_08","dpit_09", "dpit_10", "WaterColumn_day1", "WaterColumn_day5"),
                             labels=c("Daily PIT day1","Daily PIT day2",
                                        "Daily PIT day3","Daily PIT day5","Daily PIT day6","Daily PIT day7",
                                        "Daily PIT day8","Daily PIT day9","Daily PIT day10", "WaterColumn day1", "WaterColumn day5"))

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
bubble %>%
  ggplot(aes(x=Sample, y=Metabolism_2, fill=Metabolism_2, size=total_cpm))+
  geom_point(shape = 21, color = "black") +
  scale_size(range = c(1, 10)) +  # adjust point size range
  scale_fill_manual(values=pathway_color)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1, size=10),
        legend.position = "right")+
  guides(fill="none")+
  labs(title="", x="",y="", fill="Pathway", size="Counts per million") +
  annotate("rect", xmin = 9.5, xmax = 12.5, ymin = 0, ymax = Inf, 
           fill = 'white', alpha = 0.2) +
  annotate("rect", xmin = 0.5, xmax = 9.5, ymin = 0, ymax = Inf, 
           fill = "#C4A484", alpha = 0.2)+
  facet_wrap(~Taxa, nrow=6) 

ggsave("Fig_3d.pdf", width = 6.0, height = 8.0, units = "in", dpi = 600)




