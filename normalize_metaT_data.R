### PARAGON-I NetTrap metaT normalization ###
### EdgeR TMM normalization ###
### By: Isha Kalra ###
### Last updated: 02/11/2026 ###

library(edgeR)
library(tidyr)
library(dplyr)
library(reshape2)

setwd("~/Desktop/Caron_lab_research/SCOPE/PARAGON-I/manuscript/analyses/compiled_data/")

#read the full data (counts+tax+function)
df <- read.csv("metaT_net_wc_wide.csv")
#df <- metaT_wide

#remove metazoa sequences, remove genes that don't have a kegg ID, remove NA values
df2 <- df %>% 
  filter(!grepl("Metazoa", taxonomy)) %>% 
  filter(KO != "")

#remove low abundance transcripts - filter rows with less than 10 sequeneces
df3 <- df2 %>% filter(rowSums(df2[,5:19]) > 10) #584137
colnames(df3) <- c("geneID", "taxonomy", "description", "KO", "NetTrap_150_1","NetTrap_150_2", "NetTrap_150_3",
                   "NetTrap_175_1", "NetTrap_175_2", "NetTrap_200_1", "NetTrap_200_2", "NetTrap_200_3",
                   "NetTrap_300_1", "NetTrap_300_2", "NetTrap_300_3", "WaterColumn_1_1", "WaterColumn_1_2",
                   "WaterColumn_2_1", "WaterColumn_2_2")

# Set up edgeR list. Counts are salmon sample counts, genes are identifiers (i.e. geneID, KO, and taxonomy), groups are you telling edgeR which groups exist in your dataset (i.e. sample types)
dge <- DGEList(counts=df3[5:19],genes=df3[1:4], 
               group=c(rep("NetTrap_150",3), rep("NetTrap_175",2), rep("NetTrap_200",3), 
                       rep("NetTrap_300",3), rep("WaterColumn_1",2), rep("WaterColumn_2",2)))

#check the groupings
dge$samples

#normalise data based on TMM method
data <- calcNormFactors(dge, method="TMM") # TMM normalization
data$samples # Normalized library values
cpm_data <- cpm(data, normalized.lib.sizes=TRUE, log=FALSE) #obtain only CPM (not logged)
cpm_data <- as.data.frame(cpm_data)
data_CPM <- data.frame(data$genes,cpm_data)

#save file
write.csv(data_CPM,"metaT_net_wc_normalized.csv", row.names = FALSE)

##--------------------------------------------------------------------###

### only NetTrap samples ###

df3_net <- df3 %>% select(1:4, starts_with("NetTrap"))

#set up EdgeR list
dge.net <- DGEList(counts=df3_net[5:15],genes=df3_net[1:4], 
               group=c(rep("NetTrap_150",3), rep("NetTrap_175",2), rep("NetTrap_200",3), 
                       rep("NetTrap_300",3)))

#normalise
data.net <- calcNormFactors(dge.net, method="TMM")
data.net$samples 
cpm.net <- cpm(data.net, normalized.lib.sizes = TRUE, log = FALSE)
cpm.net <- as.data.frame(cpm.net)
data.net.cpm <- data.frame(data.net$genes, cpm.net)

#save file
write.csv(data.net.cpm,"metaT_net_only_normalized.csv", row.names = FALSE)
