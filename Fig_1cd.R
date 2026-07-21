### PARAGON-I manuscript Figure 1 c&d ###
### 18S rRNA diversity analyses ###
### Daily-PITS, Net-trap and Water-column samples ###
### By: Isha Kalra ###
### Last Updated: 07/21/2026 ###

library(tidyverse)
library(ggplot2)
library(microeco)
library(reshape2)
library(ggpubr)

#set seed
set.seed(123)
# set the plotting background
theme_set(theme_bw())

##load files
#asv table - without single and doubletons
asv <- read.csv("asv_nodouble.csv") #reading asv table and remove blank
asv <- asv %>%
  select(-Blank)

#taxonomy - without single and doubletons
tax <- read.csv("taxonomy_pr2.v5_nodouble.csv")
colnames(tax)[1] <- "otu_id"

#split taxonomy into levels
tax1 <- colsplit(tax$Taxon, ";", c("Domain", "Supergroup", "Division","Subdivision", "Class", "Order",
                                   "Family", "Genus", "Species"))
tax <- cbind(tax, tax1) %>% select(-Taxon)

#metadata
sample_info <- read.csv("metadata_18s_paragon1.csv") %>% 
  filter(!Sample == "Blank")

#making all column1 as rownames
row.names(asv) <- NULL
asv <- asv %>%
  tibble::column_to_rownames("otu_id")

tax <- tax %>% 
  tibble::column_to_rownames("otu_id")

sample_info <- sample_info %>%
  tibble::column_to_rownames("Sample")

#create microtable object
euks <- microtable$new(sample_table = sample_info, otu_table = asv, 
                       tax_table = tax)
euks$tidy_dataset() #8425 ASVs

##Only protist - remove metazoa and fungi ASVs
protist <- clone(euks)
#remove metazoa
protist$tax_table <- subset(protist$tax_table, Subdivision != "Metazoa")
protist$tidy_dataset()
#remove fingi
protist$tax_table <- subset(protist$tax_table, Subdivision != "Fungi")
protist$tidy_dataset() #7515 ASVs

### --------------------------------------- Figure 1c - Ordination Plots -------------------------------------- ###

# first normalize data using total sum scaling
df <- trans_norm$new(dataset = protist)
#total sum scaling method for normalisation
protist.tss <- df$norm(method = "TSS")

#calculate beta diversity
protist.tss$cal_betadiv()
t1 <- trans_beta$new(dataset = protist.tss, group = "Group", measure = "bray")

# NMDS calculation
t1$cal_ordination(method = "NMDS")

# plot NMDS with stress value
fig_1c <- t1$plot_ordination(plot_color = "SampleType",plot_shape = "Depth",
                                   plot_type = c("point", "ellipse"))

### ---------------------------------------- Figure 1d - Shannon Diversity --------------------------------------- ###
t2 <- trans_alpha$new(dataset = protist.tss, group = "Depth", by_group = "SampleType")
# return t1$data_stat
head(t2$data_stat)

# plot Shannon diversity
fig_1d <- t2$plot_alpha(measure = "Shannon", add_sig=FALSE, add = "jitter", order_x_mean = TRUE)

# calculate significant differences between sample types
t3 <- trans_alpha$new(dataset = protist.tss, group = "SampleType")
head(t3$data_stat)
t3$cal_diff(method = "KW_dunn")
