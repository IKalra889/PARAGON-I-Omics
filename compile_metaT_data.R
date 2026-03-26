### PARAGON-I NetTrap metaT compilation ###
### Salmon counts, Eggnog annotations and EukZoo taxonomy ###
### By: Isha Kalra ###
### Last updated: 02/10/2026 ###

library(tidyr)
library(dplyr)
library(reshape2)

samples <- c("sample_01","sample_02","sample_03","sample_04","sample_05",
             "sample_07","sample_08","sample_09","sample_10","sample_11","sample_12",
             "wc1_1", "wc1_2", "wc2_1", "wc2_2") #sample 06 - NetTrap at 175m rep 3 failed so no file

# read and compile salmon outputs ({SAMPLE}.sf)
full <- NULL
for(sample in samples){
  df <- read.delim(paste(sample,"sf",sep="."),header=TRUE)
  df$Sample <- sample
  df <- df[c(1,2,5:6)] #select contig name, num reads and sample
  full <- rbind(full,df)
}

#salmon counts for all samples
salmon_count <- full
colnames(salmon_count) <- c("geneID","length" ,"count", "sample")

#eggnog annotations
net_wc_eggnog <- read.delim("net_wc_eggnog.emapper.annotations", skip = 4, header= TRUE) #first 4 lines contain code metadata
colnames(net_wc_eggnog)[1] <- "query"

#extract kegg id's from eggnog file
ko <- net_wc_eggnog %>% select(query, Description, KEGG_ko)
ko.2 <- colsplit(ko$KEGG_ko, ":", c("ko", "kegg"))
ko.3 <- colsplit(ko.2$kegg, ",", c("kegg", "extra"))
ko$kegg <- ko.3$kegg
ko <- ko %>% select(query, Description, kegg)
colnames(ko) <- c("geneID", "description", "KO")

##taxonomy file mRNA - EukZoo output
tax_net_wc <- read.delim("net_wc_contigs.txt", col.names = c("geneID", "taxonomy"))

#join all the data - salmon counts, Eggnog KO IDs and taxonomy
tax.kegg <- full_join(tax_net_wc, ko)
tax.kegg.count <- inner_join(tax.kegg, salmon_count,by="geneID")

#make data wide
metaT_wide <- pivot_wider(tax.kegg.count, names_from = sample, values_from = count) #length of contigs >=300

#save data
write.csv(metaT_wide,"metaT_net_wc_wide.csv", row.names = FALSE)
write.csv(salmon_count, "salmon_counts_net_wc_compiled.csv",row.names = FALSE)
