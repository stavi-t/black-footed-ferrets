setwd("/Users/stavrt/ferret/dmclust")
library(dplyr)
library(ggplot2)
# read MACAU files - unsorted
testes_firm <- read.table("MACAU_output_firmness_testes9.assoc.clean.txt", sep = "\t", 
                          header=FALSE)
testes_firm <- setNames(testes_firm, c("id","position","n", "beta","pvalue"))
##
sperm_count <- read.table("MACAU_output_spermcount_sperm5.assoc.clean.txt", sep = "\t",
                          header=FALSE)
sperm_count <- setNames(sperm_count, c("id","position","n", "beta","pvalue"))
##
blood_litter <- read.table("MACAU_output_avglittersize_blood8.assoc.clean.txt", sep = "\t",
                           header=FALSE)
blood_litter <- setNames(blood_litter, c("id","position","n", "beta","pvalue"))
##
sperm_litter <- read.table("MACAU_output_avglittersize_sperm5.assoc.clean.txt", sep = "\t",
                           header=FALSE)
sperm_litter <- setNames(sperm_litter, c("id","position","n", "beta","pvalue"))
##
testes_litter <- read.table("MACAU_output_avglittersize_testes9.assoc.clean.txt", sep = "\t",
                            header=FALSE)
testes_litter <- setNames(testes_litter, c("id","position","n", "beta","pvalue"))
##

head(sperm_count)

# arrange scaffolds and sort within each scaffold by genomic position
# create diff column w sequential inter-C distances
testes_firm_interC <- testes_firm %>%
  arrange(position) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(diff = position - lag(position, default = first(position)))

sperm_count_interC <- sperm_count %>%
  arrange(position) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(diff = position - lag(position, default = first(position)))

blood_litter_interC <- blood_litter %>%
  arrange(position) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(diff = position - lag(position, default = first(position)))

sperm_litter_interC <- sperm_litter %>%
  arrange(position) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(diff = position - lag(position, default = first(position)))

testes_litter_interC <- testes_litter %>%
  arrange(position) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(diff = position - lag(position, default = first(position)))

# mutate by 40bp cutoff
testes_firm_interC
head(sperm_count_interC)
interC_dist <- testes_firm_interC$diff
head(interC_dist)

### adjust p values in place
#sum(p.adjust(testes_firm_interC$pvalue,method="bonferroni") < 0.05)
#p <- p.adjust(testes_firm_interC$pvalue,method="bonferroni") < 0.05
###
make_dmclust <- function(df, threshold=40, diff=6) {
  df <- as.data.frame(df)
  group = 1  
  cluster_id <- vector()
  intercyt_dist <- df[,diff]
  for (i in 1:length(intercyt_dist)) {
    if (intercyt_dist[i] >= threshold) {
      group = group + 1
    } 
    cluster_id[i] <- group
  }
  df <- cbind(df, cluster_id)
}

# run function
testes_clust <- make_dmclust(testes_firm_interC)
# 250-270k clusters for all
sperm_clust <- make_dmclust(sperm_count_interC)
blood_clust <- make_dmclust(blood_litter_interC)
sperm_als_clust <- make_dmclust(sperm_litter_interC)
testes_als_clust <- make_dmclust(testes_litter_interC)

length(unique(testes_test$cluster_id))

# check function vs. non function cluster_id column differences
sum(!(test$cluster_id == testes_firm_interC$cluster_id))

# calc bonferroni corr for each file
0.05/count(testes_clust)
0.05/count(sperm_count)
0.05/count(blood_litter)
View(blood_clust)

# filter by bonferroni corr pval
sig_testes <- testes_clust %>%
  filter(testes_test$pvalue < 6.564446e-09)

sig_sperm <- sperm_clust %>%
  filter(sperm_clust$pvalue < 5.266355e-09)

sig_blood <- blood_clust %>%
  filter(blood_clust$pvalue < 5.303714e-09)

sig_sperm_als <- sperm_als_clust %>%
  filter(sperm_als_clust$pvalue < 5.595129e-09)

sig_testes_als <- testes_als_clust %>%
  filter(testes_als_clust$pvalue < 5.483368e-09)
View(sig_blood)  
# no significant results for sperm count...OR sperm litter size!!!
# 4 significant results for blood! predictors????
# 11 significant results for testes litter size -> one large cluster!

min(sperm_clust$pvalue)
max(sperm_clust$pvalue)

View(sig_blood)

# plot pvalues
ggplot(sig_testes, aes(pvalue)) + 
  geom_histogram(color="black", fill="gray", bins=70) +
  xlab("MACAU Bonferroni p values")

# make_dmclust function
group = 1  
dmcluster <- vector()

for (i in 1:length(interC_dist)) {
  if (interC_dist[i] >= 40) {
    group = group + 1
    #print(group) 
  } 
dmcluster[i] <- group
}

# add new column to df
testes_firm_interC$cluster_id <- dmcluster  
testes_firm_interC$cluster_id
# test values are correct
length(unique(testes_firm_interC$cluster_id))
testes_firm_interC$cluster_id[testes_firm_interC$diff >= 40]

# open df in visual mode
View(testes_firm_interC)
View(sig)

# write df to file
write.table(testes_test,file="dmclusters-testes-sig.txt",
            sep="\t",dec=".",row.names=FALSE,col.names=TRUE)




