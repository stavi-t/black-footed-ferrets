setwd("/Users/stavrt/ferret/dmclust")
library(dplyr)

testes_firm <- read.table("MACAU_output_firmness_testes9.assoc.clean.txt", sep = "\t", 
                          header=FALSE)
testes_firm <- setNames(testes_firm, c("id","position","n", "beta","pvalue"))

# arrange scaffolds and sort within them by genetic position
# create diff column w sequential inter-C distances
testes_firm_interC <- testes_firm %>%
  arrange(position) %>%
  arrange(id) %>%
  group_by(id) %>%
  mutate(diff = position - lag(position, default = first(position)))

# mutate by 40bp cutoff
testes_firm_interC
interC_dist <- testes_firm_interC$diff
head(interC_dist)

# calc adjusted p values
sum(p.adjust(testes_firm_interC$pvalue,method="bonferroni") < 0.05)

#testes_firm_interC %>%
#  mutate(if_else(diff <= 40, i, i+1)) %>%

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
testes_test <- make_dmclust(testes_firm_interC)
# check function vs. non function cluster_id column differences
sum(!(test$cluster_id == testes_firm_interC$cluster_id))

View(test)
write.table(testes_test,file="dmclusters-testes.txt",
            sep="\t",dec=".",row.names=FALSE,col.names=TRUE)

# make_dmclust code
group = 1  
dmcluster <- vector()

for (i in 1:length(interC_dist)) {
  if (interC_dist[i] >= 40) {
    group = group + 1
    #print(group) 
  } 
dmcluster[i] <- group
}

# add new column to dataframe
testes_firm_interC$cluster_id <- dmcluster  
testes_firm_interC$cluster_id
# test values 
length(unique(testes_firm_interC$cluster_id))
testes_firm_interC$cluster_id[testes_firm_interC$diff >= 40]

# open dataframe in visual mode
View(testes_firm_interC)
# 270k clusters


