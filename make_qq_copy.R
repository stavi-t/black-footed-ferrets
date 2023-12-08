setwd("/Users/stavrt/ferret/")
library(dplyr)
library(lattice)

# read 5 macau outputs
dfals_B = read.table("bff_avglittersize_blood8_df.assoc_sorted.txt", header=TRUE)
dfals_S = read.table("bff_avglittersize_sperm5_df.assoc_sorted.txt", header=TRUE)
dfals_T = read.table("bff_avglittersize_testes9_df.assoc_sorted.txt", header=TRUE)
dfT = read.table("bff_firmness_testes9_df.assoc_sorted.txt", header=TRUE)
dfS = read.table("bff_spermcountperml_sperm5_df.assoc_sorted.txt", header=TRUE)

head(dfS)
#since we're usually most interested in really small p-values, 
# we generally transform the p-values by -log10 so that the smallest values 
# near zero become the larger values and are thus easier to see.

# subset MACAU p values
pvals_dfals_b <- dfals_B %>%
  select(c(pvalue))

pvals_dfals_s <- dfals_S %>%
  select(c(pvalue))

pvals_dfals_t <- dfals_T %>%
  select(c(pvalue))

pvals_dfT <- dfT %>%
  select(c(pvalue))

pvals_dfS <- dfS %>%
  select(c(pvalue))


head(pvals_dfS)

# qqplots for pvals -> need to be log transformed, using a uniform distribution
qqplot <- qqmath(~-log10(pvals_dfS$pvalue),
       distribution=function(x){-log10(qunif(1-x))},
       pch = 19,
       xlab = "Expected -log[10](p-value)", ylab = "Observed -log[10](p-value)",
       ylim = c(0, 12),
       xlim = c(0, 12),
       aspect = 1,
       panel = function(x, y, ...) {
         # Draw the QQ plot points
         panel.qqmath(x, y, ...)
         
         # Add a diagonal trendline y = x
         panel.abline(a = 0, b = 1, col = "red", lty = 3)
       })

print(qqplot)

# examine data type
head(pvals_dfals_b)
class(pvals_dfals_b$pvalue)


