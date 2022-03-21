setwd("/Users/stavrt/ferret/")
library(VennDiagram)
library(wesanderson)

# read 5 macau outputs
dfals_B = read.table("bff_avglittersize_blood8_df.assoc_sorted.txt", header=TRUE)
dfals_S = read.table("bff_avglittersize_sperm5_df.assoc_sorted.txt", header=TRUE)
dfals_T = read.table("bff_avglittersize_testes9_df.assoc_sorted.txt", header=TRUE)
dfT = read.table("bff_firmness_testes9_df.assoc_sorted.txt", header=TRUE)
dfS = read.table("bff_spermcountperml_sperm5_df.assoc_sorted.txt", header=TRUE)

# subset dfs by beta and p-values
# cytosines w p-values in 1% x beta values in 1% or 99% distribution
dfals_B_sig <- subset(dfals_B, pvalue <= quantile(pvalue, 0.01) & (beta <= quantile(beta, 0.01) | beta >= quantile(beta, 0.99)))
dfals_S_sig <- subset(dfals_S, pvalue <= quantile(pvalue, 0.01) & (beta <= quantile(beta, 0.01) | beta >= quantile(beta, 0.99)))
dfals_T_sig <- subset(dfals_T, pvalue <= quantile(pvalue, 0.01) & (beta <= quantile(beta, 0.01) | beta >= quantile(beta, 0.99)))
dfT_sig <- subset(dfT, pvalue <= quantile(pvalue, 0.01) & (beta <= quantile(beta, 0.01) | beta >= quantile(beta, 0.99)))
dfS_sig <- subset(dfS, pvalue <= quantile(pvalue, 0.01) & (beta <= quantile(beta, 0.01) | beta >= quantile(beta, 0.99)))

# write subset dfs as csv
write.table(dfals_B_sig,file="bff_dfals_blood8_sig.txt",sep="\t",dec=".",row.names=FALSE,col.names=TRUE)
write.table(dfals_S_sig,file="bff_dfals_sperm5_sig.txt",sep="\t",dec=".",row.names=FALSE,col.names=TRUE)
write.table(dfals_T_sig,file="bff_dfals_testes9_sig.txt",sep="\t",dec=".",row.names=FALSE,col.names=TRUE)
write.table(dfT_sig,file="bff_dffirm_testes9_sig.txt",sep="\t",dec=".",row.names=FALSE,col.names=TRUE)
write.table(dfS_sig,file="bff_dfcount_sperm5_sig.txt",sep="\t",dec=".",row.names=FALSE,col.names=TRUE)

# VENN DIAGRAMS: overlap of significant DMS from 5 MACAU analyses
vector_als_blood <- as.vector(dfals_B_sig$id)
vector_als_sperm <- as.vector(dfals_S_sig$id)
vector_als_testes <- as.vector(dfals_T_sig$id)
vector_count_sperm <- as.vector(dfS_sig$id)
vector_firm_testes <- as.vector(dfT_sig)

# als all - 1 overlapping
venn.diagram(x=list(vector_als_blood,vector_als_sperm,vector_als_testes),
             category.names=c("Blood","Sperm","Testes"),
             filename="MACAU_VennDiagram_ALS_OVERLAP_DMS.png",
             output=TRUE,
             # circles
             lwd=2,
             lty="blank",
             fill=c("brown4","lightsteelblue2","navajowhite3"),
             # numbers
             cex=1.5,
             fonface="bold",
             # set names
             cat.cex=2,
             cat.fontface = "bold",
             cat.default.pos="outer")

# sperm all 
venn.diagram(x=list(vector_als_sperm,vector_count_sperm),
             category.names=c("Sperm: ALS","Sperm: Count"),
             filename="MACAU_VennDiagram_SPERM_OVERLAP_DMS.png",
             output=TRUE,
             # circles
             lwd=2,
             lty="blank",
             fill=c("brown4","lightsteelblue2"),
             # numbers
             cex=1.5,
             fonface="bold",
             # set names
             cat.cex=0.5,
             cat.fontface = "bold",
             cat.default.pos="text")

# testes all - NO OVERLAP
venn.diagram(x=list(vector_als_testes,vector_firm_testes ),
             category.names=c("Testes: ALS","Testes: Firmness"),
             filename="MACAU_VennDiagram_TESTES_OVERLAP_DMS.png",
             output=TRUE,
             # circles
             lwd=2,
             lty="blank",
             fill=c("lightsteelblue2","navajowhite3"),
             # numbers
             cex=1.5,
             fonface="bold",
             # set names
             cat.cex=0.5,
             cat.fontface = "bold",
             cat.default.pos="outer")

# ***changing scaffold names to scaffold id numbers for site annt***
# omits 'N/A' DMS that don't get assigned scaffold #s from merge 
## (one way merge - omits many sites! make sure you want to do this!) ##
setwd("/Users/stavrt/ferret/site_annt/MACAU-raw/")
test_sites = read.table("MACAU_output_avglittersize_blood8.assoc.clean.TESTSITES", 
                        header=FALSE)
dim(test_sites)
dim(final_macau_out)

scaffold_key = read.table("/Users/stavrt/ferret/fake_chrs_musputfur1.key", header=FALSE)
head(test_sites$V1)
head(scaffold_key$V1)

# data file x (V1, V2), key file y (V1, V2)
merge_scaffolds <- function(test_sites, scaffold_key) {
  macau_file = read.table(test_sites, header=FALSE)
  names(scaffold_key)
  merged_macau_out = merge(macau_file, scaffold_key, by.x="V1", by.y="V1", all.x=FALSE)
  final_macau_out = merged_macau_out[,c("V2.y", "V2.x")]
  
  out_file = paste0("query_", test_sites)
  out_file = gsub(".TESTSITES", ".txt", out_file)
  write.table(final_macau_out,file=paste0("/Users/stavrt/ferret/site_annt/MACAU-raw/", out_file),
                                          sep="\t",dec=".",row.names=FALSE,col.names=FALSE)
}

test_site_list <-
  list.files(pattern=".TESTSITES")
lapply(test_site_list, merge_scaffolds, scaffold_key=scaffold_key)

merge_scaffolds(test_site_list, scaffold_key)

# merge_scaffolds("MACAU_output_avglittersize_blood8.assoc.clean.TESTSITES", scaffold_key)
# merge_scaffolds("MACAU_output_avglittersize_sperm5.assoc.clean.TESTSITES", scaffold_key)
# merge_scaffolds("MACAU_output_avglittersize_testes9.assoc.clean.TESTSITES", scaffold_key)
# merge_scaffolds("MACAU_output_firmness_testes9.assoc.clean.TESTSITES", scaffold_key)
# merge_scaffolds("MACAU_output_spermcount_sperm5.assoc.clean.TESTSITES", scaffold_key)

# checking for replicate values
#merged_macau_out[merged_macau_out$V1=="AEYP01108960",]
  

