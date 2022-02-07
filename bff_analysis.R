setwd("/Users/stavrt/ferret")
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

# changing scaffold names to scaffold id numbers for site annt
# about 100 DMS that did not get assigned scaffold #s from merge
macau_out = read.table("bff_dfals_blood8_sig.TESTSITES", header=FALSE)
dim(macau_out_5)
dim(final_macau_out)

scaffold_key = read.table("fake_chrs_musputfur1.key", header=FALSE)
head(macau_out_5$V1)
head(scaffold_key$V1)

# data file x, key file y
merged_macau_out = merge(macau_out, scaffold_key, by.x="V1", by.y="V1", all.x=TRUE)
head(merged_macau_out)
merged_macau_out[merged_macau_out$V1=="AEYP01108960",]
final_macau_out = merged_macau_out[,c("V2.y", "V2.x")]

head(final_macau_out)
write.table(final_macau_out,file="query_MACAU_dfcount_sperm5_sig_fake_chr.txt",
            sep="\t",dec=".",row.names=FALSE,col.names=FALSE)



