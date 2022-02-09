setwd("/Users/stavrt/ferret/site_annt")
library("gprofiler2")
library(dplyr)
library(stringr)
library(VennDiagram)

# g:gost tool for functional enrichment analysis
# from Raudvere 2019 update: maps a user-provided list of genes to known functional information 
# sources and detects statistically significantly enriched biological processes, pathways, 
# regulatory motifs and protein complexes. Main source of information about genes, 
# identifier types, GO terms and associations are taken from ENSEMBL.
# Inputs for this function are a vector of gene identifiers, 
# query, and the name of the corresponding organism which is 
# constructed by concatenating the first letter of the genus name and the specific epithet

# create input files from site-annt with gene IDs
# remove .1 ext from gene IDs, rows containing "inter", unannotated rows and all duplicate ID #s
df_blood_als <- read.table("output_ann_MACAU_dfals_blood8_sig_fake_chr.txt", sep = "\t", fill=TRUE)
df_sperm_als <- read.table("output_ann_MACAU_dfals_sperm5_sig_fake_chr.txt", sep = "\t", fill=TRUE)
df_testes_als <- read.table("output_ann_MACAU_dfals_testes9_sig_fake_chr.txt", sep = "\t", fill=TRUE)
df_sperm_count <- read.table("output_ann_MACAU_dfcount_sperm5_sig_fake_chr.txt", sep = "\t", fill=TRUE)
df_testes_firm <- read.table("output_ann_MACAU_dffirm_testes9_sig_fake_chr.txt", sep = "\t", fill=TRUE)

objs <- list(
  df_blood_als, 
  df_sperm_als, 
  df_testes_als , 
  df_sperm_count,
  df_testes_firm)

for(j in 1:length(objs)){
  
  result <- objs[[j]] %>%
    filter(V3 != "inter") %>%
    na.omit() %>%
    mutate(V5 = substr(V3, 0,18)) %>%
    # output gene IDs as character vectors using pull instead of select
    pull(V5) %>%
    unique()
  
  output[[j]] <- result
}

output[[1]]
class(output[[5]])

# blood ALS 10473 gene IDs
# sperm ALS 12053 gene IDs
# testes ALS 10249 gene IDs
# sperm COUNT 11628 gene IDs
# testes FIRMNESS 5150 gene IDs

blood_als_gene_IDS <- output[[1]]
sperm_als_gene_IDS <- output[[2]]
testes_als_gene_IDS <- output[[3]]
sperm_count_gene_IDS <- output[[4]]
testes_firm_gene_IDS <- output[[5]]

# how many distinct ENSEMBL IDs are there in original file?
n_before <- length(unique(df_blood_als[!grepl("inter", df_blood_als$V3),]$V3))

# remove ".1" ext from gene IDs, rows containing "inter" only, all N/A rows and duplicate values
df_blood_gene_IDS <- df_blood_als %>%
   filter(V3 != "inter") %>%
   na.omit() %>%
   mutate(V5 = substr(V3, 0,18)) %>%
   select(V5) %>%
   distinct()

class(df_blood_gene_IDS)
# dataframe

n_after <- nrow(df_blood_gene_IDS)
# sites lost due to no annotation from site annotation script
n_after-n_before

# gene IDs = query files
?gprofiler2::gost
# running gprofiler2
# output: list of result (data.frame) & metadata (list for manhattan plot and query metadata)
blood_als_gprof <- gost(query = c(blood_als_gene_IDS),
                   organism = "mpfuro",
                   domain_scope = "annotated",
                   correction_method="fdr",
                   user_threshold = 0.05,
                   numeric_ns = "ENTREZGENE_ACC",
                   sources = c("GO:BP", "GO:MF", "GO:CC"), 
                   as_short_link = FALSE)
blood_als_gprof 
blood_als_go_terms <- blood_als_gprof$result$term_id

sperm_als_gprof <- gost(query = c(sperm_als_gene_IDS),
                   organism = "mpfuro",
                   domain_scope = "annotated",
                   correction_method="fdr",
                   user_threshold = 0.05,
                   numeric_ns = "ENTREZGENE_ACC",
                   sources = c("GO:BP", "GO:MF", "GO:CC"), 
                   as_short_link = FALSE)

sperm_als_gprof 
sperm_als_go_terms <- sperm_als_gprof$result$term_id
 
testes_als_gprof <- gost(query = c(testes_als_gene_IDS),
                    organism = "mpfuro",
                    domain_scope = "annotated",
                    correction_method="fdr",
                    user_threshold = 0.05,
                    numeric_ns = "ENTREZGENE_ACC",
                    sources = c("GO:BP", "GO:MF", "GO:CC"), 
                    as_short_link = FALSE)

testes_als_gprof 
testes_als_go_terms <- testes_als_gprof$result$term_id

sperm_count_gprof <- gost(query = c(sperm_count_gene_IDS),
                     organism = "mpfuro",
                     domain_scope = "annotated",
                     correction_method="fdr",
                     user_threshold = 0.05,
                     numeric_ns = "ENTREZGENE_ACC",
                     sources = c("GO:BP", "GO:MF", "GO:CC"), 
                     as_short_link = FALSE)

sperm_count_gprof 
sperm_count_go_terms <- sperm_count_gprof$result$term_id

testes_firm_gprof <- gost(query = c(testes_firm_gene_IDS),
                     organism = "mpfuro",
                     domain_scope = "annotated",
                     correction_method="fdr",
                     user_threshold = 0.05,
                     numeric_ns = "ENTREZGENE_ACC",
                     sources = c("GO:BP", "GO:MF", "GO:CC"), 
                     as_short_link = FALSE)

testes_firm_gprof
testes_firm_go_terms <- testes_firm_gprof$result$term_id
 
# venn diagrams of GO term overlap - mpfuro (DF) genome
venn.diagram(x=list(blood_als_go_terms,sperm_als_go_terms,testes_als_go_terms),
             category.names=c("Blood","Sperm","Testes"),
             filename="GO_MACAU_B_S_T_ALS_site_annt.png",
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

venn.diagram(x=list(sperm_als_go_terms,sperm_count_go_terms),
             category.names=c("Sperm ALS","Sperm Count"),
             filename="GO_MACAU_SPERM_site_annt.png",
             output=TRUE,
             # circles
             lwd=2,
             lty="blank",
             fill=c("lightsteelblue2","lightsteelblue2"),
             # numbers
             cex=1.5,
             fonface="bold",
             # set names
             cat.cex=1.5,
             cat.fontface = "bold",
             cat.default.pos="outer")

venn.diagram(x=list(testes_als_go_terms,testes_firm_go_terms),
             category.names=c("Testes ALS","Testes Firm"),
             filename="GO_MACAU_TESTES_site_annt.png",
             output=TRUE,
             # circles
             lwd=2,
             lty="blank",
             fill=c("navajowhite3","navajowhite3"),
             # numbers
             cex=1.5,
             fonface="bold",
             # set names
             cat.cex=1.5,
             cat.fontface = "bold",
             cat.just=list(c(0,0) , c(1,0)),
             cat.default.pos="outer")

