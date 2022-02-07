# plotting SITE ANNT results of significant DMS in DF genome
setwd("/Users/stavrt/ferret/site_annt")
library(devtools)
library(ggplot2)
library(patchwork)
library(dplyr)
# load theme
names(wes_palettes)
wes_palette("FantasticFox1")

#clrs <- wes_palette(exon, intergenic, intron, promoter)
clrs <- c("#DD8D29", "#E2D200", "#46ACC8", "#B40F20")

# clean site annt outputs - select last entry in df only and replace "prom error" lines with prom
SA_blood_als_file_name <- "output_ann_MACAU_dfals_blood8_sig_fake_chr.txt"
SA_sperm_als_file_name <- "output_ann_MACAU_dfals_sperm5_sig_fake_chr.txt"
SA_testes_als_file_name <- "output_ann_MACAU_dfals_testes9_sig_fake_chr.txt"
SA_sperm_count_file_name <- "output_ann_MACAU_dfcount_sperm5_sig_fake_chr.txt"
SA_testes_firm_file_name <- "output_ann_MACAU_dffirm_testes9_sig_fake_chr.txt"

objs <- list(
  SA_blood_als_file_name, 
  SA_sperm_als_file_name, 
  SA_testes_als_file_name, 
  SA_sperm_count_file_name,
  SA_testes_firm_file_name)

output <- list()

for(j in 1:length(objs)){
  
  print(paste("Starting number", j))
  
  final_list <- vector()
  
  conn <- file(objs[[j]],open="r")
  linn <-readLines(conn)
  
  for (i in 1:length(linn)){
    line_split = strsplit(linn[i], "\t")
    final_list= c(final_list, line_split[[1]][length(line_split[[1]])])
  }
  
  close(conn)
  
  for (i in 1:length(final_list)) {
    if (final_list[i] == "promoter error: <2000 bp upstream of start site") {
      final_list[i] = "prom"}
  }
  
  result <- final_list %>%
    data.frame(x = .) %>%
    group_by(x) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    mutate(name = objs[[j]])
  
  output[[j]] <- result
  
  print(paste("Finished number", j))
}

# conn <- file(SA_blood_als_file_name,open="r")
# linn <-readLines(conn)
# for (i in 1:length(linn)){
#   line_split = strsplit(linn[i], "\t")
#   SA_blood_als= c(SA_blood_als, line_split[[1]][length(line_split[[1]])])
# }
# close(conn)
# for (i in 1:length(SA_blood_als)) {
#   if (SA_blood_als[i] == "promoter error: <2000 bp upstream of start site") {
#     SA_blood_als[i] = "prom"}
# }
# table(SA_blood_als)

p1 <- ggplot(output[[1]], aes(x = "", y = n, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = clrs) +
  labs(fill=NULL) +
  coord_polar("y", start=1.5) +
  theme_void() + theme(legend.position = "none")

p2 <- ggplot(output[[2]], aes(x = "", y = n, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = clrs) +
  labs(fill=NULL) +
  coord_polar("y", start=1.5) +
  theme_void() + theme(legend.position = "none")

p3 <- ggplot(output[[3]], aes(x = "", y = n, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = clrs) +
  labs(fill=NULL) +
  coord_polar("y", start=1.5) +
  theme_void() + theme(legend.position = "none")

p4 <- ggplot(output[[4]], aes(x = "", y = n, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = clrs) +
  labs(fill=NULL) +
  coord_polar("y", start=1.5) +
  theme_void() + theme(legend.position = "none")

p5 <- ggplot(output[[5]], aes(x = "", y = n, fill = x)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = clrs) +
  labs(fill=NULL) +
  coord_polar("y", start=1.5) +
  theme_void()

# combined plots
p_final <- (p1 | p2 | p3)/(p4 | p5) + plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 10))

p_final
ggsave(filename="site_annts_pies.pdf", path="/Users/stavrt/ferret/site_annt", device="pdf")

