library(dplyr)
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(annotatr))
suppressMessages(library(readxl))
suppressMessages(library(biomaRt))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(AnnotationHub))
suppressMessages(library(plotgardener))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(org.Hs.eg.db))

df1 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/new_en/CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016/nasser_bellenguez_chrall.genes.out", header = TRUE)
df2 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/new_en/CD14-positive_monocyte-ENCODE/nasser_bellenguez_chrall.genes.out", header = TRUE)
df3 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/new_en/CD14-positive_monocyte-Novakovic2016/nasser_bellenguez_chrall.genes.out", header = TRUE)
df4 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/new_en/CD14-positive_monocytes-Roadmap/nasser_bellenguez_chrall.genes.out", header = TRUE)
df5 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/new_en/THP-1_macrophage-VanBortle2017/nasser_bellenguez_chrall.genes.out", header = TRUE)
df6 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/new_en/THP-1_monocyte-VanBortle2017/nasser_bellenguez_chrall.genes.out", header = TRUE)
df7 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/excluding23andMe/CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016/output_chrall.genes.out", header = TRUE)
df8 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/excluding23andMe/CD14-positive_monocyte-ENCODE/output_chrall.genes.out", header = TRUE)
df9 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/excluding23andMe/CD14-positive_monocyte-Novakovic2016/output_chrall.genes.out", header = TRUE)
df10 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/excluding23andMe/CD14-positive_monocytes-Roadmap/output_chrall.genes.out", header = TRUE)
df11 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/excluding23andMe/THP-1_macrophage-VanBortle2017/output_chrall.genes.out", header = TRUE)
df12 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/excluding23andMe/THP-1_monocyte-VanBortle2017/output_chrall.genes.out", header = TRUE)
df13 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/including23andMe/CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016/output_chrall.genes.out", header = TRUE)
df14 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/including23andMe/CD14-positive_monocyte-ENCODE/output_chrall.genes.out", header = TRUE)
df15 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/including23andMe/CD14-positive_monocyte-Novakovic2016/output_chrall.genes.out", header = TRUE)
df16 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/including23andMe/CD14-positive_monocytes-Roadmap/output_chrall.genes.out", header = TRUE)
df17 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/including23andMe/THP-1_macrophage-VanBortle2017/output_chrall.genes.out", header = TRUE)
df18 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/including23andMe/THP-1_monocyte-VanBortle2017/output_chrall.genes.out", header = TRUE)
df19 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/noUKB/CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016.genes.out", header = TRUE)
df19 <- separate(df19, GENE, into = c("GENE_LOC", "NR_ID", "GENE"), sep = "\\|")
df20 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/noUKB/CD14-positive_monocyte-ENCODE.genes.out", header = TRUE)
df20 <- separate(df20, GENE, into = c("GENE_LOC", "NR_ID", "GENE"), sep = "\\|")
df21 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/noUKB/CD14-positive_monocyte-Novakovic2016.genes.out", header = TRUE)
df21 <- separate(df21, GENE, into = c("GENE_LOC", "NR_ID", "GENE"), sep = "\\|")
df22 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/noUKB/CD14-positive_monocytes-Roadmap.genes.out", header = TRUE)
df22 <- separate(df22, GENE, into = c("GENE_LOC", "NR_ID", "GENE"), sep = "\\|")
df23 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/noUKB/THP-1_macrophage-VanBortle2017.genes.out", header = TRUE)
df23 <- separate(df23, GENE, into = c("GENE_LOC", "NR_ID", "GENE"), sep = "\\|")
df24 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/wightman/noUKB/THP-1_monocyte-VanBortle2017.genes.out", header = TRUE)
df24 <- separate(df24, GENE, into = c("GENE_LOC", "NR_ID", "GENE"), sep = "\\|")
df25 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/jansen/CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016.genes.out", header = TRUE)
df26 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/jansen/CD14-positive_monocyte-ENCODE.genes.out", header = TRUE)
df27 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/jansen/CD14-positive_monocyte-Novakovic2016.genes.out", header = TRUE)
df28 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/jansen/CD14-positive_monocytes-Roadmap.genes.out", header = TRUE)
df29 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/jansen/THP-1_macrophage-VanBortle2017.genes.out", header = TRUE)
df30 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/jansen/THP-1_monocyte-VanBortle2017.genes.out", header = TRUE)
df31 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/bellenguez/bellenguez_nott_chrall.genes.out", header = TRUE)
df31 <- separate(df31, GENE, into = c("ENSG_ID", "GENE"), sep = ":")
df32 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/jansen/jansen.genes.out", header = TRUE)
df32 <- separate(df32, GENE, into = c("ENSG_ID", "GENE"), sep = ":")
df33 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/wightman/PGCALZ2sumstatsExcluding23andMe_chrall.genes.out", header = TRUE)
df33 <- separate(df33, GENE, into = c("ENSG_ID", "GENE"), sep = ":")
df34 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/wightman/PGCALZ2sumstatsIncluding23andMeMETALNoUKB_chrall.genes.out", header = TRUE)
df34 <- separate(df34, GENE, into = c("ENSG_ID", "GENE"), sep = ":")
df35 <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/wightman/PGCALZ2sumstatsIncluding23andMe_chrall.genes.out", header = TRUE)
df35 <- separate(df35, GENE, into = c("ENSG_ID", "GENE"), sep = ":")



df1 <- df1[c("GENE", "CHR", "START", "STOP", "P_MULTI")] 
names(df1)[names(df1) == "P_MULTI"] <- "P_BELLENGUEZ_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016"
df2 <- df2[c("GENE", "P_MULTI")] 
names(df2)[names(df2) == "P_MULTI"] <- "P_BELLENGUEZ_CD14_positive_monocyte_ENCODE"
df3 <- df3[c("GENE", "P_MULTI")] 
names(df3)[names(df3) == "P_MULTI"] <- "P_BELLENGUEZ_CD14_positive_monocyte_Novakovic2016"
df4 <- df4[c("GENE", "P_MULTI")] 
names(df4)[names(df4) == "P_MULTI"] <- "P_BELLENGUEZ_CD14_positive_monocytes_Roadmap"
df5 <- df5[c("GENE", "P_MULTI")] 
names(df5)[names(df5) == "P_MULTI"] <- "P_BELLENGUEZ_THP_1_macrophage_VanBortle2017"
df6 <- df6[c("GENE", "P_MULTI")] 
names(df6)[names(df6) == "P_MULTI"] <- "P_BELLENGUEZ_THP_1_monocyte_VanBortle2017"
df7 <- df7[c("GENE", "P_MULTI")] 
names(df7)[names(df7) == "P_MULTI"] <- "P_WIGHTMAN_EXCLUDE23ANDME_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016"
df8 <- df8[c("GENE", "P_MULTI")] 
names(df8)[names(df8) == "P_MULTI"] <- "P_WIGHTMAN_EXCLUDE23ANDME_CD14_positive_monocyte_ENCODE"
df9 <- df9[c("GENE", "P_MULTI")] 
names(df9)[names(df9) == "P_MULTI"] <- "P_WIGHTMAN_EXCLUDE23ANDME_CD14_positive_monocyte_Novakovic2016"
df10 <- df10[c("GENE", "P_MULTI")] 
names(df10)[names(df10) == "P_MULTI"] <- "P_WIGHTMAN_EXCLUDE23ANDME_CD14_positive_monocytes_Roadmap"
df11 <- df11[c("GENE", "P_MULTI")] 
names(df11)[names(df11) == "P_MULTI"] <- "P_WIGHTMAN_EXCLUDE23ANDME_THP_1_macrophage_VanBortle2017"
df12 <- df12[c("GENE", "P_MULTI")] 
names(df12)[names(df12) == "P_MULTI"] <- "P_WIGHTMAN_EXCLUDE23ANDME_THP_1_monocyte_VanBortle2017"
df13 <- df13[c("GENE", "P_MULTI")] 
names(df13)[names(df13) == "P_MULTI"] <- "P_WIGHTMAN_INCLUDE23ANDME_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016"
df14 <- df14[c("GENE", "P_MULTI")] 
names(df14)[names(df14) == "P_MULTI"] <- "P_WIGHTMAN_INCLUDE23ANDME_CD14_positive_monocyte_ENCODE"
df15 <- df15[c("GENE", "P_MULTI")] 
names(df15)[names(df15) == "P_MULTI"] <- "P_WIGHTMAN_INCLUDE23ANDME_CD14_positive_monocyte_Novakovic2016"
df16 <- df16[c("GENE", "P_MULTI")] 
names(df16)[names(df16) == "P_MULTI"] <- "P_WIGHTMAN_INCLUDE23ANDME_CD14_positive_monocytes_Roadmap"
df17 <- df17[c("GENE", "P_MULTI")] 
names(df17)[names(df17) == "P_MULTI"] <- "P_WIGHTMAN_INCLUDE23ANDME_THP_1_macrophage_VanBortle2017"
df18 <- df18[c("GENE", "P_MULTI")] 
names(df18)[names(df18) == "P_MULTI"] <- "P_WIGHTMAN_INCLUDE23ANDME_THP_1_monocyte_VanBortle2017"
df19 <- df19[c("GENE", "P_MULTI")] 
names(df19)[names(df19) == "P_MULTI"] <- "P_WIGHTMAN_NOUKB_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016"
df20 <- df20[c("GENE", "P_MULTI")] 
names(df20)[names(df20) == "P_MULTI"] <- "P_WIGHTMAN_NOUKB_CD14_positive_monocyte_ENCODE"
df21 <- df21[c("GENE", "P_MULTI")] 
names(df21)[names(df21) == "P_MULTI"] <- "P_WIGHTMAN_NOUKB_CD14_positive_monocyte_Novakovic2016"
df22 <- df22[c("GENE", "P_MULTI")] 
names(df22)[names(df22) == "P_MULTI"] <- "P_WIGHTMAN_NOUKB_CD14_positive_monocytes_Roadmap"
df23 <- df23[c("GENE", "P_MULTI")] 
names(df23)[names(df23) == "P_MULTI"] <- "P_WIGHTMAN_NOUKB_THP_1_macrophage_VanBortle2017"
df24 <- df24[c("GENE", "P_MULTI")] 
names(df24)[names(df24) == "P_MULTI"] <- "P_WIGHTMAN_NOUKB_THP_1_monocyte_VanBortle2017"
df25 <- df25[c("GENE", "P_MULTI")] 
names(df25)[names(df25) == "P_MULTI"] <- "P_JANSEN_CD14_positive_monocyte_treated_with_LPS_4h_Novakovic2016"
df26 <- df26[c("GENE", "P_MULTI")] 
names(df26)[names(df26) == "P_MULTI"] <- "P_JANSEN_CD14_positive_monocyte_ENCODE"
df27 <- df27[c("GENE", "P_MULTI")] 
names(df27)[names(df27) == "P_MULTI"] <- "P_JANSEN_CD14_positive_monocyte_Novakovic2016"
df28 <- df28[c("GENE", "P_MULTI")] 
names(df28)[names(df28) == "P_MULTI"] <- "P_JANSEN_CD14_positive_monocytes_Roadmap"
df29 <- df29[c("GENE", "P_MULTI")] 
names(df29)[names(df29) == "P_MULTI"] <- "P_JANSEN_THP_1_macrophage_VanBortle2017"
df30 <- df30[c("GENE", "P_MULTI")] 
names(df30)[names(df30) == "P_MULTI"] <- "P_JANSEN_THP_1_monocyte_VanBortle2017"
df31 <- df31[c("GENE", "P_MULTI")] 
names(df31)[names(df31) == "P_MULTI"] <- "P_BELLENGUEZ_Nott"
df32 <- df32[c("GENE", "P_MULTI")] 
names(df32)[names(df32) == "P_MULTI"] <- "P_JANSEN_Nott"
df33 <- df33[c("GENE", "P_MULTI")] 
names(df33)[names(df33) == "P_MULTI"] <- "P_WIGHTMAN_NOUKB_Nott"
df34 <- df34[c("GENE", "P_MULTI")]
names(df34)[names(df34) == "P_MULTI"] <- "P_WIGHTMAN_EXCLUDING23ANDME_Nott"
df35 <- df35[c("GENE", "P_MULTI")]
names(df35)[names(df35) == "P_MULTI"] <- "P_WIGHTMAN_INCLUDING23ANDME_Nott"



large_df <- merge(df1, df2, by = "GENE", all = TRUE)
large_df <- merge(large_df, df3, by = "GENE", all = TRUE)
large_df <- merge(large_df, df4, by = "GENE", all = TRUE)
large_df <- merge(large_df, df5, by = "GENE", all = TRUE)
large_df <- merge(large_df, df6, by = "GENE", all = TRUE)
large_df <- merge(large_df, df7, by = "GENE", all = TRUE)
large_df <- merge(large_df, df8, by = "GENE", all = TRUE)
large_df <- merge(large_df, df9, by = "GENE", all = TRUE)
large_df <- merge(large_df, df10, by = "GENE", all = TRUE)
large_df <- merge(large_df, df11, by = "GENE", all = TRUE)
large_df <- merge(large_df, df12, by = "GENE", all = TRUE)
large_df <- merge(large_df, df13, by = "GENE", all = TRUE)
large_df <- merge(large_df, df14, by = "GENE", all = TRUE)
large_df <- merge(large_df, df15, by = "GENE", all = TRUE)
large_df <- merge(large_df, df16, by = "GENE", all = TRUE)
large_df <- merge(large_df, df17, by = "GENE", all = TRUE)
large_df <- merge(large_df, df18, by = "GENE", all = TRUE)
large_df <- merge(large_df, df19, by = "GENE", all = TRUE)
large_df <- merge(large_df, df20, by = "GENE", all = TRUE)
large_df <- merge(large_df, df21, by = "GENE", all = TRUE)
large_df <- merge(large_df, df22, by = "GENE", all = TRUE)
large_df <- merge(large_df, df23, by = "GENE", all = TRUE)
large_df <- merge(large_df, df24, by = "GENE", all = TRUE)
large_df <- merge(large_df, df25, by = "GENE", all = TRUE)
large_df <- merge(large_df, df26, by = "GENE", all = TRUE)
large_df <- merge(large_df, df27, by = "GENE", all = TRUE)
large_df <- merge(large_df, df28, by = "GENE", all = TRUE)
large_df <- merge(large_df, df29, by = "GENE", all = TRUE)
large_df <- merge(large_df, df30, by = "GENE", all = TRUE)
large_df <- merge(large_df, df31, by = "GENE", all = TRUE)
large_df <- merge(large_df, df32, by = "GENE", all = TRUE)
large_df <- merge(large_df, df33, by = "GENE", all = TRUE)
large_df <- merge(large_df, df34, by = "GENE", all = TRUE)
large_df <- merge(large_df, df35, by = "GENE", all = TRUE)


mart = useEnsembl("ensembl", "hsapiens_gene_ensembl", version = "GRCh37")

## Query biomart to get all genes and their ensembl ID, names, coordinates
t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position', 'external_gene_name'), mart = mart, useCache = FALSE)

colnames(large_df)[colnames(large_df) == "GENE"] <- "external_gene_name"
large_df_merge <- merge(large_df, t2g, by = "external_gene_name", all.x = TRUE)

large_df_merge <- large_df_merge[grepl('^[0-9]+$', large_df_merge$chromosome_name), ]

fisher_combine <- function(p_values) {
  fisher_statistic <- -2 * sum(log(na.omit(p_values)))
  fisher_p_value <- pchisq(fisher_statistic, df = 2*length(na.omit(p_values)), lower.tail = FALSE)
  return(fisher_p_value*1e50)
}

calc_num_signif <- function(p_values) {
  prop_signif <- sum(p_values < 5e-08, na.rm = TRUE)
  return(prop_signif)
}

calc_prop_signif <- function(p_values) {
  prop_signif <- sum(p_values < 5e-08, na.rm = TRUE) / sum(!is.na(p_values))
  return(prop_signif)
}


fisher_combined_p_value <- function(p_values) {
    num_significant <- sum(p_values < 5e-08, na.rm = TRUE) 
    combined_p_value <- metap::metap(p_values, method="fisher")$p.value
      return(combined_p_value)
}

large_df_merge$num_significant <- apply(large_df_merge[,5:39], 1, calc_num_signif)
large_df_merge$prop_significant <- apply(large_df_merge[,5:39], 1, calc_prop_signif)
large_df_merge$overall_p_value <- apply(large_df_merge[,5:39], 1, fisher_combine)
large_df_merge$chromosome_name <- as.numeric(large_df_merge$chromosome_name)

unique(large_df_merge) %>% write_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/magma_output/large_table_output_chrall.genes.out"))
