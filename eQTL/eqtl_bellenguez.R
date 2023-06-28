library(dplyr)
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(annotatr))
suppressMessages(library(readxl))
library(readr)
suppressMessages(library(biomaRt))
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))
suppressMessages(library(AnnotationHub))
suppressMessages(library(plotgardener))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg38.knownGene))
suppressMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressMessages(library(org.Hs.eg.db))

for(i in 1:22){
    eqtl.sumstats = read_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/eqtl_hg19/eqtl_data_chr", i, ".tsv"))
    sumstats_bellenguez = read_tsv(paste0("/sc/arion/projects/load/users/pradha02/Projects/AIMAGMA/Bellenguez/hg19/recoded/Bellenguez_hg19_chr", i, ".tsv"))
    colnames(eqtl.sumstats)[colnames(eqtl.sumstats) == "variant_id"] <- "variant_alternate_id"
    
    merged_df <- merge(sumstats_bellenguez, eqtl.sumstats, by = "variant_alternate_id") #, all.y = TRUE   
    merged_df$z_score <- abs(merged_df$beta.y / merged_df$se)
    median(merged_df$z_score)

    #merged_df$weighted_p <- merged_df$p_value / merged_df$z_score * median(merged_df$z_score)

    unique(merged_df) %>% write_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/eqtl/bellenguez_eqtl_chr", i, ".tsv"))

}
