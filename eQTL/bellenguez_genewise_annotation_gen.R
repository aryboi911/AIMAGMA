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

chromosome_lengths = c(249250621,
                       243199373,
                       198022430,
                       191154276,
                       180915260,
                       171115067,
                       159138663,
                       146364022,
                       141213431,
                       135534747,
                       135006516,
                       133851895,
                       115169878,
                       107349540,
                       102531392,
                       90354753,
                       81195210,
                       78077248,
                       59128983,
                       63025520,
                       48129895,
                       51304566,
                       155270560,
                       59373566,
                       16571)

names(chromosome_lengths) = c(seq(1,22), c("X", "Y", "MT")) # seq(1,25)

for(i in 1:22){
  gene_chrom = paste0("chr", i)
  eqtl.sumstats = read_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/eqtl/bellenguez_eqtl_chr", i, ".tsv"))

  list_of_df <- split(eqtl.sumstats, eqtl.sumstats$molecular_trait_id)

  mart = useEnsembl("ensembl", "hsapiens_gene_ensembl", version = "GRCh37")
  ## Query biomart to get all genes and their ensembl ID, names, coordinates
  t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position'), mart = mart, useCache = FALSE)

  unique_ids <- unique(eqtl.sumstats$molecular_trait_id)
  unique_ids_df <- data.frame(molecular_trait_id = unique_ids)
  colnames(t2g)[colnames(t2g) == "ensembl_gene_id"] <- "molecular_trait_id"
  merged_df <- merge(unique_ids_df, t2g, by = "molecular_trait_id")

  for(z in 1:length(list_of_df)){
    molecular_trait_id <- names(list_of_df)[z]
    merged_row <- merged_df[merged_df$molecular_trait_id == molecular_trait_id, ]
    
    file_string <- paste(
      molecular_trait_id,
      paste0(merged_row$chromosome_name, ":", merged_row$start_position, ":", merged_row$end_position),
      paste(list_of_df[[z]]$variant_alternate_id, collapse = " "),
      sep = " "
    )
    sumstats <- list_of_df[[molecular_trait_id]]

    unique(sumstats) %>% write_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/annotations/eqtl_gene/sumstats/chr", i, "/", molecular_trait_id, ".tsv"))
  
    write(file_string, file = paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/annotations/eqtl_gene/chr", i, "/", molecular_trait_id, ".genes.annot"))
  }
}
