library(rtracklayer)
library(GenomicRanges)
library(rtracklayer)
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
library(qqman)
suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

chain_file_path <- "/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/eqtl_hg19/hg38ToHg19_tabs.over.chain"
chain <- import.chain(chain_file_path)

for(i in 1:22){
    input_data <- read_tsv(paste0("/sc/arion/projects/load/users/edwart10/projects/xQTL-STARNET-Analysis-05-24-2022/QTL-Association-Testing/cisQTL-Analysis-Workflows/TensorQTL-QTL-Association-Testing/output/10-12-2022/Association_Scan/Known_Cov_Status_BiCV/STARNET.mol_phe.resid.bed.processed_phenotype.per_chrom_STARNET.tpm.gct.BiCV.cov.", i, ".norminal.cis_long_table.txt"))
    input_data$variant_id <- gsub(pattern = "_", replacement = ":", x = input_data$variant_id)
    input_data$variant_id <- gsub(pattern = "chr", replacement = "", x = input_data$variant_id)

gr <- GRanges(seqnames = input_data$chrom,
              ranges = IRanges(start = input_data$pos, end = input_data$pos),
              strand = "*")
    gr_lifted <- liftOver(gr, chain)
    lifted_positions <- start(gr_lifted)

    input_data$pos <- as.integer(lifted_positions)
    input_data$variant_id <- gsub(pattern = "_", replacement = ":", x = input_data$variant_id)
    input_data$variant_id <- substr(input_data$variant_id, start = nchar(input_data$variant_id) - 3, stop = nchar(input_data$variant_id))
    input_data$variant_id <- paste0(i,":", input_data$pos, input_data$variant_id)
    input_data$z_score <- abs(input_data$beta / input_data$se)
    input_data$n <- input_data$n * 100
    input_data <- na.omit(input_data, cols = "pos")
    unique(input_data) %>% write_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/eqtl_hg38/eqtl_chr", i, ".tsv"))
}
