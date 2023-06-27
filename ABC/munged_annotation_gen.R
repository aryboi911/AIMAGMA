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
celltype <- c("CD14-positive_monocyte_treated_with_LPS_4h-Novakovic2016", "CD14-positive_monocyte-ENCODE", "CD14-positive_monocyte-Novakovic2016", "CD14-positive_monocytes-Roadmap", "THP-1_macrophage-VanBortle2017", "THP-1_monocyte-VanBortle2017")

for(i in 1:22){
tsv <- read_tsv(paste0("/sc/arion/projects/LOAD/Papers/2018-08-22.AD.myeloid/sumstats.gwas/output/Jansen2019_phase3.chr",i,".CPRA_b37.tsv.gz"), col_names = FALSE)
tsv <- tsv[-(1:8), ]
tsv <- tsv %>% 
  separate(X1, into = c("ID", "CHROM", "POS", "REF", "ALT", "AF", "TRAIT", "BETA", "SE", "Z", "P", "N", "N_CASES", "N_CTRLS", "OR", "OR_L95", "OR_U95", "DIR", "P_HET", "G1000_ID", "G1000_VARIANT", "DBSNP_ID", "DBSNP_VARIANT", "OLD_ID", "OLD_VARIANT"), sep = "\t")
tsv$OLD_ID <- gsub("_", ":", tsv$OLD_ID)
unique(tsv) %>% write_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/jansen/Jansen2019_phase3.chr", i, ".CPRA_b37.tsv"))
}

##
for(z in 1:6){
  ct = celltype[z]

for(i in 1:22){
gene_chrom = paste0("chr", i)
annot_tsv = paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/annotations/jansen/", ct, "/abc_magma_", gene_chrom, ".genes.annot")
sumstats_bellenguez.hg19 = read_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/jansen/Jansen2019_phase3.", gene_chrom, ".CPRA_b37.tsv"))
sumstats_bellenguez.annot.hg19 = sumstats_bellenguez.hg19 %>% 
  dplyr::mutate(end = POS+1) %>%
  dplyr::select(snp=ID, chromosome=CHROM, start=POS, end=end)

sumstats_bellenguez.annot.hg19 <- sumstats_bellenguez.annot.hg19[complete.cases(sumstats_bellenguez.annot.hg19$start, sumstats_bellenguez.annot.hg19$end), ]

## make granges object
sumstats_bellenguez.annot.hg19 = sumstats_bellenguez.annot.hg19 %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = chromosome_lengths)

abc <- read_tsv(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/annotations/new_en/", ct, "/abc_ranges_updated_", gene_chrom, ".genes.annot"))

abc.annot.hg19 = abc %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = chromosome_lengths)

bellenguez_abc.annot = annotate_regions(sumstats_bellenguez.annot.hg19, abc.annot.hg19) %>% as_tibble

bellenguez_abc.annot %<>%
  dplyr::rename(range_id = annot.range_id, gene_loc = annot.gene_loc) %>%
  distinct(range_id, gene_loc, snp) %>%
  arrange(range_id, gene_loc, snp)

bellenguez_abc.annot %<>%
  group_by(range_id, gene_loc) %>%
  summarize(snps = str_c(snp, collapse = " "))

unique(bellenguez_abc.annot) %>% write.table(annot_tsv, row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
}
}

for(z in 1:6){
  ct = celltype[z]

gene <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/bellenguez/bellenguez_nott_chr1.genes.out", header = TRUE)

for(i in 2:22){
  gene1 <- read.table(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/bellenguez/bellenguez_nott_chr", i, ".genes.out"), header = TRUE)
  gene <- rbind(gene, gene1)
}

#gene <- separate(gene, GENE, into = c("GENE_LOC", "NR_ID", "GENE"), sep = "\\|")

unique(gene) %>% write.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/magma_output/bellenguez/bellenguez_nott_chrall.genes.out", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

}
