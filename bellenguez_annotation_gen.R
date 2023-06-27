suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(annotatr))
suppressMessages(library(readxl))
suppressMessages(library(biomaRt))
suppressMessages(library(GenomicRanges))
library(rtracklayer)
library(dplyr)
library(AnnotationHub)

## chromosomes for generating granges objects
chromosome_lengths.hg19 = c(249250621,
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
                            16571) #hg19

chromosome_lengths.hg38 = c(248956422,
                            242193529,
                            198295559,
                            190214555,
                            181538259,
                            170805979,
                            159345973,
                            145138636,
                            138394717,
                            133797422,
                            135086622,
                            133275309,
                            114364328,
                            107043718,
                            101991189,
                            90338345,
                            83257441,
                            80373285,
                            58617616,
                            64444167,
                            46709983,
                            50818468,
                            156040895,
                            57227415,
                            16569) #hg38

names(chromosome_lengths.hg19) = c(seq(1,22), c("X", "Y", "MT")) # seq(1,25)
names(chromosome_lengths.hg38) = c(seq(1,22), c("X", "Y", "MT")) # seq(1,25)

epig_xlsx = "/sc/arion/projects/load/users/patelt16/projects/aimagma/071422.aimagma/data/placseq_nott_2019/Nott2019_supplement_Table_S5.xlsx"
snatac_morabito.hg38 = "/sc/arion/projects/load/users/patelt16/projects/aimagma/071422.aimagma/data/snatac_morabito_2021/morabito_peaks_hg38_093021.txt"

for(i in 1:22){

sumstats_bellenguez.hg19 = read_tsv(paste0("/sc/arion/projects/load/users/pradha02/Projects/AIMAGMA/Bellenguez/hg19/recoded/Bellenguez_hg19_chrall", i, ".tsv"))


## format sumstats for granges input
sumstats_bellenguez.annot.hg19 = sumstats_bellenguez.hg19 %>% 
  dplyr::mutate(end = base_pair_location+1) %>%
  dplyr::select(snp=variant_alternate_id, chromosome, start=base_pair_location, end=end)

## make granges object
sumstats_bellenguez.annot.hg19 = sumstats_bellenguez.annot.hg19 %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = chromosome_lengths.hg19)


nott.promoter.annot <- read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/input/epig_nott2019/nott_2019.microglia.pu1.annot")

nott.promoter.annot.hg19 = nott.promoter.annot %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = chromosome_lengths.hg19)

bellenguez_nott.annot = annotate_regions(sumstats_bellenguez.annot.hg19, nott.promoter.annot.hg19) %>% as_tibble
bellenguez_nott.annot %<>%
  dplyr::rename(range_id = annot.range_id, gene_loc = annot.gene_loc) %>%
  distinct(range_id, gene_loc, snp) %>%
  arrange(range_id, gene_loc, snp)

bellenguez_nott.annot %<>%
  group_by(range_id, gene_loc) %>%
  summarize(snps = str_c(snp, collapse = " "))

unique(bellenguez_nott.annot) %>% write.table(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/nott/annotations/bellenguez_chr", i, ".genes.annot"), row.names = FALSE, col.names = FALSE, sep = " ", quote = FALSE)
}
