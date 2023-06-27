suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
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

abc_txt = "/sc/arion/projects/load/users/patelt16/projects/aimagma/20230117.abc_scores/data/2021.Nasser_abc_scores/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt"
refseq_bed = "/sc/arion/projects/LOAD/Dado/datasets/Nasser2021/data/RefSeqCurated.170308.bed"

abc_ref <- read.table(abc_txt, header=TRUE, sep="\t")
refseq_ref <- rtracklayer::import.bed(refseq_bed) %>% as_tibble

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
##
for(z in 1:6){
  ct = celltype[z]

for(i in 11:22){
gene_chrom = paste0("chr", i)
annot_tsv = paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/annotations/new_en/", ct, "/abc_magma_", gene_chrom, ".genes.annot")
sumstats_bellenguez.hg19 = read_tsv(paste0("/sc/arion/projects/load/users/pradha02/Projects/AIMAGMA/Bellenguez/hg19/recoded/Bellenguez_hg19_", gene_chrom, ".tsv"))

abc <- abc_ref
refseq <- refseq_ref

abc %<>%
    filter(CellType == ct) %>%
    ## arrange(TargetGene, desc(ABC.Score)) %>%
    ## distinct(TargetGene, .keep_all = TRUE) %>%
    mutate(seqid = str_replace(chr, "chr", ""), strand = "*") %>%
    filter(seqid %in% str_c(1:22))

abc <- subset(abc, select = -chr)
abc <- subset(abc, select = c(seqid, start, end, strand, class, TargetGene, TargetGeneTSS, TargetGeneExpression,TargetGenePromoterActivityQuantile, TargetGeneIsExpressed, isSelfPromoter, distance, ABC.Score))

#############

abc %<>% GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = chromosome_lengths) %>% as_tibble

refseq$TargetGene <- refseq$name

refseq %<>%
    mutate(seqnames = str_replace(seqnames, "chr", "")) %>%
    filter(seqnames %in% str_c(1:22)) %>%
    unite("gene_loc", seqnames, start, end, sep = ":") %>%
    separate(TargetGene, c("TargetGene", "refseq"), sep = ";")
    
refseq <- subset(refseq, select = c(TargetGene, gene_loc, refseq))

refseq %<>%
    group_by(TargetGene) %>%
    summarize(refseqs = paste(refseq, collapse = ";"), gene_locs = paste(gene_loc, collapse = ";"), gene_loc = gene_loc[[1]]) %>%
    unite("range_id", gene_locs, refseqs, TargetGene, sep = "|", remove = FALSE)

refseq <- subset(refseq, select = c(TargetGene, range_id, gene_loc))
###

`%notin%` <- Negate(`%in%`)

stopifnot(sum(abc$TargetGene %notin% refseq$TargetGene) == 0)

abc %<>%
    left_join(refseq, by = "TargetGene")

abc$chrom <- abc$seqnames

abc <- subset(abc, select = c(range_id, chrom, start, end, gene_loc))

enhancers <- read.table(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/regulatory_elements/", gene_chrom, ".txt"), header = TRUE)
enhancers <- separate(enhancers, reserved, into = c("r", "g", "b"), sep = ",")
enhancers <- enhancers %>%
  filter(as.numeric(r) == as.numeric(g))


abc <- subset(abc, chrom == i)

for ( j in 1:nrow(abc)){
  for( k in 1:nrow(enhancers)){
    if((abc[j, "start"] > enhancers$chromStart[k] & abc[j, "end"] < enhancers$chromEnd[k])){
      if(abc[j, "start"] > enhancers$chromStart[k]){ abc[j, "start"] <- enhancers$chromStart[k]}
      if(abc[j, "end"] < enhancers$chromEnd[k]){ abc[j, "end"] <- enhancers$chromEnd[k]}
    }
  }
}

unique(abc) %>% write.table(paste0("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/annotations/new_en/", ct, "/abc_ranges_updated_", gene_chrom, ".genes.annot"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


sumstats_bellenguez.annot.hg19 = sumstats_bellenguez.hg19 %>% 
  dplyr::mutate(end = base_pair_location+1) %>%
  dplyr::select(snp=variant_alternate_id, chromosome, start=base_pair_location, end=end)

## make granges object
sumstats_bellenguez.annot.hg19 = sumstats_bellenguez.annot.hg19 %>% 
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = chromosome_lengths)


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
