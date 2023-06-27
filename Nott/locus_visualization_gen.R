##2023-01-30

############# EDIT THIS #############
gene_list <- read.table("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/gene_list.txt", header = FALSE, sep = " ")
output <- "/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/viz_out/multi3.pdf"
#####################################

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

pdf(output)

mart = useEnsembl("ensembl", "hsapiens_gene_ensembl", version = "GRCh37")

## Query biomart to get all genes and their ensembl ID, names, coordinates
t2g<-getBM(attributes=c('ensembl_gene_id',"ensembl_gene_id_version",'chromosome_name','start_position','end_position', 'external_gene_name'), mart = mart, useCache = FALSE)

colnames(gene_list)[colnames(gene_list) == "V1"] <- "ensembl_gene_id"
gene_list <- merge(gene_list, t2g, by = "ensembl_gene_id")
gene_list$chrom <- paste0("chr", gene_list$chromosome_name)


for (z in 1:nrow(gene_list)){
  gene_chrom <- gene_list[z, ]$chrom
  gene_start <- gene_list[z, ]$start_position
  gene_end <- gene_list[z, ]$end_position
  gene_name <- gene_list[z, ]$external_gene_name



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

epig_xlsx = "/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/placseq_nott_2019/Nott2019_supplement_Table_S5.xlsx"
snatac_morabito.hg38 = "/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/snatac_morabito_2021/morabito_peaks_hg38_093021.txt"

promoters = read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/nott_microglial_promoters.tsv", col_names = F)
promoters <- subset(promoters, ((X2 > gene_start & X2 < gene_end) | (X3 > gene_start & X3 < gene_end) |
                                  (X2 < gene_start & X3 > gene_end)) & X1 == gene_chrom)

if (nrow(promoters >0)){
  for(i in 1:nrow(promoters)){
dnaLoops <- read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/placseq_nott_2019/nott2019_plaqseq_microglia_chromatin_interactions_tables5.tsv")
dnaLoops <- subset(dnaLoops, ((start1 > promoters[i, ]$X2 & start1 < promoters[i, ]$X3) | (end1 > promoters[i, ]$X2 & end1 < promoters[i, ]$X3) |
                                (start2 > promoters[i, ]$X2 & start2 < promoters[i, ]$X3) | (end2 > promoters[i, ]$X2 & end2 < promoters[i, ]$X3) |
                                (start1 < promoters[i, ]$X2 & end1 > promoters[i, ]$X3) | (start2 < promoters[i, ]$X2 & end2 > promoters[i, ]$X3)) & chr1 == gene_chrom)
  }
enhancers = read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/nott_microglial_enhancers.tsv", col_names = F)
enhancers <- subset(enhancers, enhancers$X1 == gene_chrom)

filtered_en <- data.frame()

if(nrow(dnaLoops > 0)){
for (i in 1:nrow(enhancers)) {
  dna <- dnaLoops
  dna <- subset(dna, (start1 > enhancers[i,]$X2  & start1 < enhancers[i,]$X3 ) | (end1 > enhancers[i,]$X2  & end1 < enhancers[i,]$X3 ) |
                  (start2 > enhancers[i,]$X2  & start2 < enhancers[i,]$X3 ) | (end2 > enhancers[i,]$X2  & end2 < enhancers[i,]$X3 ) |
                  (start1 < enhancers[i,]$X2  & end1 > enhancers[i,]$X3 ) | (start2 < enhancers[i,]$X2  & end2 > enhancers[i,]$X3 ))
  if (nrow(dna) > 0) {
    en <- enhancers
    en <- subset(en, ((dna$start1  > X2 & dna$start1  < X3) | (dna$end1  > X2 & dna$end1  < X3) |
                        (dna$start2  > X2 & dna$start2  < X3) | (dna$end2  > X2 & dna$end2  < X3) |
                        (dna$start1  < X2 & dna$end1  > X3) | (dna$start2  < X2 & dna$end2  > X3)))
    filtered_en <- rbind(filtered_en, en)
  }
}
}
filtered_en <- unique(filtered_en)
promoters <- unique(promoters)


sumstats_bellenguez.hg19 = read_tsv(paste0("/sc/arion/projects/load/users/pradha02/Projects/AIMAGMA/Bellenguez/hg19/recoded/Bellenguez_hg19_", gene_chrom, ".tsv"))

chromo <- gene_chrom

enhancers <- filtered_en

min <- c(floor(min(promoters$X2)/100000), floor(min(enhancers$X2)/100000), floor(min(dnaLoops$start1)/100000))
start_reg <- min(min)*100000
max <- c(ceiling(max(promoters$X3)/100000), ceiling(max(enhancers$X3)/100000), ceiling(max(dnaLoops$end2)/100000))
end_reg <- max(max)*100000

sumstats.gwas = read_tsv(paste0("/sc/arion/projects/load/users/pradha02/Projects/AIMAGMA/Bellenguez/hg19/recoded/Bellenguez_hg19_", chromo, ".tsv"))
sumstats.gwas = sumstats.gwas %>% separate(variant_alternate_id, c("chrom","pos","a1","a2")) 
sumstats.gwas = sumstats.gwas %>% mutate(snp = paste0("chr",chrom,":",pos))

## ld scores from 1000 genomes LDassoc
sumstats_ld = read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/input/ld/chr11.picalm.ld.tsv")
sumstats_ld = sumstats_ld %>% 
  separate(Coord, c("chrom","pos")) %>% 
  mutate(snp = paste(chrom, pos, sep=":"))

## merge sumstats with ld scores
gwas_ld <- sumstats.gwas
gwas_ld = gwas_ld %>% mutate(chrom=paste0("chr",chrom)) %>%
  dplyr::select(snp=variant_id, chrom, pos, p=p_value)
gwas_ld$pos = as.integer(gwas_ld$pos)

max_snp <- subset(gwas_ld, pos > start_reg  & pos < end_reg & p != 0)
maxVal = (ceiling((((log10(min(max_snp$p))*-1)*1.15)/10)))*10
Yaxis = (ceiling(((log10(min(max_snp$p))*-1)/10)))*10
mostSigPos <- as.numeric(max_snp$pos[which.min(max_snp$p)])


min <- c(floor(min(promoters$X2)/100000), floor(min(enhancers$X2)/100000), floor(min(dnaLoops$start1)/100000), floor(gene_start/100000))
start_reg <- min(min)*100000
max <- c(ceiling(max(promoters$X3)/100000), ceiling(max(enhancers$X3)/100000), ceiling(max(dnaLoops$end2)/100000), ceiling(gene_end/100000))
end_reg <- max(max)*100000

wid = (end_reg - start_reg) / 333.333333
diff = wid * 1.66666

## create a plotgardener page
pageCreate(width = 7, height = 5.75, default.units = "inches", showGuides = FALSE)


region <- pgParams(chrom = chromo, 
                   assembly = "hg19",
                   chromstart = start_reg, 
                   chromend = end_reg)

## add genes track within selected region and highlight gene of interest
genes_a <- plotGenes(params = region, 
                     stroke = 1, 
                     fontsize = 6, 
                     assembly = "hg19",
                     y = 0, 
                     height = 1.0, x = 0.75, width = 6, 
                     just = c("left", "top"),
                     geneHighlights = data.frame("gene" = c(gene_name), "color" = c("#225EA8")),
                     geneBackground = "grey", default.units = "inches") 

# pagePlotRemove(genes_a)

## add bp chromosome position label at bottom
annoGenomeLabel(plot = genes_a, 
                x = 0.75, y = 0.95, #4.6, 
                fontsize = 8, 
                just = c("left", "top"), 
                default.units = "inches")

h3k27ac_peaks <- read_tsv(file = paste0("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/human_PU1nuclei_H3K27ac_epilepsy_pooled_hg19.", chromo, ".ucsc.bed"), col_names = FALSE) %>% unique(); colnames(h3k27ac_peaks) = c("chrom", "start", "end", "peak_id", "score")

nott_signal_nonoverlap <- h3k27ac_peaks %>%
  filter(start >= start_reg & end <= end_reg) %>%
  unique() %>%
  mutate(end_original = end, end = end-0.1)

signal_bar_plot <- data.frame()

for (i in seq(from = start_reg, to = end_reg, by = diff)) {
  small_df <- subset(nott_signal_nonoverlap, start > i & end < (i+diff))
  max_index <- which.max(small_df$score)
  small_df <- small_df[max_index,]
  signal_bar_plot <- rbind(signal_bar_plot, small_df)
}

signal_bar_plot$score <- ifelse(signal_bar_plot$score < 4, 0, signal_bar_plot$score)

signal <- ggplot(data = signal_bar_plot, aes(x = start, y = score)) +
  geom_bar(width = wid, stat = "identity", fill = "blue") +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

plotGG(
  plot = signal, x = 0.43, y = 1.35, width = 6.64, height = 0.65,
  just = c("left", "top"), default.units = "inches"
)

plotSegments(
  x0 = 0.75, y0 = 1.93, x1 = 6.75, y1 = 1.93,
  default.units = "inches",
  lwd = 1, lty = 1, linecolor = "blue"
)

# pagePlotRemove(signal1)



total_enhancers = read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/nott_microglial_enhancers.tsv", col_names = F)
total_enhancers <- subset(total_enhancers, total_enhancers$X1 == gene_chrom)

plotRanges(total_enhancers,
           chrom = chromo,
           assembly = "hg19",
           fill = "grey",
           boxHeight = unit(5, "mm"),
           x = 0.75, y = 2.15,
           width = 6, height = 0.5,
           just = c("left", "top"),
           params = region)

plotText(label = "Microglial H3K27Ac ChIP-seq peaks - Nott et al (2019)", 
         fontsize = 8, 
         fontcolor = c("black", "#37a7db"),
         x = 0.75, y = 1.25, 
         just = "left, top", 
         default.units = "inches")

plotText(label = "Microglia promoter/enhancer elements - Nott et al (2019)", 
         fontsize = 8, 
         fontcolor = c("black", "#37a7db"),
         x = 0.75, y = 2.05, 
         just = "left, top", 
         default.units = "inches")

if(nrow(promoters) > 0){

plotRanges(promoters,
           chrom = chromo,
           assembly = "hg19",
           fill = "blue",
           boxHeight = unit(5, "mm"),
           x = 0.75, y = 2.15,
           width = 6, height = 0.5,
           just = c("left", "top"),
           params = region)

}

if(nrow(enhancers > 0 )){
plotRanges(enhancers,
           chrom = chromo,
           assembly = "hg19",
           fill = "red",
           boxHeight = unit(5, "mm"),
           x = 0.75, y = 2.15,
           width = 6, height = 0.5,
           just = c("left", "top"),
           params = region)

}

if(nrow(dnaLoops) > 0){

dnaLoops$length <- (dnaLoops$start2 - dnaLoops$start1) / 1000
# dnaLoops <- subset(dnaLoops, ((start1 >= 85600000 & end2 <= 86100000) & chrom == "chr11"))
# dnaLoops <- subset(dnaLoops, count > 50)
dnaLoops$h <- dnaLoops$length / max(dnaLoops$length)

## Add a length column
dnaLoops$length <- (dnaLoops$start2 - dnaLoops$start1) / 1000

## plot pairs of arches
archPlot <- plotPairsArches(data = dnaLoops, 
                            params = region,
                            chrom = chromo,
                            chromstart = start_reg,
                            chromend = end_reg,
                            fill = colorby("fdr", palette = 
                                             colorRampPalette(c("dodgerblue2", "firebrick2"))),
                            linecolor = "fill",
                            archHeight = "h", alpha = 1, 
                            assembly = "hg19", 
                            width = 6, height = 0.4,
                            x = 0.75, y = 2.9,
                            just = c("left", "top"),
                            default.units = "inches", flip = TRUE)
if(nrow(dnaLoops) > 1){
annoHeatmapLegend(
  plot = archPlot, fontcolor = "black",
   x = 0.35, y = 2.85,
   width = 0.10, height = 0.5, fontsize = 0
)

 plotText(
   label = "FDR", rot = 90, x = 0.55, y = 3.1,
   just = c("center", "center"),
   fontsize = 8
 )
 plotText(
   label = format(min(dnaLoops$fdr), scientific = TRUE, digits = 3), x = 0.4, y = 3.39,
   just = c("center", "center"),
   fontsize = 5
 )

 plotText(
   label = format(max(dnaLoops$fdr), scientific = TRUE, digits = 3), x = 0.4, y = 2.81,
   just = c("center", "center"),
   fontsize = 5
 )
}
}



## text tracks

plotText(label = "Microglial chromatin interactions - Nott et al (2019)", 
         fontsize = 8, 
         fontcolor = c("black", "#37a7db"),
         x = 0.75, y = 2.75, 
         just = "left, top", 
         default.units = "inches")

plotText(label = "AD GWAS summary statistics LD - Bellenguez et al (2021)", 
         fontsize = 8, 
         fontcolor = c("black", "#37a7db"),
         x = 0.75, y = 3.4, 
         just = "left, top", 
         default.units = "inches")

## gwas sumstats
sumstats.gwas = read_tsv(paste0("/sc/arion/projects/load/users/pradha02/Projects/AIMAGMA/Bellenguez/hg19/recoded/Bellenguez_hg19_", chromo, ".tsv"))
sumstats.gwas = sumstats.gwas %>% separate(variant_alternate_id, c("chrom","pos","a1","a2")) 
sumstats.gwas = sumstats.gwas %>% mutate(snp = paste0("chr",chrom,":",pos))

## ld scores from 1000 genomes LDassoc
sumstats_ld = read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/input/ld/chr11.picalm.ld.tsv")
sumstats_ld = sumstats_ld %>% 
  separate(Coord, c("chrom","pos")) %>% 
  mutate(snp = paste(chrom, pos, sep=":"))

## merge sumstats with ld scores
gwas_ld <- sumstats.gwas
gwas_ld = gwas_ld %>% mutate(chrom=paste0("chr",chrom)) %>%
  dplyr::select(snp=variant_id, chrom, pos, p=p_value)
gwas_ld$pos = as.integer(gwas_ld$pos)

max_snp <- subset(gwas_ld, pos > start_reg  & pos < end_reg & p != 0)
maxVal = (ceiling((((log10(min(max_snp$p))*-1)*1.1)/10)))*10
Yaxis = (ceiling(((log10(min(max_snp$p))*-1)/10)))*10
mostSigRSID <- gsub("\"", "", paste(max_snp$snp[which.min(max_snp$p)]))
## highlight lead SNP from 1kg
leadSNP <- gwas_ld %>% 
  filter(grepl(paste(mostSigRSID), snp)) %>%
  dplyr::select(snp)

chr11_manhattanPlot_ld <- plotManhattan(data = gwas_ld, params = region, chromstart = start_reg,
                                        chromend = end_reg,
                                        fill = c("grey", "#37a7db"),
                                        # sigLine = TRUE, sigCol = "orange",
                                        leadSNP = list(snp = leadSNP$snp,
                                                       pch = 18,
                                                       cex = 0.75,
                                                       fill = "red",
                                                       fontsize = 8),
                                        x = 0.75, y = 3.6, range = c(0,maxVal),
                                        width = 6, height = 1.9,
                                        just = c("left", "top"),
                                        default.units = "inches")

## Annotate genome label
annoGenomeLabel(plot = chr11_manhattanPlot_ld, x = 0.75, y = 5.5,
                fontsize = 8, scale = "Mb",
                just = c("left", "top"), default.units = "inches")

## Annotate y-axis
annoYaxis(plot = chr11_manhattanPlot_ld,
          at = c(0, Yaxis * 0.2, Yaxis * 0.4, Yaxis * 0.6, Yaxis * 0.8, Yaxis),
          axisLine = TRUE, fontsize = 8)

## Plot y-axis label
plotText(label = "-log10(p-value)", x = 0.3, y = 4.55, rot = 90,
         fontsize = 8, fontface = "bold", just = "center",
         default.units = "inches")

# pagePlotRemove(chr11_manhattanPlot)
#for (i in 1:nrow(dnaLoops)) {
#  annoHighlight(
 #   plot = archPlot,
  #  chrom = chromo,
   # chromstart = dnaLoops[i,]$start1, chromend = dnaLoops[i,]$end1,
#    y = 2.9, height = 2.6, just = c("left", "top"),
 #   default.units = "inches"
  #)
#}

#for (i in 1:nrow(dnaLoops)) {
 # annoHighlight(
  #  plot = archPlot,
   # chrom = chromo,
    #chromstart = dnaLoops[i,]$start2, chromend = dnaLoops[i,]$end2,
#    y = 2.9, height = 2.6, just = c("left", "top"),
 #   default.units = "inches"
  #)
#}

print(gene_name)
}
else{
  pageCreate(width = 7, height = 5.75, default.units = "inches", showGuides = FALSE)


region <- pgParams(chrom = chromo, 
                   assembly = "hg19",
                   chromstart = start_reg, 
                   chromend = end_reg)

## add genes track within selected region and highlight gene of interest
genes_a <- plotGenes(params = region, 
                     stroke = 1, 
                     fontsize = 6, 
                     assembly = "hg19",
                     y = 0, 
                     height = 1.0, x = 0.75, width = 6, 
                     just = c("left", "top"),
                     geneHighlights = data.frame("gene" = c(gene_name), "color" = c("#225EA8")),
                     geneBackground = "grey", default.units = "inches") 

# pagePlotRemove(genes_a)

## add bp chromosome position label at bottom
annoGenomeLabel(plot = genes_a, 
                x = 0.75, y = 0.95, #4.6, 
                fontsize = 8, 
                just = c("left", "top"), 
                default.units = "inches")

h3k27ac_peaks <- read_tsv(file = paste0("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/human_PU1nuclei_H3K27ac_epilepsy_pooled_hg19.", chromo, ".ucsc.bed"), col_names = FALSE) %>% unique(); colnames(h3k27ac_peaks) = c("chrom", "start", "end", "peak_id", "score")

nott_signal_nonoverlap <- h3k27ac_peaks %>%
  filter(start >= start_reg & end <= end_reg) %>%
  unique() %>%
  mutate(end_original = end, end = end-0.1)

signal_bar_plot <- data.frame()

for (i in seq(from = start_reg, to = end_reg, by = diff)) {
  small_df <- subset(nott_signal_nonoverlap, start > i & end < (i+diff))
  max_index <- which.max(small_df$score)
  small_df <- small_df[max_index,]
  signal_bar_plot <- rbind(signal_bar_plot, small_df)
}

signal_bar_plot$score <- ifelse(signal_bar_plot$score < 4, 0, signal_bar_plot$score)

signal <- ggplot(data = signal_bar_plot, aes(x = start, y = score)) +
  geom_bar(width = wid, stat = "identity", fill = "blue") +
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        plot.margin = unit(c(0,0,0,0),"cm"))

plotGG(
  plot = signal, x = 0.43, y = 1.35, width = 6.64, height = 0.65,
  just = c("left", "top"), default.units = "inches"
)

plotSegments(
  x0 = 0.75, y0 = 1.93, x1 = 6.75, y1 = 1.93,
  default.units = "inches",
  lwd = 1, lty = 1, linecolor = "blue"
)

total_enhancers = read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/nott_microglial_enhancers.tsv", col_names = F)
total_enhancers <- subset(total_enhancers, total_enhancers$X1 == gene_chrom)

plotRanges(total_enhancers,
           chrom = chromo,
           assembly = "hg19",
           fill = "grey",
           boxHeight = unit(5, "mm"),
           x = 0.75, y = 2.15,
           width = 6, height = 0.5,
           just = c("left", "top"),
           params = region)

plotText(label = "Microglial H3K27Ac ChIP-seq peaks - Nott et al (2019)", 
         fontsize = 8, 
         fontcolor = c("black", "#37a7db"),
         x = 0.75, y = 1.25, 
         just = "left, top", 
         default.units = "inches")

plotText(label = "Microglia promoter/enhancer elements - Nott et al (2019)", 
         fontsize = 8, 
         fontcolor = c("black", "#37a7db"),
         x = 0.75, y = 2.05, 
         just = "left, top", 
         default.units = "inches")
  print(gene_name)

}
}

dev.off()
