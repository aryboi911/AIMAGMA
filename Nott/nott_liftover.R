input_data <- read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/placseq_nott_2019/nott2019_plaqseq_microglia_chromatin_interactions_tables5.tsv")
#Start1
gr <- GRanges(seqnames = input_data$chr1,
          ranges = IRanges(start = input_data$start1, end = input_data$start1),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$start1 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "start1")

#End1
gr <- GRanges(seqnames = input_data$chr1,
          ranges = IRanges(start = input_data$end1, end = input_data$end1),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$end1 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "end1")

#Start2
gr <- GRanges(seqnames = input_data$chr1,
          ranges = IRanges(start = input_data$start2, end = input_data$start2),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$start2 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "start2")

#End2
gr <- GRanges(seqnames = input_data$chr1,
          ranges = IRanges(start = input_data$end2, end = input_data$end2),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$end2 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "end2")

unique(input_data) %>% write_tsv("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/nott/nott2019_HG19_plaqseq_microglia_chromatin_interactions_tables5.tsv")




input_data <- read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/input/epig_nott2019/nott_2019.microglia.pu1.annot")

#start
input_data$chrom <- paste0("chr", input_data$chrom)
gr <- GRanges(seqnames = input_data$chrom,
          ranges = IRanges(start = input_data$start, end = input_data$start),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$start <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "start")

#end
gr <- GRanges(seqnames = input_data$chrom,
          ranges = IRanges(start = input_data$end, end = input_data$end),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$end <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "end")

#gene_loc
input_data <- input_data %>%
  separate(gene_loc, into = c("chr", "st", "en"), sep = ":")

input_data$st <- as.numeric(input_data$st)
input_data$en <- as.numeric(input_data$en)

#st
gr <- GRanges(seqnames = input_data$chrom,
          ranges = IRanges(start = input_data$st, end = input_data$st),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$st <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "st")

#en
gr <- GRanges(seqnames = input_data$chrom,
          ranges = IRanges(start = input_data$en, end = input_data$en),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$en <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "en")

input_data$gene_loc <- paste0(input_data$chr, ":", input_data$st, ":", input_data$en)
input_data <- subset(input_data, select = -c(chr,en,st))


unique(input_data) %>% write_tsv("/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/nott/nott_2019_HG19.microglia.pu1.annot")

input_data <- read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/nott_microglial_promoters.tsv", col_names = F)

gr <- GRanges(seqnames = input_data$X1,
          ranges = IRanges(start = input_data$X2, end = input_data$X2),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$X2 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "X2")

gr <- GRanges(seqnames = input_data$X1,
          ranges = IRanges(start = input_data$X3, end = input_data$X3),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$X3 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "X3")

unique(input_data) %>% write.table(file = "/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/nott/nott_microglial_promoters_HG19.tsv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)


input_data <- read_tsv("/sc/arion/projects/load/users/patelt16/projects/aimagma/20220714.aimagma/data/h3k27ac_nott_2019/nott_microglial_enhancers.tsv", col_names = F)

gr <- GRanges(seqnames = input_data$X1,
          ranges = IRanges(start = input_data$X2, end = input_data$X2),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$X2 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "X2")

gr <- GRanges(seqnames = input_data$X1,
          ranges = IRanges(start = input_data$X3, end = input_data$X3),
          strand = "*")
gr_lifted <- liftOver(gr, chain)
lifted_positions <- start(gr_lifted)

input_data$X3 <- as.integer(lifted_positions)
input_data <- na.omit(input_data, cols = "X3")

unique(input_data) %>% write.table(file = "/sc/arion/projects/load/users/pradha04/projects/aimagma/abc_magma/data/nott/nott_microglial_enhancers_HG19.tsv", sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

