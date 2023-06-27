
df.magma.gwas <- read.table("/Users/aryanpradhan/Downloads/large_table_output_chrall.genes (1).out", header = TRUE)
#df.magma.gwas <- separate(df.magma.gwas, GENE, into = c("ENSG_ID", "GENE"), sep = ":")

df.magma.gwas$sig <- -log10(df.magma.gwas$overall_p_value)

df.magma.gwas <- df.magma.gwas %>%
  mutate(overall_p_value = ifelse(num_significant > 10 & overall_p_value < 1e-320, 1e-320, overall_p_value))


## manhattan plot for all chroms
genomewide = -log10(0.05 / nrow(df.magma.gwas))
suggestive = -log10(0.1 / nrow(df.magma.gwas))

stats <- df.magma.gwas %>%
  dplyr::rename(SNP = external_gene_name, BP = start_position, P = "overall_p_value") %>%
  mutate(CHR = as.integer(chromosome_name), BP = as.integer(BP))

stats %>%
  dplyr::rename(chromosome = CHR, position = BP) %>%
  arrange(chromosome, position)# %>%
# write_tsv(paste0(fname, ".loczoom.tsv"))

stats <- stats[stats$P > 0, ]

max <- ceiling(-log10(min(stats$P)) * 0.11)*10

manhattan <- stats %>%
  ggman(snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "P",
        relative.positions = TRUE, sigLine = 5,
        lineColour = "blue") +
  geom_hline(yintercept = -log10(5e-8), color = "red") +
  geom_point(aes(color=as.factor(CHR), size = prop_significant), alpha=0.5) +
  geom_point(aes(color=as.factor(CHR), size = prop_significant), alpha=1, size=2) + # Smaller opaque points
  scale_color_manual(values = rep(c("blue", "skyblue"), 22 )) +
  scale_size_continuous(range = c(1, 15)) + # adjust size range as needed
  coord_cartesian(ylim = c(0, 350)) +
  theme_bw() + # customize theme
  ggtitle("Top Overall Genes") +  # Add a title
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

manhattan

## subset SNPs with -log10P > 8
stats.sig <- stats[-log10(stats$P) > 8, ]
stats.sig.top <- stats.sig %>% 
  group_by(CHR) %>%
  arrange(P, .by_group = TRUE) %>%
  top_n(-3, P) %>%
  filter(-log10(P) > 8)

stats.sig.top <- subset(stats.sig.top, num_significant > 5)
stats.sig.top <- subset(stats.sig.top, CHR != 19)

filtered_genes <- stats %>%
  arrange(CHR, P, desc(prop_significant), desc(num_significant)) %>%
  group_by(CHR) %>%
  slice(1)


stats.sig.top <- rbind(stats.sig.top, filtered_genes)
stats.sig.top <- unique(stats.sig.top)
#stats.sig.top <- subset(stats.sig.top, prop_significant > 0.8)

## select specific genes to highlight
# stats.sig.picalm = stats.sig %>% filter(grepl("PICALM", SNP))
# stats.sig.eed = stats.sig %>% filter(grepl("EED", SNP))
# stats.sig.top <- rbind(stats.sig.top, stats.sig.picalm, stats.sig.eed)

ggmanLabel(manhattan, labelDfm = as.data.frame(stats.sig.top), snp = "SNP", label = "SNP")
