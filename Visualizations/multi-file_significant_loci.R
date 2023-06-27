
df.magma.gwas <- read.table("/Users/aryanpradhan/Downloads/large_table_output_chrall.genes.out", header = TRUE)
#df.magma.gwas <- separate(df.magma.gwas, GENE, into = c("ENSG_ID", "GENE"), sep = ":")

df.magma.gwas <- df.magma.gwas %>%
  mutate(overall_p_value = ifelse(num_significant > 10 & overall_p_value < 1e-300, 1e-300, overall_p_value))


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


stats <- stats %>%
  mutate(BP_REGION = cut(BP, breaks = seq(0, max(BP), by = 5e6), include.lowest = TRUE))

stats <- stats %>%
  group_by(CHR, BP_REGION) %>%
  arrange(P, desc(num_significant)) %>%
  summarise(
    P = first(P),  # smallest P-value in the region
    SNP = first(SNP),  # SNP corresponding to smallest P-value
    BP = mean(BP),  # average BP value to represent the region
    sum_num_significant = sum(num_significant, na.rm = TRUE)  # sum of num_significant for the region
  )

stats <- stats %>%
  mutate(sum_num_significant = ifelse(sum_num_significant > 225, 225, sum_num_significant))

stats$p_num <- 10^(-stats$sum_num_significant)




manhattan <- stats %>%
  ggman(snp = "SNP", bp = "BP", chrom = "CHR", pvalue = "p_num",
        relative.positions = TRUE, sigLine = 5,
        lineColour = "blue") +
  geom_hline(yintercept = 20, color = "red") +
  geom_point(aes(color=as.factor(CHR), size = sum_num_significant), alpha=0.5) +
  geom_point(aes(color=as.factor(CHR), size = sum_num_significant), alpha=1, size=2) + # Smaller opaque points
  scale_color_manual(values = rep(c("blue", "skyblue"), 22 )) +
  scale_size_continuous(range = c(1, 1)) + # adjust size range as needed
  coord_cartesian(ylim = c(0, 250)) +
  ylab("# Significant Hits in Locus") +  
  theme_bw() + # customize theme
  ggtitle("Top Significant Loci (225 Cap)") +  # Add a title
  theme(legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

manhattan

stats <- subset(stats, sum_num_significant > 19)

ggmanLabel(manhattan, labelDfm = as.data.frame(stats), snp = "SNP", label = "SNP")
