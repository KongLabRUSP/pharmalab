# |---------------------------------------------------------------------------------|
# | Project:  Rodica's NGS                                                          |
# | Script:   DNA vs. RNA                                                           |
# | Author:   Davit Sargsyan                                                        |
# | Created:  08/25/2018                                                            |
# | Modified: 08/28/2018: plot al 5 comparisons at once; use updated data from Renyi| 
# |                       with distance by gene, not by transcript.                 |
# |---------------------------------------------------------------------------------|
sink(file = "tmp/log_rodica_ngs_v2.txt")

# Header----
require(data.table)
require(ggplot2)

# Read data----
# RNA----
rna <- fread("rodica_ngs/data/Heatmap Output (3).csv")[,  1:6]
names(rna)[1] <- "gene"
rna

# DNA----
dna <- fread("rodica_ngs/data/combined_fr5c5_Yen_raw_by_gene.csv")

# Rename variables----
colnames(dna)[colnames(dna) %in% c("X11",
                                   "X12",
                                   "X13",
                                   "X14",
                                   "X15",
                                   "X16",
                                   "geneId")] <- c("11.Control",
                                                "12.UA",
                                                "13.SFN",
                                                "14.RA",
                                                "15.UA+RA",
                                                "16.SFN+RA",
                                                "gene")

# In DNA data, keep only the genes selected from RNA expressions----
dna <- dna[gene %in% rna$gene, ]
dna

length(unique(dna$gene))
# 97

length(unique(rna$gene))
# 109

# Differences----
dna$RA.vs.C.1 <- 100*(dna$`14.RA` - dna$`11.Control`)

dna$SFN.vs.C.1 <- 100*(dna$`13.SFN` - dna$`11.Control`)
dna$RA.SFN.vs.RA.1 <- 100*(dna$`16.SFN+RA` - dna$`14.RA`)

dna$UA.vs.C.1 <- 100*(dna$`12.UA` - dna$`11.Control`)
dna$RA.UA.vs.RA.1 <- 100*(dna$`15.UA+RA` - dna$`14.RA`)

# Separate genes with meaningfull (>10%) differences----
gene.keep <- unique(dna$gene[abs(dna$RA.vs.C.1) >= 10 |
                               abs(dna$SFN.vs.C.1) >= 10 |
                               abs(dna$RA.SFN.vs.RA.1) >= 10 |
                               abs(dna$UA.vs.C.1) >= 10 |
                               abs(dna$RA.UA.vs.RA.1) >= 10])
gene.keep
# 51 genes

dna <- dna[gene %in% gene.keep, ]

dna[, distRank := rev(rank(distanceToTSS)),
    by = gene]

dna

# Transform to Long format----
dt1 <- melt.data.table(data = dna,
                       id.vars = c("gene",
                                   "CpG",
                                   "annotation",
                                   "distanceToTSS",
                                   "distRank"),
                       measure.vars = c("RA.vs.C.1",
                                        "SFN.vs.C.1",
                                        "RA.SFN.vs.RA.1",
                                        "UA.vs.C.1",
                                        "RA.UA.vs.RA.1"),
                       variable.name = "Treatment",
                       value.name = "DNA")
dt1$Treatment <- as.character(dt1$Treatment)

dt1$annotation[substr(dt1$annotation, 1, 4) == "Exon"] <- "Exon"
dt1$annotation[substr(dt1$annotation, 1, 6) == "Intron"] <- "Intron"
dt1$annotation[substr(dt1$annotation, 1, 8) == "Promoter"] <- "Promoter"
dt1$annotation[substr(dt1$annotation, 1, 4) == "Down"] <- "Downstream"
dt1$annotation <- factor(dt1$annotation)

summary(dt1$CpG)

dt1$reg <- "5 to 10"
dt1$reg[dt1$CpG > 10] <- "11 to 20"
dt1$reg[dt1$CpG > 20] <- ">20"
dt1$reg <- factor(dt1$reg,
                      levels = c("5 to 10",
                                 "11 to 20",
                                 ">20"))
summary(dt1)

# RNA data Long format----
dt2 <- melt.data.table(data = rna,
                       id.vars = "gene",
                       measure.vars = c("RA.vs.C.1",
                                        "SFN.vs.C.1",
                                        "RA.SFN.vs.RA.1",
                                        "UA.vs.C.1",
                                        "RA.UA.vs.RA.1"),
                       variable.name = "Treatment",
                       value.name = "RNA")
dt2$Treatment <- as.character(dt2$Treatment)
dt2

# Merge DNA with RNA----
dt1 <- merge(dt1,
             dt2,
             by = c("gene",
                    "Treatment"))
dt1

# Isolate genes----
for (i in 1:length(unique(dna$gene))) {
  gX <- unique(dt1$gene)[i]
  dna.gX <- dt1[dt1$gene %in% gX, ]
  dna.gX$y0 <- 0
  
  dna.gX$Treatment <- paste(dna.gX$Treatment,
                            " (RNA = ",
                            round(dna.gX$RNA, 3),
                            ")",
                            sep = "")
  
  p1 <- ggplot(dna.gX,
               aes(x = distRank,
                   y = DNA)) +
    facet_wrap(.~ Treatment,
               scales = "free_y",
               ncol = 1) +
    geom_rect(aes(xmin = -Inf,
                  xmax = Inf,
                  ymin = -Inf,
                  ymax = -10),
              fill = "red",
              alpha = 0.1) +
    geom_rect(aes(xmin = -Inf,
                  xmax = Inf,
                  ymin = 10,
                  ymax = Inf),
              fill = "green",
              alpha = 0.1) +
    geom_hline(yintercept = 0) +
    geom_hline(yintercept = c(10,
                              -10),
               linetype = "dashed") +
    geom_segment(aes(x = distRank,
                     y = y0,
                     xend = distRank,
                     yend = DNA)) + 
    geom_point(aes(x = distRank,
                   y = DNA,
                   fill = annotation,
                   size = reg),
               shape = 21) +
    # geom_rect(aes(xmin = -Inf,
    #               xmax = Inf,
    #               ymin = -10,
    #               ymax = 10),
    #           fill = "white",
    #           alpha = 0.1) +
    ggtitle(paste("Gene:",
                  gX)) +
    scale_x_continuous("Distance from TSS",
                       breaks = dna.gX$distRank,
                       labels = dna.gX$distanceToTSS) +
    scale_y_continuous("% Methylation") +
    scale_fill_manual("Region",
                      values = c("Distal Intergenic" = "purple",
                                 "Exon" = "blue",
                                 "Intron" = "white",
                                 "Promoter" = "brown",
                                 "3' UTR" = "black",
                                 "5' UTR" = "yellow",
                                 "Downstream" = "orange")) +
    scale_size_manual("Number of CpG-s",
                      values = c("5 to 10" = 5,
                                 "11 to 20" = 6,
                                 ">20" = 7)) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(plot.title = element_text(hjust = 0.5),
          #legend.position = "top",
          axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  p1
  tiff(filename = paste("tmp/",
                        gX,
                        ".tiff",
                        sep = ""),
       height = 8,
       width = 8,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  print(p1)
  graphics.off()
}

sessionInfo()
sink()
