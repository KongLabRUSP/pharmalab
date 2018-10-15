# |---------------------------------------------------------------------------------|
# | Project:  Rodica's NGS                                                          |
# | Script:   DNA vs. RNA                                                           |
# | Author:   Davit Sargsyan                                                        |
# | Created:  08/25/2018                                                            |
# | Modified: 08/28/2018: plot al 5 comparisons at once; use updated data from Renyi| 
# |                       with distance by gene, not by transcript.                 |
# |---------------------------------------------------------------------------------|
# sink(file = "tmp/log_rodica_ngs_v2.txt")

# Header----
require(data.table)
require(ggplot2)
require(ggdendro)

# Read data----
# RNA----
rna <- fread("rodica_ngs/data/Heatmap Output (3).csv")[,  1:6]
names(rna)[1] <- "gene"
rna

# Heatmap----
# Distances----
tmp <- as.matrix(rna[, -1])
rownames(tmp) <- rna$gene
geneDist <- dist(tmp)

# Make dendrogram data----
dhc <- as.dendrogram(hclust(d = geneDist),
                     horiz = TRUE)
ddata <- dendro_data(dhc, 
                     type = "rectangle")

# Segment data----
dtp1 <- segment(ddata)
# Hitmap data----
dtp2 <- melt.data.table(rna,
                        id.vars = 1,
                        measure.vars = 2:6,
                        variable.name = "Comparison",
                        value.name = "dRNA")
dtp2$gene <- factor(dtp2$gene,
                    levels = ddata$labels$label)
dtp2

p1 <- ggplot(data = dtp2) +
  coord_polar("y",
              start = 0,
              direction = -1) +
  geom_tile(aes(x =  as.numeric(Comparison),
                y = gene, 
                fill = dRNA),
            color = "white") +
  geom_text(data = dtp2[Comparison == "UA.vs.C.1", ],
            aes(x = rep(4.75,
                        nlevels(gene)),
                y = gene,
                angle = 90 + seq(from = 0,
                                 to = 360,
                                 length.out = nlevels(gene))[as.numeric(gene)],
                label = unique(gene)),
            hjust = 0) +
  scale_fill_gradient2(low = "red", 
                       high = "green", 
                       mid = "grey", 
                       midpoint = 0, 
                       name = "Legend") +
  scale_y_discrete("",
                   expand = c(0, 0)) +
  ggtitle("Delta RNA Expression in Differentially Expressed Genes") + 
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_segment(data = dtp1,
               aes(x = -sqrt(y) + 0.5,
                   y = x, 
                   xend = -sqrt(yend) + 0.5,
                   yend = xend),
               size = 1) 
p1

tiff(filename = "tmp/heatmap.tiff",
     height = 10,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# PCA----
tmp <- t(tmp)
m1 <- prcomp(tmp)
m1
summary(m1)

# Biplot while keep only the most important variables (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2
# Scores, i.e. points (df.u)
dt.scr <- data.table(m1$x[, choices])
# Add grouping variable
dt.scr$grp <- rownames(tmp)
dt.scr

# Loadings, i.e. arrows (df.v)
dt.rot <- as.data.frame(m1$rotation[, choices])
dt.rot$feat <- rownames(dt.rot)
dt.rot <- data.table(dt.rot)
dt.rot

dt.load <- melt.data.table(dt.rot,
                           id.vars = "feat",
                           measure.vars = 1:2,
                           variable.name = "pc",
                           value.name = "loading")
dt.load$feat <- factor(dt.load$feat,
                       levels = unique(dt.load$feat))
# Plot loadings
p0 <- ggplot(data = dt.load,
             aes(x = feat,
                 y = loading)) +
  facet_wrap(~ pc,
             nrow = 2) +
  geom_bar(stat = "identity") +
  ggtitle("PC Loadings") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p0

# Axis labels
u.axis.labs <- paste(colnames(dt.rot)[1:2], 
                     sprintf('(%0.1f%% explained var.)', 
                             100*m1$sdev[choices]^2/sum(m1$sdev^2)))
u.axis.labs

# Keep only a few variables with high loadings in PC1 (94% of variation)----
gene.rank<- rank(x = abs(dt.rot$PC1)) < 10

var.keep.ndx <- which(dt.rot$feat %in% dt.rot$feat[gene.rank])

p1 <- ggplot(data = dt.rot[var.keep.ndx,],
             aes(x = PC1,
                 y = PC2)) +
  # coord_equal() +
  geom_point(data = dt.scr,
             aes(fill = grp),
             shape = 21,
             size = 3,
             alpha = 0.5) +
  geom_segment(aes(x = 0,
                   y = 0,
                   xend = 400*PC1,
                   yend = 10*PC2),
               arrow = arrow(length = unit(1/2, 'picas')),
               color = "black") +
  geom_text(aes(x = 450*PC1,
                y = 11*PC2,
                label = dt.rot$feat[var.keep.ndx]),
            size = 3,
            hjust = 0.5) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_discrete(name = "Comparison") +
  # scale_fill_discrete(name = grp) +
  ggtitle("Biplot of Differentially Expressed Genes") +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 10))

tiff(filename = "tmp/pca_biplot.tiff",
     height = 5.5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

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

# Starburst plot----
for (i in 1:length(unique(dt1$Treatment))) {
  cmp <- unique(dt1$Treatment)[i]
  
  
  tmp <- dt1[dt1$Treatment == cmp, ]
  tmp
  
  g1 <- tmp[DNA >= 10 & 
              RNA <= -0.5 &
              annotation == "Promoter"]
  g1
  length(unique(g1$gene))
  
  g2 <- tmp[DNA <= -10 & 
              RNA >= 0.5 &
              annotation == "Promoter"]
  g2
  length(unique(g2$gene))
  
  p1 <- ggplot(data = tmp,
               aes(x = DNA,
                   y = RNA,
                   fill = annotation)) +
    geom_point(alpha = 0.7,
               size = 2,
               shape = 21) +
    geom_text(data = unique(tmp[gene %in% unique(g1$gene) &
                                  annotation == "Promoter",
                                c("gene",
                                  "annotation",
                                  "RNA")]),
              aes(x = 40,
                  y = RNA,
                  label = gene),
              color = "blue",
              size = 2) +
    geom_text(data = unique(tmp[gene %in% unique(g2$gene) &
                                  annotation == "Promoter", 
                                c("gene",
                                  "annotation",
                                  "RNA")]),
              aes(x = -40,
                  y = RNA,
                  label = gene),
              color = "blue",
              size = 2) +
    geom_hline(yintercept = c(-0.5, 0.5),
               linetype = "dashed") +
    geom_vline(xintercept = c(-10, 10),
               linetype = "dashed") +
    scale_x_continuous("DNA Methylation Difference(%)",
                       breaks = seq(-50, 30, 10)) +
    scale_y_continuous("RNA Expression Difference (log2)",
                       breaks = seq(-5, 10, 1)) +
    ggtitle(cmp) +
    # scale_fill_manual("Region",
    #                   values = c("Promoter" = "green",
    #                              "5' UTR" = "white",
    #                              "Body" = "blue",
    #                              "3' UTR" = "grey",
    #                              "Downstream" = "red")) +
    theme(plot.title = element_text(hjust = 0.5))
  p1
  tiff(filename = paste("tmp/",
                        cmp,
                        ".tiff",
                        sep = ""),
       height = 10,
       width = 10,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  print(p1)
  graphics.off()
}

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

# sessionInfo()
# sink()