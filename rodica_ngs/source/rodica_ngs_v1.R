# |----------------------------------------------------------------------------------|
# | Project: Rodica's NGS                                                            |
# | Script: DNA vs. RNA                                                              |
# | Author: Davit Sargsyan                                                           |
# | Created: 08/25/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_rodica_ngs_v1.txt")

# Header----
require(data.table)
require(ggplot2)
require(readxl)

# Read data----
# RNA----
rna <- fread("rodica_ngs/data/Heatmap Output (3).csv")[,  1:6]
names(rna)[1] <- "gene"
rna

# DNA----
dna <- read_xlsx(path = "rodica_ngs/data/combined_fr5c5_Yen_anno.xlsx",
                 sheet = 1)
dna <- data.table(dna)
dna$`11.Control` <- as.numeric(as.character(dna$`11.Control`))
dna$`12.UA` <- as.numeric(as.character(dna$`12.UA`))
dna$`13.SFN` <- as.numeric(as.character(dna$`13.SFN`))
dna$`14.RA` <- as.numeric(as.character(dna$`14.RA`))
dna$`15.UA+RA` <- as.numeric(as.character(dna$`15.UA+RA`))
dna$`16.SFN+RA` <- as.numeric(as.character(dna$`16.SFN+RA`))

# In DNA data, keep only the genes selected from RNA expressions----
dna <- dna[gene %in% rna$gene, ]
dna

length(unique(dna$gene))
# 99

length(unique(rna$gene))
# 109

# Differences----
dna$diff_UA_Ctrl <- 100*(dna$`12.UA` - dna$`11.Control`)
dna$mean_UA_Ctrl <- 100*(dna$`12.UA` + dna$`11.Control`)/2
# dna <- dna[!is.na(dna$diff_UA_Ctrl), ]

plot(dna$diff_UA_Ctrl ~ dna$mean_UA_Ctrl)

# Separate genes with meaningfull (>10%) differences----
gene.keep <- unique(dna$gene[abs(dna$diff_UA_Ctrl) >= 10])
gene.keep
# 25 genes

dna <- dna[gene %in% gene.keep, ]
dna

# Merge with RNA----
dt1 <- merge(dna,
             rna,
             by = "gene")
dt1$feature[substr(dt1$feature, 1, 4) == "Exon"] <- "Exon"
dt1$feature[substr(dt1$feature, 1, 6) == "Intron"] <- "Intron"
dt1$feature[substr(dt1$feature, 1, 8) == "Promoter"] <- "Promoter"
dt1$feature <- factor(dt1$feature)
table(dt1$feature)

summary(dt1$CpG)

dt1$cpg.grp <- "5 to 10"
dt1$cpg.grp[dt1$CpG > 10] <- "11 to 20"
dt1$cpg.grp[dt1$CpG > 20] <- ">20"
dt1$cpg.grp <- factor(dt1$cpg.grp,
                      levels = c("5 to 10",
                                 "11 to 20",
                                 ">20"))

# Isolate genes----
for (i in 1:length(unique(dna$gene))) {
  gX <- unique(dt1$gene)[i]
  dna.gX <- dt1[dt1$gene %in% gX, ]
  dna.gX[, x := rev(rank(distance)),
         by = gene]
  dna.gX$y0 <- 0
  
  p1 <- ggplot(dna.gX,
               aes(x = x,
                   y = diff_UA_Ctrl)) +
    facet_wrap(.~ gene,
               scales = "free_x",
               nrow = 1) +
    geom_hline(yintercept = c(10, 
                              0, 
                              -10),
               linetype = c("dashed", 
                            "solid", 
                            "dashed")) +
    geom_segment(aes(x = x,
                     y = y0,
                     xend = x,
                     yend = diff_UA_Ctrl)) + 
    geom_point(aes(x = x,
                   y = diff_UA_Ctrl,
                   fill = feature,
                   size = cpg.grp),
               shape = 21) +
    ggtitle(paste("Gene: ",
                  gX,
                  "\nRNA log2 Difference =",
                  round(dna.gX$UA.vs.C.1[1],
                        3))) +
    scale_x_continuous("Distance from TSS",
                       breaks = dna.gX$x,
                       labels = dna.gX$distance) +
    scale_y_continuous("UA - Control, % Methylation") +
    scale_fill_manual("Region",
                      values = c("Distal Intergenic" = "red",
                                 "Exon" = "blue",
                                 "Intron" = "white",
                                 "Promoter" = "green",
                                 "3' UTR" = "black",
                                 "5' UTR" = "yellow")) +
    scale_size_manual("Number of CpG-s",
                      values = c("5 to 10" = 5,
                                 "11 to 20" = 6,
                                 ">20" = 7)) +
    guides(fill = guide_legend(override.aes = list(size = 7))) +
    theme(plot.title = element_text(hjust = 0.5),
          legend.position = "top",
          axis.text.x = element_text(angle = 45,
                                     hjust = 1))
  
  tiff(filename = paste("tmp/",
                        gX,
                        ".tiff",
                        sep = ""),
       height = 5,
       width = 10,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  print(p1)
  graphics.off()
}

# sessionInfo()
# sink()