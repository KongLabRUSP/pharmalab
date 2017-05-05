# Project: MethylSeq Data Gene Selection, Yue's New Samples
# Author: Davit Sargsyan
# Created: 04/26/2017
# Source: http://stackoverflow.com/questions/8713994/venn-diagram-proportional-and-color-shading-with-semi-transparency
#*****************************************************
require(data.table)
require(VennDiagram)
require(gridExtra)

# Source: file:///C:/R/R-3.3.2/library/ChIPseeker/doc/ChIPseeker.html
# Source: http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html

# source("https://bioconductor.org/biocLite.R")
# biocLite("ChIPseeker")
# biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
# biocLite("org.Mm.eg.db")
# biocLite("ReactomePA")

require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
require(ReactomePA)

# Parameters----
minCpG <- 5
nLabCol <- 4

#*****************************************************
# Load Methyl-seq Data----
# 1. Negative Control vs. AOM+DSS
f1 <- "D:/Yue_MethylSeq_Processed_03202017/Results/results_yue_new_2.ctrl_vs_2.aom.dss_all.p_mindiff.0.1.csv"
peakAnno1 <- annotatePeak(peak = f1, 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                          annoDb = "org.Mm.eg.db")
dt1 <- data.table(as.data.frame(peakAnno1@anno@elementMetadata@listData))

dt1[, p.fdr := p.adjust(p = Control..Exptl.pval,
                        method = "fdr")]
setkey(dt1)
dt1

# 2. AOM+DSS vs. AOM+DSS+Cur
f2 <- "D:/Yue_MethylSeq_Processed_03202017/Results/results_yue_new_2.aom.dss_vs_2.aom.dss.curc_all.p_mindiff.0.1.csv"
peakAnno2 <- annotatePeak(peak = f2, 
                          tssRegion = c(-3000, 3000), 
                          TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                          annoDb = "org.Mm.eg.db")
dt2 <- data.table(as.data.frame(peakAnno2@anno@elementMetadata@listData))
dt2

dt2[, p.fdr := p.adjust(p = Control..Exptl.pval,
                        method = "fdr")]
setkey(dt2)
dt2

#*****************************************************
# Venn Diagrams----
# NOTE: 
# in Contrast1: diff = (AOM+DSS) - (Negative Control)
# in Contrast2: diff = (AOM+DSS+Cur) - (AOM+DSS)

# 1a. Upregulation by (AOM+DSS)---- 
aom.dss.up <- dt1[Control..Exptl.diff < 0 & 
                    p.fdr <= 0.01 &
                    CpG >= minCpG, ]

# List of  upregulated genes
l1 <- as.character(unique(aom.dss.up$SYMBOL))
l1 <- l1[!is.na(l1)]
l1 <- l1[order(l1)]
l1

# 1b. Downregulation by (AOM+DSS+Cur)----
aom.dss.cur.down <- dt2[Control..Exptl.diff > 0 & 
                          p.fdr <= 0.01 &
                          CpG >= minCpG,]

# List of  upregulated genes
l2 <- as.character(unique(aom.dss.cur.down$SYMBOL))
l2 <- l2[!is.na(l2)]
l2 <- l2[order(l2)]
l2

# 1c. Venn Diagram: Upregulation by (AOM+DSS) and Downregulation by (AOM+DSS+Cur)----
xx <- unique(c(l1, l2))
lbls <- l1[l1 %in% l2]
lbls <- lbls[order(lbls)]
if (floor(length(lbls)/nLabCol) < length(lbls)/nLabCol) {
  lbls <- c(lbls,
            rep("", 
                nLabCol*ceiling(length(lbls)/nLabCol) - length(lbls))) 
}
lbls <- matrix(lbls,
               ncol = nLabCol)

tiff(filename = "tmp/venn_yue_new_up_down_chipseeker_cpg5over.tiff",
     height = 5,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
p1 <- venn.diagram(x = list(A = which(xx %in% l1),
                            B = which(xx %in% l2)),
                   filename = NULL,
                   fill = c("red", "green"),
                   alpha = c(0.5, 0.5),
                   compression = "lzw+p",
                   main = "Yue's New Samples, Two Per Treatment, CpG > 4",
                   sub = "A (Upregulation): (AOM+DSS) < (Negative Control)\nB (Downlregulation): (AOM+DSS+Cur) > (AOM+DSS)")
grid.arrange(gTree(children = p1),
             #tableGrob(lbls),
             nrow = 1)
graphics.off()

#*****************************************************
# 2a. Downregulation by (AOM+DSS)---- 
aom.dss.down <- dt1[Control..Exptl.diff > 0 & 
                      p.fdr <= 0.01 &
                      CpG >= minCpG, ]

# List of  upregulated genes
l3 <- as.character(unique(aom.dss.down$SYMBOL))
l3 <- l3[!is.na(l3)]
l3 <- l3[order(l3)]
l3

# 2b. Upregulation by (AOM+DSS+Cur)----
aom.dss.cur.up <- dt2[Control..Exptl.diff < 0 & 
                        p.fdr <= 0.01 &
                        CpG >= minCpG, ]

# List of  upregulated genes
l4 <- as.character(unique(aom.dss.cur.up$SYMBOL))
l4 <- l4[!is.na(l4)]
l4 <- l4[order(l4)]
l4

# 2c. Venn Diagram: Downregulation by (AOM+DSS) and Upregulation by (AOM+DSS+Cur)----
xx <- unique(c(l3, l4))
lbls <- l3[l3 %in% l4]
lbls <- lbls[order(lbls)]
if (floor(length(lbls)/nLabCol) < length(lbls)/nLabCol) {
  lbls <- c(lbls,
            rep("", 
                nLabCol*ceiling(length(lbls)/nLabCol) - length(lbls))) 
}
lbls <- matrix(lbls,
               ncol = nLabCol)

tiff(filename = "tmp/venn_yue_new_down_up_chipseeker_cpg5over.tiff",
     height = 5,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
p1 <- venn.diagram(x = list(A = which(xx %in% l3),
                            B = which(xx %in% l4)),
                   filename = NULL,
                   fill = c("green", "red"), 
                   alpha = c(0.5, 0.5),
                   compression = "lzw+p",
                   main = "Yue's New Samples, Two Per Treatment",
                   sub = "A (Downlregulation): (AOM+DSS) > (Negative Control)\nB (Upregulation): (AOM+DSS+Cur) < (AOM+DSS)")
grid.arrange(gTree(children = p1),
             #tableGrob(lbls),
             nrow = 1)
graphics.off()

# CHECKPOINT
l1[l1 == "Tnf"]
l2[l2 == "Tnf"]

l3[l3 == "Tnf"]
l4[l4 == "Tnf"]

#*****************************************************
# Save----
write.csv(dt1, 
          file = "tmp/negctr_aomdss_new.csv",
          row.names = FALSE)
write.csv(dt2, 
          file = "tmp/aomdss_aomdsscur_new.csv",
          row.names = FALSE)
out <- list(negctr_aomdss_new = dt1,
            aomdss_aomdsscur_new = dt2)
save(out,
     file = "tmp/yue_new_genes.RData")