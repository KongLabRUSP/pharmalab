# Project: MethylSeq Data Gene Selection, Yue's Old Samples
# Author: Davit Sargsyan
# Created: 04/26/2017
# Source: http://stackoverflow.com/questions/8713994/venn-diagram-proportional-and-color-shading-with-semi-transparency
#*****************************************************
getwd()
require(data.table)
require(VennDiagram)
require(gridExtra)

#*****************************************************
# Load gene data----
genes <- fread("methylation seq/data/genes.csv", header = FALSE)

dt.g <- subset(genes, select = c(1, 3:5, 9))
names(dt.g) <- c("chr", "reg", "start", "stop", "gene")

setkey(dt.g, start)
lvls <- unique(dt.g$chr)[order(as.numeric(substr(unique(dt.g$chr), 4, 6)))]
dt.g$chr <- factor(dt.g$chr,
                   levels = lvls)

dt.g$reg <- factor(dt.g$reg,
                   levels = c("exon",
                              "CDS",
                              "start_codon",
                              "stop_codon"))
levels(dt.g$reg)

dt.g$gene <- substr(dt.g$gene, 11, nchar(dt.g$gene) - 2)
dt.g$gene <- factor(dt.g$gene, levels = unique(dt.g$gene))

# Remove duplicates
dt.g <- unique(dt.g)
dt.g

# Start and end of the genes
dt.g$start
dt.g$stop

dt.g[, gstart := min(start, stop), 
     by = list(gene, chr)]
dt.g[, gstop := max(start, stop), 
     by = list(gene, chr)]
setkey(dt.g, gstart)
dt.g

gg <- unique(subset(dt.g,
                    select = c(1, 5:7)))
setkey(gg, gstart, gene)
gg$gene <- as.character(gg$gene)
gg
levels(gg$chr)
gg$chr <- as.character(gg$chr)

# Include extra 20kb on each side of the gene
extra.start = 20000
extra.stop = 20000

#*****************************************************
# Load Methyl-seq Data----
# 1. Negative Control vs. AOM+DSS
dt1 <- fread("D:/Yue_MethylSeq_Processed_03202017/Results/results_yue_old_2.ctrl_vs_2.aom.dss_all.p_mindiff.0.1.csv")
dt1[, p.fdr := p.adjust(p = `Control->Exptl:pval`,
                        method = "fdr")]
setkey(dt1)
dt1

# 2. AOM+DSS vs. AOM+DSS+Cur
dt2 <- fread("D:/Yue_MethylSeq_Processed_03202017/Results/results_yue_old_2.aom.dss_vs_2.aom.dss.curc_all.p_mindiff.0.1.csv")
dt2[, p.fdr := p.adjust(p = `Control->Exptl:pval`,
                        method = "fdr")]
setkey(dt2)
dt2

#*****************************************************
# Venn Diagrams----
# NOTE: 
# in Contrast1: diff = (AOM+DSS) - (Negative Control)
# in Contrast2: diff = (AOM+DSS+Cur) - (AOM+DSS)

# 1a. Upregulation by (AOM+DSS)---- 
aom.dss.up <- dt1[`Control->Exptl:diff` < 0 & 
                    p.fdr <= 0.01,
                  c(1:4, 7, 9), with = FALSE]
# Map genes
aom.dss.up$gene <- NA
l1 <- list()
for(i in 1:nrow(aom.dss.up)) {
  out <- gg[chr == aom.dss.up$chr[i] &
              gstart <= (min(aom.dss.up$start[i], aom.dss.up$end[i]) + extra.start) & 
              gstop >= (max(aom.dss.up$start[i], aom.dss.up$end[i])) - extra.stop, ]
  if (nrow(out) > 0){
    l1[[i]] <- out$gene
    aom.dss.up$gene[i] <- paste(out$gene, collapse = "***")
  }  
}
aom.dss.up

# List of  upregulated genes
l1 <- unique(do.call("c", l1))
l1 <- l1[order(l1)]
l1

# 1b. Downregulation by (AOM+DSS+Cur)----
aom.dss.cur.down <- dt2[`Control->Exptl:diff` > 0 & 
                          p.fdr <= 0.01,
                        c(1:4, 7, 9), with = FALSE]
# Map genes
aom.dss.cur.down$gene <- NA
l2 <- list()
for(i in 1:nrow(aom.dss.cur.down)) {
  out <- gg[chr == aom.dss.cur.down$chr[i] &
              gstart <= (min(aom.dss.cur.down$start[i], aom.dss.cur.down$end[i]) + extra.start) & 
              gstop >= (max(aom.dss.cur.down$start[i], aom.dss.cur.down$end[i])) - extra.stop, ]
  if (nrow(out) > 0){
    l2[[i]] <- out$gene
    aom.dss.cur.down$gene[i] <- paste(out$gene, collapse = "***")
  }  
}
aom.dss.cur.down

# List of  upregulated genes
l2 <- unique(do.call("c", l2))
l2 <- l2[order(l2)]
l2

# 1c. Venn Diagram: Upregulation by (AOM+DSS) and Downregulation by (AOM+DSS+Cur)----
xx <- unique(c(l1, l2))
lbls <- l1[l1 %in% l2]
lbls <- lbls[order(lbls)]

tiff(filename = "tmp/venn_yue_old_up_down.tiff",
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
                   main = "Yue's Old Samples, Two Per Treatment",
                   sub = "A (Upregulation): (AOM+DSS) < (Negative Control)\nB (Downlregulation): (AOM+DSS+Cur) > (AOM+DSS)")
grid.arrange(gTree(children = p1),
             tableGrob(matrix(c(lbls,
                                rep("",
                                    75 - length(lbls))),
                              ncol = 5)),
             nrow = 1)
graphics.off()

#*****************************************************
# 2a. Downregulation by (AOM+DSS)---- 
aom.dss.down <- dt1[`Control->Exptl:diff` > 0 & 
                      p.fdr <= 0.01,
                    c(1:4, 7, 9), with = FALSE]
# Map genes
aom.dss.down$gene <- NA
l4 <- list()
for(i in 1:nrow(aom.dss.down)) {
  out <- gg[chr == aom.dss.down$chr[i] &
              gstart <= (min(aom.dss.down$start[i], aom.dss.down$end[i]) + extra.start) & 
              gstop >= (max(aom.dss.down$start[i], aom.dss.down$end[i])) - extra.stop, ]
  if (nrow(out) > 0){
    l4[[i]] <- out$gene
    aom.dss.down$gene[i] <- paste(out$gene, collapse = "***")
  }  
}
aom.dss.down

# List of  upregulated genes
l4 <- unique(do.call("c", l4))
l4 <- l4[order(l4)]
l4

# 2b. Upregulation by (AOM+DSS+Cur)----
aom.dss.cur.up <- dt2[`Control->Exptl:diff` < 0 & 
                        p.fdr <= 0.01,
                      c(1:4, 7, 9), with = FALSE]
# Map genes
aom.dss.cur.up$gene <- NA
l5 <- list()
for(i in 1:nrow(aom.dss.cur.up)) {
  out <- gg[chr == aom.dss.cur.up$chr[i] &
              gstart <= (min(aom.dss.cur.up$start[i], aom.dss.cur.up$end[i]) + extra.start) & 
              gstop >= (max(aom.dss.cur.up$start[i], aom.dss.cur.up$end[i])) - extra.stop, ]
  if (nrow(out) > 0){
    l5[[i]] <- out$gene
    aom.dss.cur.up$gene[i] <- paste(out$gene, collapse = "***")
  }  
}
aom.dss.cur.up

# List of  upregulated genes
l5 <- unique(do.call("c", l5))
l5 <- l5[order(l5)]
l5

# 2c. Venn Diagram: Downregulation by (AOM+DSS) and Upregulation by (AOM+DSS+Cur)----
xx <- unique(c(l4, l5))
lbls <- l4[l4 %in% l5]
lbls <- lbls[order(lbls)]

tiff(filename = "tmp/venn_old_down_up.tiff",
     height = 5,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
p1 <- venn.diagram(x = list(A = which(xx %in% l4),
                            B = which(xx %in% l5)),
                   filename = NULL,
                   fill = c("green", "red"), 
                   alpha = c(0.5, 0.5),
                   compression = "lzw+p",
                   main = "Yue's old Samples, Two Per Treatment",
                   sub = "A (Downlregulation): (AOM+DSS) > (Negative Control)\nB (Upregulation): (AOM+DSS+Cur) < (AOM+DSS)")
grid.arrange(gTree(children = p1),
             tableGrob(matrix(c(lbls,
                                rep("", 3)),
                              ncol = 4)),
             nrow = 1)
graphics.off()

l1[l1 == "Tnf"]
l2[l2 == "Tnf"]

l4[l4 == "Tnf"]
l5[l5 == "Tnf"]

#*****************************************************
#*****************************************************
# OLD CODE----
# 2. Positive Contrast1 + Negative Contrast2---- 
# i.e. upregulation by AOM+DSS and downregulation by AOM+DOW+Cur
dd.down <- merge(dt1[`Control->Exptl:diff` > 0,
                     c(1:4, 7, 9), with = FALSE],
                 dt2[`Control->Exptl:diff` < 0,
                     c(1:4, 7, 9), with = FALSE],
                 by = names(dt1)[1:4],
                 all = FALSE)
dd.down
dd.down.sig <- subset(dd.down, 
                      p.fdr.y <= 0.05 & 
                        p.fdr.x <= 0.05)
dd.down.sig

# Map genes
dd.down.sig$gene <- NA
ldown <- list()
for(i in 1:nrow(dd.down.sig)) {
  out <- gg[gg$chr == dd.down.sig$chr[i] &
              gstart <= (min(dd.down.sig$start[i], dd.down.sig$end[i]) + extra.start) & 
              gstop >= (max(dd.down.sig$start[i], dd.down.sig$end[i])) - extra.stop, ]
  if (nrow(out) > 0){
    ldown[[i]] <- out$gene
    dd.down.sig$gene[i] <- paste(out$gene, collapse = "***")
  }  
}
dd.down.sig

# List of  downregulated genes
ldown <- unique(do.call("c", ldown))
ldown <- ldown[order(ldown)]
ldown

#*****************************************************
# Save----
dd.up.sig
lup
lup[substr(lup, 1, 1) == "T"]

dd.down.sig
ldown

out <- list(dd.up.sig = dd.up.sig,
            lup = lup,
            dd.down.sig = dd.down.sig,
            ldown = ldown)

save(out,
     file = "tmp/yue_old_genes.RData")
