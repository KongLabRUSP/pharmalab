# Project: MethylSeq Data Gene Selection, Yue's Old Samples
# Author: Davit Sargsyan
# Created: 04/26/2017
#*****************************************************
getwd()
require(data.table)

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
dt1 <- fread("D:/David_MethylSeq_Processed_04042017/Results/results_yue_old_2.ctrl_vs_2.aom.dss_all.p_mindiff.0.1.csv")
dt1[, p.fdr := p.adjust(p = `Control->Exptl:pval`,
                        method = "fdr")]
setkey(dt1)
dt1

# NOTE: 
# in Contrast1: diff = (AOM+DSS) - (Negative Control)
# in Contrast2: diff = (AOM+DSS+Cur) - (AOM+DSS)

# 1. Negative Contrast1 + Positive Contrast2---- 
# i.e. downregulation by AOM+DSS and upregulation by AOM+DOW+Cur
dd.up <- merge(dt1[`Control->Exptl:diff` < 0,
                   c(1:4, 7, 9), with = FALSE],
               dt2[`Control->Exptl:diff` > 0,
                   c(1:4, 7, 9), with = FALSE],
               by = names(dt1)[1:4],
               all = FALSE)
dd.up
dd.up.sig <- subset(dd.up, 
                    p.fdr.y <= 0.05 & 
                      p.fdr.x <= 0.05)
dd.up.sig

# Map genes
dd.up.sig$gene <- NA
lup <- list()
for(i in 1:nrow(dd.up.sig)) {
  out <- gg[gg$chr == dd.up.sig$chr[i] &
              gstart <= (min(dd.up.sig$start[i], dd.up.sig$end[i]) + extra.start) & 
              gstop >= (max(dd.up.sig$start[i], dd.up.sig$end[i])) - extra.stop, ]
  if (nrow(out) > 0){
    lup[[i]] <- out$gene
    dd.up.sig$gene[i] <- paste(out$gene, collapse = "***")
  }  
}
dd.up.sig

# List of  upregulated genes
lup <- unique(do.call("c", lup))
lup <- lup[order(lup)]
lup

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
