# Gene data is from NCBI mice mm9 ~/Annotation/Genes/gene.gtf file
require(data.table)
library(ggplot2)
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

hist(dt.g$gstop - dt.g$gstart, 100)
summary(dt.g$gstop - dt.g$gstart)

dt.g[dt.g$gstop - dt.g$gstart > 1000000, ]

# CHECK: separate TNF
dt.tnf <- droplevels(subset(dt.g,
                            chr == "chr17" &
                              gene == "Tnf"))
dt.tnf

ggplot(dt.tnf) +
  geom_rect(aes(xmin = start,
                xmax = stop,
                ymin = -as.numeric(dt.tnf$reg)/4,
                ymax = as.numeric(dt.tnf$reg)/4,
                fill = reg,
                alpha = 0.5)) +
  geom_segment(aes(x = min(start[reg == "start_codon"],
                           stop[reg == "start_codon"]),
                   xend = max(start[reg == "stop_codon"],
                              stop[reg == "stop_codon"]), 
                   y = 0, 
                   yend = 0), 
               size = 2,
               arrow = arrow(length = unit(0.5, "cm"))) +
  geom_hline(yintercept = 0) +
  scale_x_continuous("Location") + 
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  ggtitle("Chromosome 17: Tnf")

#*****************************************************
# Load alligned read-outs----
# dt1 <- fread("methylation seq/data/results_yue_03312017.csv")
# dt1 <- fread("methylation seq/data/results_yue_old_samples_only_all.csv")

# dt1 <- fread("methylation seq/results/kong_meeting_041017/results_yue_2X2_all.p_mindiff.0.1.csv")
dt1 <- fread("methylation seq/results/kong_meeting_041017/results_david_2X2_all.p_mindiff.0.1.csv")
setkey(dt1, start)
dt1

gg <- unique(subset(dt.g,
             select = c(1, 5:7)))
setkey(gg, gstart, gene)
gg$gene <- as.character(gg$gene)
gg

unique(dt1$chr)
levels(gg$chr)

# If the DMR region is completely within the gene region, assign that gene name to the region
dt1$gene <- NA
 
gg$chr <- as.character(gg$chr)

extra.start = 20000
extra.stop = 20000
for(i in 1:nrow(dt1)) {
  out <- gg[gg$chr == dt1$chr[i] &
              gstart <= (min(dt1$start[i], dt1$end[i]) + extra.start) & 
              gstop >= (max(dt1$start[i], dt1$end[i])) - extra.stop, ]
  if (nrow(out) > 0) dt1$gene[i] <- paste(out$gene, collapse = "***")
}
dt1
sum(is.na(dt1$gene))/nrow(dt1)
dt1 <- subset(dt1, !is.na(gene))
write.csv(dt1, file = "tmp/dt1.csv")


subset(dt.g, chr == "chr17" & gene %in% c("Tnf", "Sult1c1"))

subset(dt.g, chr == "chr2" & gene %in% c("Mir684-1"))

#*****************************************************
# Plot----
tmp <- droplevels(subset(dt1, 
                         chr == levels(chr)[1]))
tmp$ystart <- rep(0, nrow(tmp))
tmp$yend <- rep(1, nrow(tmp))

tmp <- melt.data.table(tmp,
                       id.vars = c("chr", "start", "end", "CpG", "ystart", "yend"),
                       variable.name = "experiment",
                       value.name = "nreadns")
tmp$experiment <- as.character(tmp$experiment)

tmp$sample <- substr(tmp$experiment, 1, nchar(tmp$experiment) - 2)
tmp$sample <- factor(tmp$sample,
                     levels = unique(tmp$sample))

tmp$what <- substr(tmp$experiment, nchar(tmp$experiment), nchar(tmp$experiment))
tmp$what <- factor(tmp$what,
                   levels = c("X", "N"))
tmp.X <- subset(tmp, what == "X", select = -c(7, 10))
tmp.N <- subset(tmp, what == "N", select = -c(7, 10))
tmp <- merge(tmp.X, tmp.N, 
             by = c("chr", "start", "end", "CpG", "ystart", "yend", "sample"))
tmp$pct.meth <- 100*tmp$nreadns.x/tmp$nreadns.y
hist(tmp$pct.meth, 100)
tmp

ggplot(tmp) +
  facet_wrap(~ sample,
             ncol = 1) +
  geom_segment(aes(x = start,
                   xend = start,
                   y = ystart,
                   yend = pct.meth),
               col = "green") +
  scale_x_continuous("Starting Location of Regions") +
  scale_y_continuous("Percent CpG Methylated") +
  ggtitle(paste("Chromosome:", tmp$chr[1]))

ggplot(tmp) +
  facet_wrap(~ sample,
             ncol = 1) +
  geom_rect(aes(xmin = start,
                xmax = end,
                ymin = ystart,
                ymax = pct.meth),
            col = "green") +
  scale_x_continuous("Starting Location of Regions") +
  scale_y_continuous("Percent CpG Methylated") +
  ggtitle(paste("Chromosome:", tmp$chr[1]))


#*****************************************************
# Genes----
genes <- fread("data/genes.csv", header = FALSE)
genes

dt.g <- subset(genes, select = c(1, 3:5, 9))
names(dt.g) <- c("chr", "reg", "start", "stop", "gene")

setkey(dt.g, start)
lvls <- unique(dt.g$chr)[order(as.numeric(substr(unique(dt.g$chr), 4, 6)))]
dt.g$chr <- factor(dt.g$chr,
                   levels = lvls)

dt.g$reg <- factor(dt.g$reg)

dt.g$gene <- substr(dt.g$gene, 11, nchar(dt.g$gene) - 2)
dt.g$gene <- factor(dt.g$gene, levels = unique(dt.g$gene))

summary(dt.g)