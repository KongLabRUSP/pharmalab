# Gene ddata is from NCBI mice mm9 ~/Annotation/Genes/gene.gtf file

require(data.table)
library(ggplot2)
library(Cairo)
#*****************************************************
genes <- fread("methylation seq/data/genes.csv", header = FALSE)

dt.g <- subset(genes, select = c(1, 3:5, 9))
names(dt.g) <- c("chr", "reg", "start", "stop", "gene")

setkey(dt.g, start)
lvls <- unique(dt.g$chr)[order(as.numeric(substr(unique(dt.g$chr), 4, 6)))]
dt.g$chr <- factor(dt.g$chr,
                   levels = lvls)

dt.g$reg <- factor(dt.g$reg)
levels(dt.g$reg)

dt.g$gene <- substr(dt.g$gene, 11, nchar(dt.g$gene) - 2)
dt.g$gene <- factor(dt.g$gene, levels = unique(dt.g$gene))
dt.g


gn <- subset(dt.g,
             chr == levels(chr)[i] &
               gene == "Xkr4")
gn

ggplot(gn) +
  geom_rect(aes(xmin = start,
                xmax = stop,
                ymin = rep(0, length(start)),
                ymax = rep(100, length(start)),
                fill = reg))

#*****************************************************
# Read-outs, alligned
dt1 <- fread("methylation seq/data/combined_ye_s1-4_03202017.csv")
dt1

setkey(dt1, start)
lvls <- unique(dt1$chr)[order(as.numeric(substr(unique(dt1$chr), 4, 6)))]
dt1$chr <- factor(dt1$chr,
                  levels = lvls)

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