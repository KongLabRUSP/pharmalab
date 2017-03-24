require(data.table)
require(ggplot2)

dt1 <- fread("methylation seq/results.csv")
dt1
setkey(dt1, start)
dt1$chr <- factor(dt1$chr)


for (i in 1:nlevels(dt1$chr)) {
  tmp <- droplevels(subset(dt1, chr == levels(chr)[i]))
  tmp
  min(tmp$start)
  max(tmp$end)
  tmp$ystart <- rep(0, nrow(tmp))
  tmp$yend <- rep(1, nrow(tmp))
  tmp$red.meth <- "green"
  tmp$red.meth[tmp$`Control->Exptl:diff` < 0] <- "red"
  tmp$red.meth[tmp$`Control->Exptl:pval` > 0.01] <- "black"
  tmp
  
  sum(tmp$`Control->Exptl:pval` <= 0.01)/nrow((tmp))
  
  # ggplot(tmp) +
  #   geom_rect(aes(xmin = start/10^6,
  #                 xmax = end/10^6,
  #                 ymin = ystart,
  #                 ymax = CpG))
  
  tiff(filename = paste("methylation seq/plots/plot", i, ".tiff", sep = ""),
       height = 5,
       width = 10,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  
  p1 <- ggplot(tmp) +
    geom_segment(aes(x = start,
                     xend = start,
                     y = ystart,
                     yend = CpG),
                 col = tmp$red.meth) +
    scale_x_continuous("Starting Location of Regions") +
    scale_y_continuous("Total CpG") +
    ggtitle(paste("Chromosome:", tmp$chr[1]))
  print(p1)
  
  graphics.off()
}











dt2 <- data.table(xstart = c(0, 3, 10),
                  xend = c(1, 5, 20),
                  ystart = rep(0, 3),
                  yend = rep(1, 3))
ggplot(dt2) +
  geom_rect(aes(xmin = xstart,
            xmax = xend,
            ymin = ystart,
            ymax = yend),
            fill = 1:3)


df <- data.frame(
  x1 = c(0, 3, 10),
  y = sample(10, 20, replace = TRUE)
)
ggplot(df, aes(xmin = x, xmax = x + 1, ymin = y, ymax = y + 2)) +
  geom_rect()

