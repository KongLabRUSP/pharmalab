# Project: Epigenome of radiation induced mouse leukemia
# Protocol: [PROTO201800004](https://eiacuc.rutgers.edu/eIACUC/sd/Rooms/DisplayPages/LayoutInitial?Container=com.webridge.entity.Entity[OID[FCBC8146BFF7334A970ACCD252866C8F]])
# Script: Sample size calculation
# Investigator: Wenji Li
# Author: Davit Sargsyan
# Date: 03/05/2018

# * Describe step-by-step the procedure for administering the substance:
# Mice will be placed in a steel chamber and placed within the small animal irradiator
# Grp1,2,3 (5 mice) will receive 100 cGy whole body radiation harvested on day 1, 7 
# and 10 d respectively total 15 mice
# Grp4,5 and 6 (5 mice) will receive 200 cGy whole body radiation and harvested 
# similarly total 15 mice
# Grp7,8 and 9 (5 mice ) will receive 300 cGy whole body radiation and harvested 
# similarly total 15 mice total 45 mice
# each group will have a control non irradiated group total 45 mice
# each of the above groups will receive tha natural compund and irradiation or 
# no irradiation for another 90 mice grand total 180 mice plus 20 extras for 200 mice

require(data.table)
require(pwr)
require(bit64)
require(ggplot2)

# Power Curves----
std <- seq(10, 60, by = 10)
std
delta <- matrix(rep(seq(10, 90, by = 10),
                    length(std)),
                nrow = length(std),
                byrow = TRUE)
delta
h <- delta/std
h
colnames(h) <- rep(seq(10, 
                       90, 
                       by = 10))
rownames(h) <- paste("SD=", 
                     std,
                     "%",
                     sep = "")
h
n <- seq(3, 5, by = 1)
n

# Aplha = 0.05
res <- list()
for(i in 1:nrow(h)) {
  out <- list()
  for(j in 1:ncol(h)) {
    out[[j]] <- pwr.2p.test(h = h[i, j],
                            n = n,
                            sig.level = 0.05,
                            alternative = "g")$power
  }
  res[[i]] <- data.table(mu = rep(colnames(h), 
                                  each = length(n)),
                         n = rep(n,
                                 ncol(h)), 
                         power = do.call("c", out))
}
dt1 <- data.table(std = rep(rownames(h),
                            each = nrow(res[[1]])),
                  do.call("rbind",
                          res))
dt1$n <- factor(dt1$n)

p1 <- ggplot(dt1) +
  facet_wrap(~ std,
             nrow = 2) +
  geom_line(aes(x = mu,
                y = power,
                group = n,
                colour = n),
            size = 1) + 
  geom_hline(yintercept = 0.8,
             linetype = "dashed") +
  scale_x_discrete("Mean Differences (%)") + 
  scale_y_continuous("Power",
                     breaks = seq(0, 
                                  1, 
                                  by = 0.1)) + 
  ggtitle("") +
  guides(colour = guide_legend(title = "Number of Animals per Group")) +
  theme(plot.title = element_text(hjust = 0.5))
p1
tiff(filename = "figure1_power_wenji.tiff",
     height = 6,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

t1 <- dcast.data.table(dt1,
                       n + mu ~ std,
                       value.var = "power")
t1
write.csv(t1,
          file = "table1_power_wenji.csv",
          row.names = FALSE)