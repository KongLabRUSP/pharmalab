# |----------------------------------------------------------------------------------|
# | Project: Curcumin PK/PD                                                          |
# | Script: Analysis of PK/PD data                                                   |
# | Scientist: Yuquing (Anne) Yang, Lujing Wang                                      |
# | Author: Davit Sargsyan                                                           |
# | Created: 02/01/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_lujing_pkpd_v1.txt")

# Header----
require(data.table)
require(ggplot2)

setwd("lujing_pkpd")
# Source: '~/docs/PD data.xlsx' and '~/docs/Wenji's new method PK data.xlsx'
dt1 <- fread("data/lujing_pkpd_data_comb.csv")
dt1$Time <- factor(dt1$Time,
                   levels = unique(dt1$Time))
dt1

# PartI: exposure----
dt2 <- subset(dt1,
              Treatment == "CUR+LPS")
dt2

dt.exp <- melt.data.table(data = dt2,
                          id.vars = 2,
                          measure.vars = 6:8,
                          value.name = "Curcumin(ng/mL)")
dt.exp

dt.exp$hh <- dt.exp$Hours + 10^c(rep(0.01, 12),
                              rep(0.02, 12),
                              rep(0.03, 12))

p1 <- ggplot(dt.exp,
             aes(x = Hours,
                 y = `Curcumin(ng/mL)`)) +
  geom_point(size = 3,
             shape = 21,
             fill = "red",
             alpha = 0.5) +
  geom_hline(yintercept = c(10, 50),
             linetype = "dashed") +
  scale_y_continuous(breaks = c(10, 50, 100, 250, 500, 1000),
                     limits = c(0, 1000)) +
  scale_x_continuous(breaks = unique(dt.exp$Hours)[-c(2:4)]) +
  ggtitle("Curcumin Concentration in Blood, Original Scale") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0,
                                   hjust = 0.5))
p1

tiff(filename = "tmp/cur_orig.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

p2 <- ggplot(dt.exp,
             aes(x = hh,
                 y = `Curcumin(ng/mL)`)) +
  geom_point(size = 3,
             shape = 21,
             fill = "red",
             alpha = 0.5) +
  geom_hline(yintercept = c(10, 50),
             linetype = "dashed") +
  scale_y_log10(limits = c(1, 1000),
                breaks = c(1, 10, 50, 100, 1000)) +
  scale_x_log10("Hours",
                breaks = dt.exp$hh[13:24],
                labels = unique(dt.exp$Hours)) +
  ggtitle("Curcumin Concentration in Blood, Semi-Log Axes") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p2
tiff(filename = "tmp/cur_semi-log_semi-log.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Predict elimination rate (20 min to 1 hour)----
tmp <- subset(dt.exp,
              Hours <= 1)
tmp$logconc <- log(tmp$`Curcumin(ng/mL)`)
tmp

# Predict elimination
m1 <- lm(logconc ~ Hours,
         data = tmp[tmp$Hours >= 0.33])
summary(m1)
tmp$elim.prd <- predict(m1,
                        newdata = list(Hours = tmp$Hours))

# Predict distribution
tmp$dist <- tmp$logconc - tmp$elim.prd
m2 <- lm(dist ~ Hours,
         data = tmp[tmp$Hours < 0.33])
summary(m2)
tmp$dist.prd <- predict(m2,
                        newdata = list(Hours = tmp$Hours))
tmp$dist[tmp$dist.prd < 0] <- NA
tmp$dist.prd[tmp$dist.prd < 0] <- NA

p3 <- ggplot(tmp,
       aes(x = Hours,
           y = logconc)) +
  geom_point(size = 3,
             shape = 21,
             fill = "red",
             alpha = 0.5) +
  geom_hline(yintercept = log10(c(10, 50)),
             linetype = "dotted") +
  geom_line(aes(x = Hours,
                y = elim.prd),
            linetype = "dashed",
            color = "blue") +
  geom_point(aes(x = Hours,
                 y = dist),
             size = 3,
             shape = 21,
             fill = "blue") +
  geom_line(aes(x = Hours,
                y = dist.prd),
            linetype = "dashed",
            color = "red") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = unique(dt.exp$Hours)) +
  scale_y_continuous("Log(Concentration(ng/mL))") +
  ggtitle("Curcumin Concentration in Blood, Log Concentration") +
  theme(plot.title = element_text(hjust = 0.5))
p3

tiff(filename = "tmp/cur_log_y_dist_elim.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()

# Predict total concentration----
foo <- function(m1, m2){
  tm <- seq(0, 1, 0.01)
  # point of connection
  x <- -m2$coefficients[1]/m2$coefficients[2]
  y <- predict(m1,
               newdata = list(Hours = x))
  prd <- predict(m1,
                 newdata = list(Hours = tm))
  prd[tm < x] <- prd[tm < x] + predict(m2,
                                       newdata = list(Hours = tm[tm < x]))
  return(data.table(Hours = tm,
                    Predicted = prd))
}

prd <- foo(m1, m2)
p4 <- ggplot(data = tmp,
             aes(x = Hours,
                 y = logconc)) +
  geom_point(size = 3,
             shape = 21,
             fill = "red",
             alpha = 0.5) +
  geom_line(data = prd,
            aes(x = Hours,
                y = Predicted),
            linetype = "dashed",
            color = "blue") +
  scale_x_continuous(limits = c(0, 1),
                     breaks = unique(dt.exp$Hours)) +
  scale_y_continuous("Log(Concentration(ng/mL))") +
  ggtitle("Curcumin Concentration in Blood, Log Concentration") +
  theme(plot.title = element_text(hjust = 0.5))
p4
tiff(filename = "tmp/cur_log_y_dist_elim_pred.tiff",
     height = 5,
     width = 5,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4)
graphics.off()

# PK parameters----
# a. C(0)
mean(tmp$`Curcumin(ng/mL)`[tmp$Hours == 0])
# almost same as predicted
exp(m1$coefficients[1] + m2$coefficients[1])

sd(tmp$`Curcumin(ng/mL)`[tmp$Hours == 0])

# b. k1 (distribution phase) (1/min)
-m2$coefficients[2]

# c. k2 (elimination phase) (1/min)
-m1$coefficients[2]

# Half-life (elimination) (min)
-60*log(2)/m1$coefficients[2]

# Gene expressions----
tmp <- melt.data.table(data = dt1,
                       id.vars = 1:3,
                       measure.vars = 4:5,
                       variable.name = "Gene",
                       value.name = "Expression")
tmp <- tmp[!is.na(tmp$Expression), ]

p5 <- ggplot(tmp,
       aes(x = Time,
           y = Expression,
           group = Treatment,
           fill = Treatment)) +
  facet_wrap(~ Gene,
             nrow = 1,
             scales = "free_y") +
  geom_line(position = position_dodge(0.3)) +
  geom_point(size = 3,
             shape = 21,
             alpha = 0.8,
             position = position_dodge(0.3)) +
  ggtitle("Gene Expressions After Exposure to Curcumin") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p5
tiff(filename = "tmp/cur_gene_expr.tiff",
     height = 5,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p5)
graphics.off()

# Expression inhibition (Cur+LPS vs LPS)----
dt3 <- subset(dt1,
              dt1$Treatment == "LPS")
dt3
dt2$d_iNOS <- dt3$iNOS/dt2$iNOS
dt2$d_TNFa <- dt3$TNFa/dt2$TNFa
dt2

ddt <- melt.data.table(data = dt2,
                       id.vars = 1:2,
                       measure.vars = 9:10,
                       value.name = "LPS/Cur+LPS",
                       variable.name = "Gene")
ddt <- ddt[!is.na(ddt$`LPS/Cur+LPS`), ]
ddt

p6 <- ggplot(ddt,
             aes(x = Time,
                 y = `LPS/Cur+LPS`,
                 group = Gene,
                 fill = Gene)) +
  facet_wrap(~ Gene,
             # scales = "free_y",
             nrow = 1) +
  geom_line(position = position_dodge(0.3)) +
  geom_point(size = 3,
             shape = 21,
             alpha = 0.8,
             position = position_dodge(0.3)) +
  geom_hline(yintercept = 1,
             linetype = "dashed") + 
  stat_smooth() +
  ggtitle("Gene Expressions Ratios After Exposure to Curcumin\nLOESS Smoother") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p6

tiff(filename = "tmp/cur_gene_expr_diff.tiff",
     height = 5,
     width = 12,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p6)
graphics.off()

# sink()