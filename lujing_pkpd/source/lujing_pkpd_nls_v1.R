# |----------------------------------------------------------------------------------|
# | Project: Curcumin PK/PD                                                          |
# | Script: Analysis of PK/PD data                                                   |
# | Scientist: Yuquing (Anne) Yang, Lujing Wang                                      |
# | Author: Davit Sargsyan                                                           |
# | Created: 02/01/2018                                                              |
# | Modified: 02/03/2018: replaced manual calculations with nls                      |
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

# Estimate PK parameters with 'nls'----
dt.exp$logConc <- log(dt.exp$`Curcumin(ng/mL)`)
plot(dt.exp$logConc ~ dt.exp$Hours)
abline(h = 1, 
       lty = 2)

# Remove very small values: leave <= 3 hours data 
dt.exp <- subset(dt.exp,
              Hours <= 3 & 
                !is.na(logConc),
              select = c(1, 5))
plot(dt.exp)

# Initial conditions
init <- list(A = 3,
             ka = 0.2,
             B = 4,
             kb = 6)
m1 <- nls(logConc ~ A*exp(-ka*Hours) + B*exp(-kb*Hours),
          data = dt.exp,
          start = init,
          control = nls.control(maxiter = 1000, 
                                minFactor = 1/4096),
          na.action = na.omit)
summary(m1)

# Predicted data----
prd <- data.table(Predicted = predict(m1,
                                      newdata = list(Hours = seq(0, 3, 0.01))),
                  Hours = seq(0, 3, 0.01))
prd

p3 <- ggplot() +
  geom_hline(yintercept = log(c(10, 50)),
             linetype = "dotted") +
  geom_line(data = prd,
            aes(x = Hours,
                y = Predicted),
            linetype = "dashed",
            color = "blue",
            size = 1.5) +
  geom_point(data = dt.exp,
             aes(x = Hours,
                 y = logConc),
             size = 3,
             shape = 21,
             fill = "red",
             alpha = 0.5) +
  scale_x_continuous(limits = c(0, 3),
                     breaks = unique(dt.exp$Hours)) +
  scale_y_continuous("Log(Concentration(ng/mL))") +
  ggtitle("Curcumin: Observed Concentration in Blood (Red Dots)\nand Concentration Predicted by NLS (Blue Dahsed Line)") +
  theme(plot.title = element_text(hjust = 0.5))
p3

tiff(filename = "tmp/cur_log_y_nls_pred.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p3)
graphics.off()

# AUC: cumulative exposure over time----
prd <- data.table(Predicted = exp(predict(m1,
                                      newdata = list(Hours = unique(dt.exp$Hours)))),
                  Hours = unique(dt.exp$Hours))
prd$dTime <- prd$Hours - c(0, prd$Hours[-7])
prd$muConc <- (prd$Predicted + c(0, prd$Predicted[-7]))/2
prd$auc <- prd$muConc*prd$dTime
prd$cumConc <- cumsum(prd$auc)
prd

p4 <- ggplot(data = prd,
             aes(x = Hours,
                 y = cumConc)) +
  geom_line() +
  geom_point(size = 3,
             shape = 21,
             fill = "red",
             alpha = 0.5) +
  scale_x_continuous(limits = c(0, 3),
                     breaks = unique(dt.exp$Hours)) +
  scale_y_continuous("Total Exposure (ng/mL)",
                     breaks = seq(0, 110, 10)) +
  ggtitle("Curcumin: Prediceted Cumulative Systemic Exposure") +
  theme(plot.title = element_text(hjust = 0.5))
p4

tiff(filename = "tmp/cur_total_exposure.tiff",
     height = 5,
     width = 6,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p4)
graphics.off()

# Part II: Gene expressions----
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

# Part III: predict iNOS over time with nls----
init <- list(E0 = 0.5, # baseline effect
             Emax = 8, # maximum effect (set to 10 for LPS-only)
             ET50 = 1, # time of 50% response
             IT50 = 1.5, # time of 50% inhibition
             kin = 3, # rate of response
             kout = 2) # rate of inhibition

dt2 <- dt1[dt1$Treatment == "CUR+LPS", c(2, 4)]
# dt2 <- dt1[dt1$Treatment == "LPS", c(2, 4)]
plot(dt2)

dt2 <- log(dt2 + 1)
dt2
plot(dt2)

dt2$pred <- init$E0 + (init$Emax - init$E0)/
  (1 + exp((init$ET50 - dt2$Hours)*init$kin)) - 
  (init$Emax - init$E0)/
  (1 + exp((init$IT50 - dt2$Hours)*init$kout))
dt2
lines(dt2$pred ~ dt2$Hours)

m1 <- nls(iNOS ~ E0 + (Emax - E0)/
            (1 + exp((ET50 - Hours)*kin)) - 
            (Emax - E0)/
            (1 + exp((IT50 - Hours)*kout)),
          data = dt2,
          start = init,
          control = nls.control(maxiter = 1000, 
                                minFactor = 1/4096),
          na.action = na.omit)
s1 <- summary(m1)
s1
exp(s1$coefficients[, 1])


# Predicted data----
prd <- data.table(Predicted = predict(m1,
                                      newdata = list(Hours = log(seq(1, 24, 0.01)))),
                  Hours = log(seq(1, 24, 0.01)))
prd

p7 <- ggplot() +
  geom_line(data = prd,
            aes(x = Hours,
                y = Predicted),
            linetype = "dashed",
            color = "blue",
            size = 1.5) +
  geom_point(data = dt2,
             aes(x = Hours,
                 y = iNOS),
             size = 3,
             shape = 21,
             fill = "red",
             alpha = 0.5) +
  scale_y_continuous("Log(iNOS)") +
  scale_x_continuous("Log(Hours)") +
  ggtitle("Curcumin: Observed (Red Dots) vs. Predicted (Blue Line) iNOS") +
  theme(plot.title = element_text(hjust = 0.5))
p7

tiff(filename = "tmp/cur_inos_pred.tiff",
     height = 5,
     width = 7,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p7)
graphics.off()

# sink()