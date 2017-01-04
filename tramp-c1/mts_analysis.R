setwd("C:\\Users\\dsargsy\\Desktop\\PhD\\Lab")
require(data.table)
require(reshape2)
require(ricf)

d1 <- fread("DS_TC1_MTS_Plate1_24h_04202016.csv", sep = ",")
plot(d1$Value ~ d1$Treatment)

# MTS Controls
mts <- subset(d1, Treatment == "MTS")
plot(mts$Value,
     pch = 16,
     col = mts$Reading,
     main = "MTS Control",
     ylab = "Value")
abline(h = median(mts$Value), 
       lty = 2)
text(x = 50,
     y = median(mts$Value), 
     labels = paste("Median =",
                    round(median(mts$Value),
                          3)))
legend("topleft",
       legend = c("Reading1", "Reading2"),
       col = 1:2,
       pch = 16)

# Medium Controls
ctrl <- subset(d1, Treatment == "None")
ctrl
plot(ctrl$Value,
     pch = 16,
     col = ctrl$Reading,
     main = "Medium Control",
     ylab = "Value")
abline(h = median(ctrl$Value), 
       lty = 2)
text(x = 12,
     y = median(ctrl$Value), 
     labels = paste("Median =",
                    round(median(ctrl$Value),
                          3)))
legend("topleft",
       legend = c("Reading1", "Reading2"),
       col = 1:2,
       pch = 16)

# Treatment
trt <- subset(d1, !(Treatment %in% c("MTS", "None")))
trt$Treatment[trt$Treatment == "DMSO"] <- 9
plot(trt$Value ~ trt$Treatment,
     col = trt$Reading)
legend("topleft",
       legend = c("Reading1", "Reading2"),
       col = 1:2,
       pch = 1)

# Analysis
trt$dilut <- 100000/2^as.numeric(as.character(trt$Treatment))
m1 <- f.robicf(data.frame(conc = log10(trt$dilut),
                    Value = trt$Value),
               ecf = c(20, 80))
summary(m1)
plot(m1)

dilut.nM <- rev(c(format(100000/2^(0:8),
                         scientific = TRUE, 
                         digits = 3),
                  "DMSO Only"))

plot(m1,
     xlab = "Curcumin Concentration (nM)",
     smooth  = TRUE,
     col.point = trt$Reading,
     xmarks = dilut.nM)
abline(h = median(mts$Value), lty = 2, col = "blue")

# Separate readings
d2 <- fread("DS_TC1_MTS_Plate1_Readings_24h_04202016.csv", sep = ",")
trt2 <- subset(d2, !(Treatment %in% c("MTS", "None")))
trt2$Treatment[trt2$Treatment == "DMSO"] <- 9

trt2$dilut <- 100000/2^as.numeric(as.character(trt2$Treatment))
m2 <- f.robicf(data.frame(conc = log10(trt2$dilut),
                          H24R1 = trt2$H24R1,
                          H24R2 = trt2$H24R2),
               ecf = c(20, 80))
summary(m2)
plot(m2,
     xlab = "Curcumin Concentration (nM)",
     smooth  = TRUE,
     xmarks = dilut.nM)
lines(m2,
      xlab = "Curcumin Concentration (nM)",
      smooth  = TRUE,
      xmarks = dilut.nM)

##############################################
##############################################
# ALL PLATES
plates <- fread("TC1_Curcumin_MTS_Plates_Exp1.csv", sep = ",", header = TRUE)
maps <- fread("TC1_Curcumin_MTS_Map_Exp1.csv", sep = ",", header = TRUE)


tmp <- subset(plates, Plate == 1 & Reading == 1)

d1 <- reshape(data = plates, 
              varying = list(names(plates)[5:ncol(plates)]),
              direction = "long",
              v.names = "Value")

d1$Treatment <- reshape(data = maps, 
                        varying = list(names(maps)[4:ncol(maps)]),
                        direction = "long",
                        v.names = "Treatment")$Treatment
d1$time <- NULL

# Separate treated cells
trt <- subset(d1, !(Treatment %in% c("None", "MTS")))
plot(trt$Value ~ trt$Treatment, 
     col = as.numeric(factor(d1$Hour)),
     xlab = "Concentration",
     ylab = "Value")

dd <- dcast(trt, Row + Treatment ~ Reading + Hour, value.var = "Value")
dd$Treatment[dd$Treatment == "DMSO"] <- 9

# Analysis
dd$dilut <- 100000/2^as.numeric(as.character(dd$Treatment))
m1 <- f.robicf(data.frame(conc = log10(dd$dilut),
                          dd[, c(3, 6, 4, 7, 5, 8)]),
               ecf = c(20, 80))
summary(m1)

dilut.nM <- rev(c(format(100000/2^(0:8),
                         scientific = TRUE, 
                         digits = 3),
                  "DMSO Only"))

plot(m1,
     xlab = "Curcumin Concentration (nM)",
     smooth  = TRUE,
     xmarks = dilut.nM)

lines(m1,
      smooth = TRUE,
      xlab = "Log of Curcumin Concentration")