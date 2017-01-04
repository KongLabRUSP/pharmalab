setwd("C:\\Users\\dsargsy\\Desktop\\PhD\\Lab")
# http://stackoverflow.com/questions/27661325/unable-to-load-rjava-on-r
# Sys.setenv(JAVA_HOME="C:\\Program Files (x86)\\Java\\jre1.8.0_73")
# Sys.getenv("JAVA_HOME")
# Sys.getenv("JAVA_HOME")
# require(xlsx)
require(data.table)
require(reshape2)
require(knitr)
require(ricf)

dilut.nM <- rev(c(format(100000/2^(0:8),
                         scientific = TRUE, 
                         digits = 3),
                  "DMSO Only"))

# Read data in
# Files
files <- dir("MTS2_04292016_05022016\\data")
files

plates <- list()
for (i in 1:12) {
  plates[[i]] <- fread(paste("MTS2_04292016_05022016\\data",
                             files[i], 
                             sep = "\\"),
                       sep = ",",
                       header = TRUE,
                       skip = 26,
                       nrows = 8)
}
names(plates) <- files[1:12]
plates <- do.call("rbind", plates)

plates <- data.table(Hour = rep(c(24, 48, 72), each = 32),
                     Plate = rep(rep(1:2, each = 16), 3),
                     Reading = rep(rep(1:2, each = 8), 6),
                     Repl = rep(c(0, rep(1:2, each = 3), 0), 12),
                     plates)
names(plates)[5:ncol(plates)] <- c("Row", paste("C", 1:12, sep = ""))
plates

maps <- fread("MTS2_04292016_05022016\\data\\TC1_MTS_Exp2_Map_05022016.csv",
              sep = ",", 
              header = TRUE)

d1 <- reshape(data = plates, 
              varying = list(names(plates)[6:ncol(plates)]),
              direction = "long",
              v.names = "Value")

d1$Concentration <- reshape(data = maps, 
                        varying = list(names(maps)[4:ncol(maps)]),
                        direction = "long",
                        v.names = "Concentration")$Concentration

d1$Treatment <- rep(c("Medium",
                      "Curcumin",
                      "Isoliquiritigenin",
                      "Liquiritigenin",
                      "Curcumin",
                      "Isoliquiritigenin",
                      "Liquiritigenin",
                      "Medium"),
                    12)
d1$time <- NULL
d1

# Median of background, i.e. medium only, per plate
d1[, bgrnd := median(Value[Concentration == "None"]),
   by = list(Hour, Plate, Reading)]

# Subtract background
d1[, no.bgrnd := Value - bgrnd]
summary(d1)

# Separate treated cells
trt <- subset(d1, Concentration != "None")
trt$Concentration[trt$Concentration == "DMSO"] <- 9

# Remove Reading2
trt <- subset(trt, Reading == 1)
trt$Reading <- NULL

trt <- droplevels(trt)

plot(trt$Value ~ jitter(rev(as.numeric(trt$Concentration))), 
     col = as.numeric(factor(trt$Hour)),
     pch = as.numeric(factor(trt$Treatment)),
     xlab = "Concentration (nM)",
     ylab = "Read-out",
     main = "All Plates, Raw Data (24, 48 and 72 Hour)",
     xaxt = "n")
axis(side = 1,
     at = 0:9,
     labels = FALSE,
     las = 2)
text(x = 0:9,
     y = 0.2,
     labels = dilut.nM,
     srt = 30,
     pos = 2,
     xpd = 2)
legend("topright", 
       legend = paste(unique(trt$Hour), "Hour"),
       col = unique(as.numeric(factor(trt$Hour))),
       pch = 16)
legend("topleft", 
       legend = unique(trt$Treatment),
       pch = unique(as.numeric(factor(trt$Treatment))))

# dd <- dcast(trt, Plate + Repl + Concentration ~ Hour + Treatment + Reading,
#             value.var = "Value")
dd <- dcast(trt, Plate + Repl + Concentration ~ Hour + Treatment,
            value.var = "Value")
head(dd)
# trt.names <- c("Hour24/Reading1",
#                "Hour48/Reading1",
#                "Hour72/Reading1",
#                "Hour24/Reading2",
#                "Hour48/Reading2",
#                "Hour72/Reading2")
# names(dd)[5:ncol(dd)] <- trt.names

# Analysis
dd$dilut <- 100000/2^as.numeric(as.character(dd$Concentration))
# m1 <- f.robicf(data.frame(conc = log10(dd$dilut)[dd$Drug == "Curcumin"],
#                           dd[dd$Drug == "Curcumin", c(5, 8, 6, 9, 7, 10)]))
m1 <- f.robicf(data.frame(Log10Conc = log10(dd$dilut),
                          dd[, 4:12]))
s1 <- summary(m1, export = TRUE)
s1
lines(m1,
      smooth = TRUE,
      xlab = "Log of Curcumin Concentration (nM)",
      ylab = "Read-out",
      main = "All Plates, Fitted Curves (24, 48 and 72 Hour)")
plot(m1, smooth = TRUE)

###########################################################
##Normalized Data
# Percent DMSO = 100*(Treatment - Median(Plate Background))/
# (Row DMSO - Median(Plate Background))

# trt[, DMSO := no.bgrnd[Concentration == 9],
#     by = list(Plate, Reading, Row)]
trt[, DMSO := no.bgrnd[Concentration == 9],
    by = list(Plate, Row)]
trt[,Normalized := 100*no.bgrnd/DMSO]

trt <- subset(trt, Concentration != 9)

plot(trt$Normalized ~ jitter(rev(as.numeric(trt$Concentration))), 
     col = as.numeric(factor(trt$Hour)),
     pch = as.numeric(factor(trt$Treatment)),
     xlab = "Concentration (nM)",
     ylab = "Read-out",
     main = "All Plates, Raw Data (24, 48 and 72 Hour)",
     xaxt = "n")
axis(side = 1,
     at = 0:9,
     labels = FALSE,
     las = 2)
text(x = 0:9,
     y = 0.2,
     labels = dilut.nM,
     srt = 30,
     pos = 2,
     xpd = 2)
legend("topright", 
       legend = paste(unique(trt$Hour), "Hour"),
       col = unique(as.numeric(factor(trt$Hour))),
       pch = 16)
legend("topleft", 
       legend = unique(trt$Treatment),
       pch = unique(as.numeric(factor(trt$Treatment))))

# dd <- dcast(trt, Plate + Repl + Concentration ~ Hour + Treatment + Reading,
#             value.var = "Normalized")
dd <- dcast(trt, Plate + Repl + Concentration ~ Hour + Treatment,
            value.var = "Normalized")
head(dd)

# Analysis
dd$dilut <- 100000/2^as.numeric(as.character(dd$Concentration))


m1 <- f.robicf(data.frame(Log10Conc = log10(dd$dilut),
                          dd[, 4:12]),
               ecf = c(60, 70))
s1 <- summary(m1, export = TRUE)
s1

ndx <- !do.call("c", lapply(s1, is.null))

t1 <- data.frame(Treatment = rep(names(s1)[ndx], each = 8),
                 do.call("rbind", s1[ndx]))
t1
write.csv(t1, "MTS2_04292016_05022016\\estimates.csv")
nrow(do.call("rbind", s1[ndx]))
lines(m1,
      smooth = TRUE,
      xlab = "Log10 of Curcumin Concentration (nM)",
      ylab = "%DMSO",
      main = "All Plates, Fitted Curves (24, 48 and 72 Hour)")

ciplot(m1,
       xlab = "Log10 of Concentration (nM)")

par(mfrow=c(1, 3))
plot.robicf(m1,
            xlab = "Curcumin Concentration (nM)",
            ylab = "% DMSO",
            smooth  = TRUE,
            xmarks = dilut.nM[-1],
            col.point = rep(rep(c("blue", "red"),
                                each = 2), 
                            9),
            main = paste(unique(trt$Treatment),
                         "\n",
                         rep(c(24, 48, 72), each = 3), 
                         "Hours"))
graphics.off()

##########################################
trt[, Mean := mean(Normalized),
    by = list(Hour, Concentration, Treatment)]
trt[, SD := sd(Normalized),
    by = list(Hour, Concentration, Treatment)]
trt[, Time := paste(Hour, "Hours")]

# Source:
# http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/
ggplot(trt, 
       aes(x = rev(Concentration),
           y = Mean,
           fill = Treatment)) + 
  facet_wrap(~ Time,
             ncol = 1) +
  geom_bar(position = position_dodge(),
           stat="identity") +
  geom_errorbar(aes(ymax = Mean + SD,
                    ymin = Mean - SD),
                width =.2,
                position = position_dodge(.9)) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  scale_x_discrete("Concentration (nM)",
                   breaks = 0:8,
                   labels = rev(unique(dd$dilut))) +
  scale_y_continuous("%DMSO") + 
  ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
  guides(fill = guide_legend(title = "Treatment",
                               title.position = "top",
                               nrow = 1))

ggplot(trt,
       aes(x = rev(Concentration),
           y = Mean,
           colour = Treatment,
           group = Treatment)) + 
  facet_wrap(~ Time,
             ncol = 1) +
  geom_line(position = position_dodge(0.3),
            size = 1) +
  geom_point(position = position_dodge(0.3),
             size = 3) +
  geom_errorbar(aes(ymax = Mean + SD,
                    ymin = Mean - SD),
                width =.4,
                size = 1,
                position = position_dodge(0.3)) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  scale_x_discrete("Concentration (nM)",
                   breaks = 0:8,
                   labels = rev(unique(dd$dilut))) +
  scale_y_continuous("%DMSO") + 
  ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
  guides(fill = guide_legend(title = "Treatment",
                             title.position = "top",
                             nrow = 1))

####################################
# Bring in Experiment1 data
files.e1 <- dir("MTS1_04182016_04222016\\data")
files.e1

plates.e1 <- list()
for (i in 1:6) {
  plates.e1[[i]] <- fread(paste("MTS1_04182016_04222016\\data",
                             files.e1[i], 
                             sep = "\\"),
                       sep = ",",
                       header = TRUE,
                       skip = 26,
                       nrows = 8)
}
names(plates.e1) <- files.e1[1:6]
plates.e1 <- do.call("rbind", plates.e1)

plates.e1 <- data.table(Hour = rep(c(24, 48, 72), each = 16),
                     Plate = rep(rep(1, each = 16), 3),
                     Reading = rep(rep(1:2, each = 8), 3),
                     Repl = rep(c(0, 1:6, 0), 6),
                     plates.e1)
names(plates.e1)[5:ncol(plates.e1)] <- c("Row", paste("C", 1:12, sep = ""))
plates.e1

maps.e1 <- fread("MTS1_04182016_04222016\\data\\TC1_MTS_Exp1_Map_04222016.csv",
              sep = ",", 
              header = TRUE)

d2 <- reshape(data = plates.e1, 
              varying = list(names(plates.e1)[6:ncol(plates.e1)]),
              direction = "long",
              v.names = "Value")

d2$Concentration <- reshape(data = maps.e1, 
                            varying = list(names(maps.e1)[4:ncol(maps.e1)]),
                            direction = "long",
                            v.names = "Concentration")$Concentration

d2$Treatment <- rep(c("Medium",
                      rep("Curcumin EX1",
                          6),
                      "Medium"),
                    6)
d2$time <- NULL
d2

# Median of background, i.e. medium only, per plate
d2[, bgrnd := median(Value[Concentration == "None"]),
   by = list(Hour, Plate, Reading)]

# Subtract background
d2[, no.bgrnd := Value - bgrnd]
summary(d2)

# Separate treated cells
trt.e1 <- subset(d2, !(Concentration %in% c("None", "MTS")))
trt.e1$Concentration[trt.e1$Concentration == "DMSO"] <- 9

# Remove Reading2
trt.e1 <- subset(trt.e1, Reading == 1)
trt.e1$Reading <- NULL

trt.e1[, DMSO := no.bgrnd[Concentration == 9],
    by = list(Plate, Row)]
trt.e1[,Normalized := 100*no.bgrnd/DMSO]

trt.e1 <- subset(trt.e1, Concentration != 9)

trt.e1 <- droplevels(trt.e1)

trt
trt.e1

comb <- rbindlist(list(trt.e1, trt[, -c(13:15), with = FALSE]))

comb.w <- dcast(comb, Plate + Repl + Concentration ~ Hour + Treatment,
            value.var = "Normalized")
head(comb.w)

# Analysis
comb.w$dilut <- 100000/2^as.numeric(as.character(comb.w$Concentration))


m3 <- f.robicf(data.frame(Log10Conc = log10(comb.w$dilut),
                          comb.w[, 4:12]))
s3 <- summary(m3, export = TRUE)
s3

ciplot.robicf(m3, 
              main = "95% C.I. for IC50 Estimates",
              )

# comb[, Mean := mean(Normalized),
#     by = list(Hour, Concentration, Treatment)]
# comb[, SD := sd(Normalized),
#     by = list(Hour, Concentration, Treatment)]
# comb[, Time := paste(Hour, "Hours")]

summ <- aggregate(comb$Normalized,
                  by = list(Hour = comb$Hour, 
                            Concentration = comb$Concentration,
                            Treatment = comb$Treatment),
                  FUN = mean)
summ$sd <- aggregate(comb$Normalized,
                     by = list(Hour = comb$Hour, 
                               Concentration = comb$Concentration,
                               Treatment = comb$Treatment),
                     FUN = sd)$x
summ$Time <- paste(summ$Hour, "Hours")
head(summ)
# Source:
# http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/
ggplot(summ, 
       aes(x = rev(Concentration),
           y = x,
           fill = Treatment)) + 
  facet_wrap(~ Time,
             ncol = 1) +
  geom_bar(position = position_dodge(),
           stat = "identity") +
  geom_errorbar(aes(ymax = x + sd,
                    ymin = x - sd),
                width =.2,
                position = position_dodge(.9)) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  scale_x_discrete("Concentration (nM)",
                   breaks = 0:8,
                   labels = rev(dilut.nM[-1])) +
  scale_y_continuous("%DMSO") + 
  ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
  guides(fill = guide_legend(title = "Treatment",
                             title.position = "top",
                             nrow = 1))

ggplot(summ,
       aes(x = rev(Concentration),
           y = x,
           colour = Treatment,
           group = Treatment)) + 
  facet_wrap(~ Time,
             ncol = 1) +
  geom_line(position = position_dodge(0.3),
            size = 1) +
  geom_point(position = position_dodge(0.3),
             size = 3) +
  geom_errorbar(aes(ymax = x + sd,
                    ymin = x - sd),
                width =.4,
                size = 1,
                position = position_dodge(0.3)) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  scale_x_discrete("Concentration (nM)",
                   breaks = 0:8,
                   labels = rev(dilut.nM[-1])) +
  scale_y_continuous("%DMSO") + 
  ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
  guides(fill = guide_legend(title = "Treatment",
                             title.position = "top",
                             nrow = 1))
####################################
# NOTES:
# 1. Did not redose at 48 Hour.
# 2. Hook effect was observed at low doses. A model that can describe the 
#    hook might be more appropriate, especially if the hook effect becomes more pronounced.
# 3. Smaller intervals between the doses can be beneficial as we can get 
#    smoother curves around IC50 values and decrease standard errors of the estimates.
# 4. Starting dose for ISO (100 microM) was too low. As a result, I missed bottom
#    part of the curve. Starting concentration needs to be at least X16 of that (X20 = 2miliM)
# 4. Starting dose for LIQ (100 microM) was too low. As a result, I only observed the upper
#    part of the curve. NO estimate of starting concentration can be obtained from the data.