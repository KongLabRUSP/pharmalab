setwd("C:\\Users\\dsargsy\\Desktop\\PhD\\Lab")
require(data.table)
require(reshape2)
require(knitr)
require(ricf)

dilut.GEN.nM <- rev(c(format(200000/2^(0:8),
                         scientific = TRUE, 
                         digits = 3),
                  "DMSO Only"))
dilut.DEI.nM <- rev(c(format(75000/2^(0:8),
                             scientific = TRUE, 
                             digits = 3),
                      "DMSO Only"))

# Read data in
# Files
files <- dir("MTS3_05162016_05192016\\data")
files

plates <- list()
for (i in 1:2) {
  plates[[i]] <- fread(paste("MTS3_05162016_05192016\\data",
                             files[i], 
                             sep = "\\"),
                       sep = ",",
                       header = TRUE,
                       skip = 26,
                       nrows = 8)
}
names(plates) <- files[1:2]
plates <- do.call("rbind", plates)

plates <- data.table(Hour = rep(c(24, 72), each = 8),
                     Plate = rep(rep(1, each = 8), 2),
                     Reading = rep(rep(1, each = 8), 2),
                     Repl = rep(c(0, rep(1:3, each = 2), 0), 2),
                     plates)
names(plates)[5:ncol(plates)] <- c("Row", paste("C", 1:12, sep = ""))
plates

maps <- fread("MTS3_05162016_05192016\\data\\TC1_MTS_Exp3_Map_05162016.csv",
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
                      "Genistein",
                      "Daidzein",
                      "Genistein",
                      "Daidzein",
                      "Daidzein",
                      "Genistein",
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

###########################################################
##Normalized Data
# Percent DMSO = 100*(Treatment - Median(Plate Background))/
# (Row DMSO - Median(Plate Background))

trt[, DMSO := no.bgrnd[Concentration == 9],
    by = list(Plate, Row)]
trt[,Normalized := 100*no.bgrnd/DMSO]

trt <- subset(trt, Concentration != 9)

dd <- dcast(trt, Plate + Repl + Concentration ~ Hour + Treatment,
            value.var = "Normalized")
head(dd)

# SAVE DD
saved.dd <- copy(dd)

#ANALYSIS I: Daidzein
dd.dei <- subset(saved.dd, select = c(1:3, 4, 6))
dd.dei$dilut <- 200000/2^as.numeric(as.character(dd.dei$Concentration))

m1 <- f.robicf(data.frame(Log10Conc = log10(dd.dei$dilut),
                          dd.dei[, 4:5]),
               ecf = c(60, 70))
s1 <- summary(m1, export = TRUE)
s1

ndx <- !do.call("c", lapply(s1, is.null))

t1 <- data.frame(Treatment = rep(names(s1)[ndx], each = 8),
                 do.call("rbind", s1[ndx]))
t1
write.csv(t1, "MTS3_05162016_05192016\\estimates.dei.csv")
nrow(do.call("rbind", s1[ndx]))
lines(m1,
      smooth = TRUE,
      xlab = "Log10 of Deidzein Concentration (nM)",
      ylab = "%DMSO",
      main = "All Plates, Fitted Curves (24 and 72 Hour)")

ciplot(m1,
       xlab = "Log10 of Concentration (nM)")

par(mfrow=c(1, 2))
plot.robicf(m1,
            xlab = "Deidzein Concentration (nM)",
            ylab = "% DMSO",
            smooth  = TRUE,
            xmarks = dilut.DEI.nM[-1],
            main = paste("Deidzein",
                         c(24, 72), 
                         "Hours"))
graphics.off()


#ANALYSIS II: Genistein
dd.gen <- subset(saved.dd, select = c(1:3, 5, 7))
dd.gen$dilut <- 75000/2^as.numeric(as.character(dd.gen$Concentration))

m2 <- f.robicf(data.frame(Log10Conc = log10(dd.gen$dilut),
                          dd.gen[, 4:5]),
               ecf = c(60, 70))
s2 <- summary(m2, export = TRUE)
s2

ndx <- !do.call("c", lapply(s2, is.null))

t2 <- data.frame(Treatment = rep(names(s2)[ndx], each = 8),
                 do.call("rbind", s2[ndx]))
t2
write.csv(t2, "MTS3_05162016_05192016\\estimates.gen.csv")
nrow(do.call("rbind", s2[ndx]))
lines(m2,
      smooth = TRUE,
      xlab = "Log10 of Genistein Concentration (nM)",
      ylab = "%DMSO",
      main = "All Plates, Fitted Curves (24 and 72 Hour)")

ciplot(m2,
       xlab = "Log10 of Concentration (nM)")

par(mfrow=c(1, 2))
plot.robicf(m2,
            xlab = "Genistein Concentration (nM)",
            ylab = "% DMSO",
            smooth  = TRUE,
            xmarks = dilut.DEI.nM[-1],
            main = paste("Deidzein",
                         c(24, 72), 
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
       aes(x = as.numeric(rev(Concentration)),
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
  scale_x_continuous("Concentration (nM)",
                   breaks = c(0:8 - 0.2, 0:8 + 0.2),
                   labels = c(dilut.DEI.nM[-1], dilut.GEN.nM[-1])) +
  scale_y_continuous("%DMSO") + 
  ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
  guides(fill = guide_legend(title = "Treatment",
                               title.position = "top",
                               nrow = 1))

ggplot(trt,
       aes(x = as.numeric(rev(Concentration)),
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
  scale_x_continuous("Concentration (nM)",
                     breaks = c(0:8 - 0.2, 0:8 + 0.2),
                     labels = c(dilut.DEI.nM[-1], dilut.GEN.nM[-1])) +
  scale_y_continuous("%DMSO") + 
  ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
  guides(fill = guide_legend(title = "Treatment",
                             title.position = "top",
                             nrow = 1))
####################################
# NOTES:
# 1. Did not redose at 48 Hour.