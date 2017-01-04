setwd("C:\\svn_mirror\\PhD\\Lab")
require(data.table)
require(reshape2)
require(knitr)
require(ricf)

dilut.nM <- rev(c(format(90000/2^(0:8),
                         scientific = TRUE, 
                         digits = 3),
                  "DMSO Only"))

# Read data in
# Files
files <- dir("MTS_08012016_Echinatin\\data")
files

plates <- list()
for (i in 2:7) {
  plates[[i - 1]] <- fread(paste("MTS_08012016_Echinatin\\data",
                             files[i], 
                             sep = "\\"),
                       sep = ",",
                       header = TRUE,
                       skip = 26,
                       nrows = 8)
}
names(plates) <- files[2:7]
plates <- do.call("rbind", plates)

plates <- data.table(Hour = rep(rep(c(24, 72, 120), each = 8), 2), 
                     Plate = rep(1:2, each = 8*3),
                     Repl = rep(c(0, rep(1:3, each = 2), 0), 6),
                     plates)
names(plates)[4:ncol(plates)] <- c("Row", paste("C", 1:12, sep = ""))
plates

# Map
maps <- fread("MTS_08012016_Echinatin\\data\\DS_TC1_Echinatin_MTS_Map_Exp5n6.csv",
              sep = ",", 
              header = TRUE)

d1 <- reshape(data = plates, 
              varying = list(names(plates)[5:ncol(plates)]),
              direction = "long",
              v.names = "Value")

d1$Concentration <- do.call("c", lapply(maps[, -1, with = FALSE], function(a) {rep(a, 6)}))

d1$Treatment <- rep(c("Medium",
                      "Echinatin",
                      "Curcumin",
                      "Echinatin",
                      "Curcumin",
                      "Echinatin",
                      "Curcumin",
                      "Medium"),
                    6)
d1$time <- NULL
d1

# Median of background, i.e. medium only, per plate
d1[, bgrnd := median(Value[Concentration == "MTS"]),
   by = list(Hour, Plate)]

# Subtract background
d1[, no.bgrnd := Value - bgrnd]
summary(d1)

# Separate treated cells
trt <- subset(d1, Concentration != "MTS")
trt$Concentration[trt$Concentration == "DMSO"] <- 0
trt <- droplevels(trt)

###########################################################
##Normalized Data
# Percent DMSO = 100*(Treatment - Median(Plate Background))/
# (Row DMSO - Median(Plate Background))

trt[, DMSO := no.bgrnd[Concentration == 0],
    by = list(Plate, Row)]
trt[,Normalized := 100*no.bgrnd/DMSO]

trt <- subset(trt, Concentration != 0)
trt$Concentration <- as.numeric(trt$Concentration)
summary(trt)

# Wide data
# Echinatin
dd.echin <- dcast(droplevels(subset(trt, 
                                   Treatment == "Echinatin")),
                 Plate + Repl + Concentration ~ Hour + Treatment,
                 value.var = "Normalized")
head(dd.echin)

# Curcumin
dd.curc <- dcast(droplevels(subset(trt, 
                                   Treatment == "Curcumin")),
                 Plate + Repl + Concentration ~ Hour + Treatment,
                 value.var = "Normalized")
head(dd.curc)


# SAVE DD
saved.dd <- list(echin = copy(dd.echin),
                 curc = copy(dd.curc))

#ANALYSIS I: Echinatin
m2 <- f.robicf(data.frame(Log10Conc = log10(dd.echin$Concentration),
                          dd.echin[, 4:6]),
               ecf = c(30, 40))
s2 <- summary(m2, export = TRUE)
s2
lines(m2, 
      smooth = TRUE)

par(mfrow=c(1, 3))
plot.robicf(m2,
            xlab = "Concentration (nM)",
            ylab = "% DMSO",
            smooth  = TRUE,
            xmarks = dilut.nM[-1],
            main = paste("Echinatin",
                         c(24, 48, 72), 
                         "Hours"),
            rot.x = 20)
graphics.off()

#ANALYSIS II: Curcumin
m3 <- f.robicf(data.frame(Log10Conc = log10(dd.curc$Concentration),
                          dd.curc[, 4:6]),
               ecf = c(30, 40))
s3 <- summary(m3, export = TRUE)
s3
lines(m3, 
      smooth = TRUE)
par(mfrow=c(1, 3))
plot.robicf(m3,
            xlab = "Concentration (nM)",
            ylab = "% DMSO",
            smooth  = TRUE,
            xmarks = dilut.nM[-1],
            main = paste("Curcumin",
                         c(24, 48, 72), 
                         "Hours"),
            rot.x = 20)
graphics.off()


ndx <- !do.call("c", lapply(s2, is.null))

t2 <- data.frame(Treatment = rep(names(s2)[ndx], each = 8),
                 Parameter = rep(rownames(s2[[2]]), 3),
                 do.call("rbind", s2[ndx]))
t2
write.csv(t2, "MTS4_06062016_06082016\\estimates.gen.csv")
nrow(do.call("rbind", s2[ndx]))
lines(m2,
      smooth = TRUE,
      xlab = "Log10 of Echinatin Concentration (nM)",
      ylab = "%DMSO",
      main = "All Plates, Fitted Curves (24, 48 and 72 Hour)")

ciplot(m2,
       xlab = "Log10 of Concentration (nM)")

par(mfrow=c(1, 3))
plot.robicf(m2,
            xlab = "Concentration (nM)",
            ylab = "% DMSO",
            smooth  = TRUE,
            xmarks = dilut.nM[-1],
            main = paste("Echinatin",
                         c(24, 48, 72), 
                         "Hours"),
            rot.x = 20)
graphics.off()


# #ANALYSIS II: Curcumin
# dd.dai <- subset(saved.dd, select = c(1:3, 4, 6, 8))
# dd.dai$dilut <- 90000/2^as.numeric(as.character(dd.dai$Concentration))
# 
# m1 <- f.robicf(data.frame(Log10Conc = log10(dd.dai$dilut),
#                           dd.dai[, 4:6]),
#                ecf = c(60, 70))
# s1 <- summary(m1, export = TRUE)
# s1
# 
# ndx <- !do.call("c", lapply(s1, is.null))
# 
# t1 <- data.frame(Treatment = rep(names(s1)[ndx], each = 8),
#                  Parameter = rep(rownames(s1[[2]]), 2),
#                  do.call("rbind", s1[ndx]))
# t1
# write.csv(t1, "MTS4_06062016_06082016\\estimates.dai.csv")
# nrow(do.call("rbind", s1[ndx]))
# lines(m1,
#       smooth = TRUE,
#       xlab = "Log10 of Genistein Concentration (nM)",
#       ylab = "%DMSO",
#       main = "All Plates, Fitted Curves (24, 48 and 72 Hour)")
# 
# ciplot(m1,
#        xlab = "Log10 of Concentration (nM)")
# 
# par(mfrow=c(1, 3))
# plot.robicf(m1,
#             xlab = "Daidzein Concentration (nM)",
#             ylab = "% DMSO",
#             smooth  = TRUE,
#             xmarks = dilut.nM[-1],
#             main = paste("Deidzein",
#                          c(24, 48, 72), 
#                          "Hours"))
# graphics.off()
# 
# ##########################################
# trt[, Mean := mean(Normalized),
#     by = list(Hour, Concentration, Treatment)]
# trt[, SD := sd(Normalized),
#     by = list(Hour, Concentration, Treatment)]
# trt[, Time := paste(Hour, "Hours")]
# 
# # Source:
# # http://www.cookbook-r.com/Graphs/Axes_(ggplot2)/
# ggplot(trt, 
#        aes(x = as.numeric(rev(Concentration)),
#            y = Mean,
#            fill = Treatment)) + 
#   facet_wrap(~ Time,
#              ncol = 1) +
#   geom_bar(position = position_dodge(),
#            stat="identity") +
#   geom_errorbar(aes(ymax = Mean + SD,
#                     ymin = Mean - SD),
#                 width =.2,
#                 position = position_dodge(.9)) +
#   theme(legend.position = "top",
#         axis.text.x = element_text(angle = 45, 
#                                    hjust = 1)) +
#   scale_x_continuous("Concentration (nM)",
#                    breaks = c(0:8 - 0.2, 0:8 + 0.2),
#                    labels = c(dilut.nM[-1], dilut.nM[-1])) +
#   scale_y_continuous("%DMSO") + 
#   ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
#   guides(fill = guide_legend(title = "Treatment",
#                                title.position = "top",
#                                nrow = 1))
# 
# ggplot(trt,
#        aes(x = as.numeric(rev(Concentration)),
#            y = Mean,
#            colour = Treatment,
#            group = Treatment)) + 
#   facet_wrap(~ Time,
#              ncol = 1) +
#   geom_line(position = position_dodge(0.3),
#             size = 1) +
#   geom_point(position = position_dodge(0.3),
#              size = 3) +
#   geom_errorbar(aes(ymax = Mean + SD,
#                     ymin = Mean - SD),
#                 width =.4,
#                 size = 1,
#                 position = position_dodge(0.3)) +
#   theme(legend.position = "top",
#         axis.text.x = element_text(angle = 45, 
#                                    hjust = 1)) +
#   scale_x_continuous("Concentration (nM)",
#                      breaks = c(0:8 - 0.2, 0:8 + 0.2),
#                      labels = c(dilut.nM[-1], dilut.nM[-1])) +
#   scale_y_continuous("%DMSO") + 
#   ggtitle("%DMSO Average +/- SD by Time, Treatment and Concentration") +
#   guides(fill = guide_legend(title = "Treatment",
#                              title.position = "top",
#                              nrow = 1))
# ####################################
# # NOTES:
# # 1. Did not redose at 48 Hour.