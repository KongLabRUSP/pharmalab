# Header----
# Project: Echinatin Screening with HepG2-C8 Cells
# Author: Davit Sargsyan
# Created: 11/15/2016
# Modified: 01/03/2016
#****************************************************
require(data.table)
require(ggplot2)

# Data----
# 1. BCA----
bca <- fread("Echinatin/data/HepG2C8_SFN_Ech_BCA_11282016.csv",
             sep = ",",
             header = TRUE)

# BCA Maps----
bca.map <- fread("Echinatin/data/Map_HepG2C8_SFN_Ech_BCA_11282016.csv",
                 sep = ",",
                 header = TRUE)

# 2. Luciferase----
luc <- fread("Echinatin/data/HepG2C8_SFN_Ech_Luciferace_11282016.csv",
             sep = ",",
             header = TRUE)

#****************************************************
# Protein concentration in the samples----
# measured against known Albumine concentration
bca.long <- melt.data.table(data = bca,
                            id.vars = "Row",
                            variable.name = "Column",
                            value.name = "BCA")
bca.map.long <- melt.data.table(data = bca.map,
                                id.vars = "Row",
                                variable.name = "Column",
                                value.name = "Treatment")
dt1 <- merge(bca.long,
             bca.map.long,
             by = c("Row",
                    "Column"))
dt1$Treatment <- droplevels(factor(dt1$Treatment, 
                                   levels = unique(bca.map.long$Treatment)))
setkey(dt1, Treatment)

dt1.1 <- subset(dt1, !(Treatment %in% c("BCA_Only", "H2O")))
dt1.1[, rep := 1:.N, by = Treatment]

dt1.w <- dcast.data.table(data = dt1.1, 
                          value.var = "BCA",
                          Treatment ~ rep)
dt1.w
write.csv(dt1.w, file = "tmp/dt1.w.csv")

# Subtract background
bg <- mean(dt1$BCA[dt1$Treatment == "BCA_Only"])
dt1 <- subset(dt1, Treatment != "BCA_Only")
dt1$BCA <- dt1$BCA - bg

# Separate Albumine
alb <- subset(dt1, substring(Treatment, 1, 8) == "Albumine")
alb$Conc <- substring(alb$Treatment, 10, nchar(as.character(alb$Treatment)))
alb$Conc <- as.numeric(gsub(pattern =  "ug/mL",
                            replacement = "",
                            x = alb$Conc))
setkey(alb, Conc)
plot(alb$BCA ~ alb$Conc)
alb

alb[, mu := mean(BCA,
                 na.rm = TRUE), 
    by = Conc]
write.csv(alb, file = "tmp/alb.csv")

# Other than albumine
smpl <- subset(dt1, substring(as.character(Treatment), 1, 8) != "Albumine")
setkey(smpl, Treatment)
smpl

smpl[, mu := mean(BCA),
     by = Treatment]

mu.smpl <- unique(smpl[, c(4, 5), with = FALSE])
mu.smpl

#****************************************************
# Nonlinear fit----
m1 <- nls(BCA ~ a*Conc/(1+b*Conc),
          start=list(a = 1,b = 0),
          data = alb)

pp <- data.table(x = seq(from = 0, to = 2000, length = 50))
pp$y <- predict(m1, list(Conc = pp$x))

plot(alb$BCA ~ alb$Conc)
lines(x = pp$x, y = pp$y, col = "green")

# Inversed estimate
m2 <- nls(Conc ~ BCA/(a - b*BCA),
          start=list(a = 1,b = 0),
          data = alb)
ppp <- data.table(y = seq(from = 0, to = 3.75, length = 50))
ppp$x <- predict(m2, list(BCA = ppp$y))

lines(x = ppp$x, y = ppp$y, col = "red")

mu.smpl$pConc <- predict(m2, list(BCA = mu.smpl$mu))
points(mu.smpl$mu ~ mu.smpl$pConc, pch = 16)

#****************************************************
# Concentration Plot----
ggplot() +
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  scale_x_continuous("Protein Concentration (ug/mL)",
                     breaks = unique(alb$Conc),
                     labels = unique(alb$Conc),
                     limits = c(-10, 1000)) +
  scale_y_continuous(limits = c(-0.1, 2)) +
  geom_point(data = alb,
             aes(x = Conc,
                 y = BCA,
                 group = Column),
             size = 3,
             shape = 1) +
  geom_line(data = alb,
            aes(x = Conc,
                y = mu),
            size = 1,
            col = "black") +
  geom_point(data = alb,
             aes(x = Conc,
                 y = mu),
             size = 3,
             col = "black",
             shape = 18) + 
  geom_line(data = ppp,
            aes(x = x,
                y = y),
            size = 1,
            col = "red") + 
  geom_point(data = mu.smpl,
             aes(x = pConc,
                 y = mu,
                 colour = Treatment),
             size = 6,
             shape = 20)

#****************************************************
# Protein expression----
dt2 <- merge(mu.smpl, luc, by = "Treatment")
dt2$muRLU <- rowMeans(dt2[, 4:6, with = FALSE],
                      na.rm = TRUE)
dt2$SD <- apply(dt2[, 4:6, with = FALSE]/dt2$pConc,
                MARGIN = 1,
                FUN = sd,
                na.rm = TRUE)
dt2$RLUp_ug_mL <- dt2$muRLU/dt2$pConc
# Percent DMSO
dt2$Pct_DMSO <- 100*dt2$RLUp_ug_mL/(dt2$RLUp_ug_mL[dt2$Treatment == "DMSO_0.1%"])
dt2$Treatment <- factor(dt2$Treatment, levels = c("Echinatin_40uM",
                                                  "Echinatin_20uM",
                                                  "Echinatin_10uM",
                                                  "SFN_10uM",
                                                  "SFN_5uM",
                                                  "SFN_2.5uM",
                                                  "DMSO_0.1%",
                                                  "TSA_100uM_5AZA_500uM",
                                                  "No_Treatment"))

setkey(dt2, Treatment)
dt2
write.csv(dt2, 
          file = "tmp/dt2.csv",
          row.names = FALSE)

#****************************************************
# Expression Plot----
ggplot(dt2,
       aes(x = Treatment,
           y = RLUp_ug_mL,
           fill = Treatment)) +
  geom_errorbar(aes(ymin = RLUp_ug_mL, 
                    ymax = RLUp_ug_mL + SD),
                colour = "black", 
                width = .1) +
  geom_bar(position = position_dodge(),
           stat="identity") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") +
  scale_y_continuous("RLU/ug/mL") +
  scale_fill_manual(values = c(rep(c("green",
                                     "red"), 
                                   each = 3), 
                               "blue",
                               "pink",
                               "black"))

ggplot(dt2,
       aes(x = Treatment,
           y = Pct_DMSO,
           fill = Treatment)) +
  geom_bar(position = position_dodge(),
           stat="identity") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") +
  scale_y_continuous("% DMSO") +
  scale_fill_manual(values = c(rep(c("green",
                                     "red"), 
                                   each = 3), 
                               "blue",
                               "pink",
                               "black"))