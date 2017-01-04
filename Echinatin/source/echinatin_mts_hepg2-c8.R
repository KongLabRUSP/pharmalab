# Header----
# Project: SFN and Echinatin Effect on HepG2-C8 Cells
# Author: Davit Sargsyan
# Created: 11/15/2016
#****************************************************
require(data.table)
require(ggplot2)

# Data
# 1. MTS
mts24 <- fread("Echinatin/data/HepG2C8_SFN_Ech_24h_MTS_09262016.csv",
             sep = ",",
             header = TRUE)
mts24

mts48 <- fread("Echinatin/data/HepG2C8_SFN_Ech_48h_MTS_09272016.csv",
             sep = ",",
             header = TRUE)
mts48

mts72 <- fread("Echinatin/data/HepG2C8_SFN_Ech_72h_MTS_09282016.csv",
               sep = ",",
               header = TRUE)
mts72

mts.map <- fread("Echinatin/data/Map_HepG2C8_SFN_Ech_MTS_09262016.csv",
             sep = ",",
             header = TRUE)
mts.map

#****************************************************
# Efficacy----
# 24 hours
# mts.long <- melt.data.table(data = mts24,
#                             id.vars = "Row",
#                             variable.name = "Column",
#                             value.name = "MTS")

# # 48 hours
# mts.long <- melt.data.table(data = mts48,
#                             id.vars = "Row",
#                             variable.name = "Column",
#                             value.name = "MTS")

#72 hours
mts.long <- melt.data.table(data = mts72,
                            id.vars = "Row",
                            variable.name = "Column",
                            value.name = "MTS")

# Maping file
mts.map.long <- melt.data.table(data = mts.map,
                                id.vars = "Row",
                                variable.name = "Column",
                                value.name = "Treatment")
dt0 <- merge(mts.long,
             mts.map.long,
             by = c("Row",
                    "Column"))
setkey(dt0, Treatment)
dt0

bg <- mean(dt0$MTS[dt0$Treatment == "MTS"])
dt0 <- subset(dt0, Treatment != "MTS")
dt0$MTS <- dt0$MTS - bg

dmso <- mean(dt0$MTS[dt0$Treatment == "DMSO_01%"])

dt0 <- subset(dt0, Treatment != "DMSO_01%")
dt0$pct.dmso <- 100*dt0$MTS/dmso

dt0$trt <- dt0$conc <- NA

dt0$trt[substr(dt0$Treatment, 1, 3) == "SFN"] <- 
  substr(dt0$Treatment,
         1, 
         3)[substr(dt0$Treatment, 1, 3) == "SFN"]
dt0$conc[substr(dt0$Treatment, 1, 3) == "SFN"] <- 
  as.numeric(substr(dt0$Treatment, 
                    5, 
                    nchar(dt0$Treatment) - 2)[substr(dt0$Treatment, 1, 3) == "SFN"])


dt0$trt[substr(dt0$Treatment, 1, 3) != "SFN"] <- 
  substr(dt0$Treatment, 1, 9)[substr(dt0$Treatment, 1, 3) != "SFN"]
dt0$conc[substr(dt0$Treatment, 1, 3) != "SFN"] <- 
  as.numeric(substr(dt0$Treatment, 
                    11, 
                    nchar(dt0$Treatment) - 2)[substr(dt0$Treatment, 1, 3) != "SFN"])
setkey(dt0, conc)
setkey(dt0, trt)
dt0[, muMTS:= mean(pct.dmso), 
    by = Treatment]

dt0$lConc <- log(dt0$conc)
dt0

#****************************************************
# Efficacy plot----
ggplot(dt0) +
  facet_grid(.~ trt) +
  scale_y_continuous("% DMSO", limits = c(0, 120)) +
  geom_point(aes(x = lConc,
                 y = pct.dmso,
                 colour = Column,
                 group = Column),
             size = 3) +
  geom_line(aes(x = lConc,
                y = muMTS),
            size = 1,
            col = "black") +
  theme(axis.text.x = element_text(angle = 45,
                                   hjust = 1),
        legend.position = "none") +
  scale_x_continuous("Drug Concentration (uM)",
                     breaks = unique(dt0$lConc),
                     labels = round(unique(dt0$conc/1000), 2))