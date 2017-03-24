require(data.table)
require(ggplot2)

dt1 <- fread("DIM/data/2017-02-28_185952 DDD.csv")
dt1

# dt1 <- unique(dt1[, c(3, 4, 6:7)])
dt1 <- unique(dt1[, 3:5])
dt1$`Sample Name` <- factor(dt1$`Sample Name`)
dt1$`Target Name` <- factor(dt1$`Target Name`,
                            levels = unique(dt1$`Target Name`))
dt1$CT <- as.numeric(as.character(dt1$CT))
summary(dt1)

# Plot
ggplot(data = dt1) +
  scale_x_discrete("Treatment") + 
  scale_y_continuous("Readout") + 
  ggtitle("Title") +
  facet_wrap(~ `Target Name`,
             ncol = 2) +
  geom_boxplot(aes(x = `Sample Name`,
                   y = CT,
                   outlier.shape = NA)) +
  geom_point(aes(x = `Sample Name`,
                 y = CT),
             size = 3,
             alpha = 0.6,
             position = position_dodge(0.3))
 
  guides(colour = guide_legend(title = "ID",
                               title.position="top",
                               nrow = 1)) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))