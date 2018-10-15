
install.packages("gplots")
install.packages("corrplot")
library(gplots)
library(corrplot)

dat <- read.csv("SFNUA/comprehensive table.csv")

dim(dat[dat$RA.vs.C!=dat$SFN.vs.C,])
dim(dat[dat$RA.vs.C!=dat$UA.vs.C,])
dim(dat[dat$RA.vs.C!=dat$RA.SFN.vs.RA,])
dim(dat[dat$RA.vs.C!=dat$RA.UA.vs.RA,])
dim(unique(dat))
length(unique(dat$RA.vs.C))


min <- apply(dat[, 7:11], 1, FUN=min)
max <- apply(dat[, 7:11], 1, FUN=max)
hist(max-min)
sum(log(max/min)<2)
plot(min, max)
abline(a=0, b=1)
abline(a=1, b=1)
abline(a=2, b=1)

# cluster data
clustdata <- dat[, 7:11]
rownames(clustdata) <- dat$RA.vs.C

# filtering
rowmax <- apply(abs(clustdata), 1, max)
summary(rowmax)
clustdata2 <- clustdata[rowmax > 4, ]

dist <- dist(clustdata2) # method="man" # is a bit better
hc_sec <- hclust(dist, method = "complete")
plot(hc_sec, main="", xlab="", sub="", ylab="", frame.plot=T, hang=-5, cex=0.8)


# heatmap
heatmapdata <- as.matrix(clustdata2)
hm <- heatmap.2(heatmapdata, cexCol=0.6, cexRow = 0.4)
hm.output <- t(hm$carpet)
hm.output <- as.data.frame(hm.output[nrow(hm.output):1, ])
hm.output
write.csv(hm.output, "SFNUA/Heatmap Output.csv")


# compare rows:
sum(dat$RA.vs.C != dat$SFN.vs.C)
sum(dat$RA.vs.C != dat$UA.vs.C)
sum(dat$RA.vs.C != dat$RA.SFN.vs.RA)
sum(dat$RA.vs.C != dat$RA.UA.vs.RA)

# too big to plot
plot(RA.vs.C ~ SFN.vs.C, data=dat)

dat$diff <- dat$RA.SFN.vs.RA.1 - dat$RA.UA.vs.RA.1
dat$direction <- dat$diff>0
dat$direction[dat$diff==0] = 0
dat$direction[dat$diff<0] = -1
table(dat$direction, useNA = "ifany")


