# Source: file:///C:/R/R-3.3.2/library/ChIPseeker/doc/ChIPseeker.html
# Source: http://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html

source("https://bioconductor.org/biocLite.R")
biocLite("ChIPseeker")
biocLite("TxDb.Mmusculus.UCSC.mm9.knownGene")
biocLite("org.Mm.eg.db")
biocLite("ReactomePA")

require(ChIPseeker)
require(TxDb.Mmusculus.UCSC.mm9.knownGene)
require(ReactomePA)
# library(help = TxDb.Mmusculus.UCSC.mm9.knownGene)
# library(help = "ChIPseeker")

f1 <- "methylation seq/results/combined_yue_03312017.csv"
# peak <- readPeakFile(f1)
# covplot(peak, weightCol="CpG")

peakAnno <- annotatePeak(peak = f1, 
                         tssRegion = c(-3000, 3000), 
                         TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                         annoDb = "org.Mm.eg.db")

dt1 <- peakAnno@anno@elementMetadata@listData
require(data.table)
dt1 <- data.table(as.data.frame(dt1))
dt1

# 5.1 Visualize Genomic Annotation----
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)
plotDistToTSS(peakAnno, 
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

# 7.1.1 Average profiles----
promoter <- getPromoters(TxDb = TxDb.Mmusculus.UCSC.mm9.knownGene,
                         upstream = 3000, 
                         downstream = 3000)
tagMatrixList <- lapply(f1, 
                        getTagMatrix, 
                        windows = promoter)
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))
# plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")




#*******************************************
files <- getSampleFiles()
print(files)

peak <- readPeakFile(files[[1]])
peak

covplot(peak)
