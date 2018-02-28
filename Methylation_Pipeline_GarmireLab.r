
source("https://bioconductor.org/biocLite.R")
biocLite("ChAMP")

source("http://bioconductor.org/biocLite.R")
biocLite(c("minfi","ChAMPdata","Illumina450ProbeVariants.db","sva","IlluminaHumanMethylation450kmanifest",
           "limma","RPMM","DNAcopy","preprocessCore","impute","marray","wateRmelon","goseq","plyr","GenomicRanges","RefFreeEWAS","qvalue","isva","doParallel","bumphunter","quadprog","shiny","shinythemes","plotly","RColorBrewer","DMRcate","dendextend",
           "IlluminaHumanMethylationEPICmanifest","FEM","matrixStats","missMethyl","combinat"))

library("ChAMP")

library("doParallel")
detectCores()

sample_annotation=read.csv(file="6_chips_EPIC/sample_annotation_EPIC_R01.csv",
                          colClasses=c('factor', 'factor', 'factor','factor','numeric','factor','factor','factor',
                                      'numeric','factor','factor','factor'),check.names=F)
head(sample_annotation)
table(sample_annotation$Parity)
table(sample_annotation$Delivery_route)
table(sample_annotation$Parity)
table(sample_annotation$Sample_Group)
table(sample_annotation$Alcohol)

# myLoad_M_Texas <- champ.load(directory="2_chips_EPIC_Texas",arraytype="EPIC",filterBeads=T,methValue="M")# 2 Chips from Texas 
# myLoad_B_Texas <- champ.load(directory="2_chips_EPIC_Texas",arraytype="EPIC",filterBeads=T,methValue="B")# 2 Chips from Texas

# myLoad_M_2 <- champ.load(directory="2_chips_EPIC",arraytype="EPIC",filterBeads=T,methValue="M")# The first 2 chips from Jashwa 
# myLoad_B_2 <- champ.load(directory="2_chips_EPIC",arraytype="EPIC",filterBeads=T,methValue="B")# The first 2 chips from Jashwa

myLoad_B <- champ.load(directory="6_chips_EPIC",arraytype="EPIC",filterBeads=T,methValue="B")# All six chips beta values
#myLoad_M <- champ.load(directory="6_chips_EPIC",arraytype="EPIC",filterBeads=T,methValue="M")# All six chips M values

# myLoad_M_2_oth <- champ.load(directory="2_chips_EPIC_other2",arraytype="EPIC",filterBeads=T,methValue="M")# The OTHER 2 chips run with last 4 chips 
# myLoad_B_2_oth <- champ.load(directory="2_chips_EPIC_other2",arraytype="EPIC",filterBeads=T,methValue="B")# The OTHER 2 chips run with last 4 chips

#             champ.load(directory = getwd(),
#                             methValue="B",
#                             filterDetP=TRUE,
#                             detSamplecut=0.1, #sample cut
#                             detPcut=0.01,     #probs cut
#                             removeDetP = 0,
#                             filterBeads=TRUE,
#                             beadCutoff=0.05,  #beads
#                             filterNoCG=TRUE,
#                             filterSNPs=TRUE,
#                             population=NULL,  #for snp filtration 
#                             filterMultiHit=TRUE,
#                             filterXY=TRUE,
#                             arraytype="450K")



myLoad_B$pd

dim(myLoad_B$beta)
head(myLoad_B$beta)

champ.QC(beta = myLoad_B$beta, 
                  pheno=myLoad_B$pd$Sample_Group,
                  mdsPlot=TRUE,
                  densityPlot=T,
                  dendrogram=F,
                  PDFplot=TRUE,
                  Rplot=TRUE,
                  Feature.sel="None",
                  resultsDir="./CHAMP_QCimages/") 


library(RColorBrewer)

numPositions=1000
b <- myLoad_B$beta
#b <- myNorm
o <- order(-matrixStats::rowVars(b))[1:numPositions]
    d <- dist(t(b[o, ]),method = "euclidean")
    fit <- cmdscale(d)
sampGroups <- as.factor(myLoad_B$pd$Sample_Group)
numGroups <- length(levels(sampGroups))
legendNCol <- numGroups
xlim <- range(fit[, 1]) * 1.2
ylim <- range(fit[, 2]) * 1.2
pal = brewer.pal(8, "Dark2")
col <- pal[sampGroups]
main <- sprintf("MDS\n%d most variable positions", 
                numPositions)
legendPos="bottomleft"
plot(fit[, 1], fit[, 2], col = col, pch = 2, xlim = xlim, 
            ylim = ylim, xlab = "", ylab = "", main = main)
legend(legendPos, legend = levels(sampGroups), ncol = legendNCol, 
            text.col = pal[1:numGroups])

#prepare the data
bb=cbind(t(b),as.factor(myLoad_B$pd$Sample_Group))
colnames(bb)[ncol(bb)]="label"
bb=data.frame(bb)
bb$label=as.factor(bb$label)
levels(bb$label)=c("control","case")

#PCA for only 1000 cpgs
library(ggplot2)
library(ggfortify)
bbb=bb[,-ncol(bb)]
autoplot(prcomp(bbb[,o]),colour = 'label',data=bb)

#install.packages("Rtsne")
library(Rtsne)
# tSNE for only 1000
 iris_matrix <- as.matrix(bbb[,o])
 set.seed(42) # Set a seed if you want reproducible results
 tsne_out <- Rtsne(iris_matrix,perplexity=4,max_iter = 1000) # Run TSNE

tsne_plot <- data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2], col = bb$label)
ggplot(tsne_plot) + geom_point(aes(x=x, y=y, color=col))

# # Show the objects in the 2D tsne representation
 #plot(tsne_out$Y,col=iris_unique$Species)

#getwd()
#ls()
list.files(path = getwd())

myNorm <- champ.norm(beta=myLoad_B$beta,
                    rgSet=myLoad_B$rgSet,
                    mset=myLoad_B$mset,
                    resultsDir="./CHAMP_Normalization/",
                    method="PBC",  #BMIQ
                    plotBMIQ=TRUE,
                    arraytype="EPIC",
                    cores=50)



#source('champ.SVD1.R')



champ.SVD(beta = myNorm,
                   rgSet=myLoad_B$rgSet,
                   pd=myLoad_B$pd,
                   RGEffect=FALSE,
                   PDFplot=TRUE,
                   Rplot=TRUE,
                   resultsDir="./CHAMP_SVDimages/")


  myCombat <- champ.runCombat(beta=myNorm,
                         pd=myLoad_B$pd,
                         variablename="Sample_Group",
                         batchname=c("Slide"),
                         logitTrans=TRUE)


champ.SVD(beta = myCombat,
                   rgSet=myLoad_B$rgSet,
                   pd=myLoad_B$pd,
                   RGEffect=FALSE,
                   PDFplot=TRUE,
                   Rplot=TRUE,
                   resultsDir="./CHAMP_SVDimages/")

myDMP <- champ.DMP(beta = myNorm,pheno=myLoad_B$pd$Sample_Group,adjPVal = 0.05,
                   adjust.method = "BH",arraytype = "EPIC")
head(myDMP)

myDMR <- champ.DMR(beta=myNorm,pheno=myLoad_B$pd$Sample_Group,
                   arraytype = "EPIC", method="Bumphunter",cores=10)
head(myDMR$BumphunterDMR)



qc <- getQC(myLoad_B$mset)
dim(qc)
plotQC(qc)




