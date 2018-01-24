library(KEGGgraph)
library(Rgraphviz)
library(hgu133plus2.db)
library("SPIA")
library("RBGL")
library(igraph)
library(annotate)


data(colorectalcancer)
x <- hgu133plus2ENTREZID
top$ENTREZ<-unlist(as.list(x[top$ID]))
top<-top[!is.na(top$ENTREZ),]
top<-top[!duplicated(top$ENTREZ),]
tg1<-top[top$adj.P.Val<0.1,]
DE_Colorectal=tg1$logFC
names(DE_Colorectal)<-as.vector(tg1$ENTREZ)
ALL_Colorectal=top$ENTREZ



deKID <- translateGeneID2KEGGID(names(DE_Colorectal))
allKID <- translateGeneID2KEGGID(ALL_Colorectal)






res=spia(de=DE_Colorectal,all=ALL_Colorectal,organism="hsa",nB=2000,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)
aa <- ggplot(someStats, aes(x = someStats$V1 , y = someStats$V2)) + geom_point() + geom_smooth(span = 0.1)




#some plotting










### Test for the null distribution


source("junkForProcess/causalDisruption2.R")
pathList <- read.csv("junkForProcess/pathlists.csv", header = F)
pathList <- as.vector(pathList)
for( i in pathList){
temp <- paste(i,".xml", sep="")
temp <- paste("0", temp, sep="")
}
pathList <- temp

pathList
pathListLocal <- paste("data-raw/",pathList,sep = "")
keggRefGenes <- geneKEGGListExtractor(pathListLocal)
pathListLocal[163]
pathList <- pathList[-163]


pathRefClean <- pathList

pathRefClean

oraMiss <- list()
cdistMiss <- list()



for(j in 1:100){

  deKIDRand <- sample(allKID,500, replace = F)


  randTest <- list()


  for ( i in pathRefClean ) {
    randTest <- rbind( randTest ,unlist(processPathway(i, deKIDRand, allKID, keggRefGenes,1000 )))
  }


  randTest <- as.data.frame(randTest)

  randTestClean <-randTest[randTest$V6 != 0,]
  randTestClean$Adj.pvalue <- p.adjust(randTestClean$V10, method = "BY")
  randTestClean$Adj.ORA <- p.adjust(randTestClean$V5, method = "BY")

  #randTestClean[randTestClean$V8 < 0.05,]
  #randTestClean[randTestClean$V5 < 0.05,]
  #sum(randTestClean$Adj.pvalue < 0.05)
  #sum(randTestClean$Adj.ORA < 0.05)


  oraMiss[j]   <-  sum(randTestClean$Adj.ORA < 0.05)
  cdistMiss[j] <-  sum(randTestClean$Adj.pvalue < 0.05)
  print("round Done")


}

unlist(oraMiss)
cdistMiss

falseError <- unlist(cdistMiss)
hist(falseError)
sum(falseError)/(100 * 137)

length(pathRefClean)
hist(as.numeric(randTestClean$V8), xlab="p-values",
     main="Histogram of Cdist Pvalues", prob=TRUE,
     cex.lab=1.5, cex.axis=1.5, cex.main=2, cex.sub=2)
hist(as.numeric(randTestClean$V5))


library(ggplot2)
ggplot(as.data.frame(as.numeric(unlist(randTestClean$V8))), aes(as.data.frame(as.numeric(unlist(randTestClean$V8))))) + geom_histogram(binwidth=0.1)



randTestClean






res2 <- list()
for ( i in pathListLocal ) {
  res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,500 )))
}


res2
res3 <- res2[which(res2[,10] != 0),]
res3 <- as.data.frame(res3)
res3$Adj.Pvalue <- p.adjust(res3$causalDisturbance, method = "BY")
res4 <- res3[res3$Adj.Pvalue < 0.05 , ]
res3$Adj.ORA <- p.adjust(res3$fTestRes, method = "BY")
res5 <- res3[res3$Adj.ORA < 0.05, ]









res2 <- list()
for ( i in pathList ) {
res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,5000 )))
}
res2
res3 <- res2[which(res2[,10] != 0),]
res3 <- as.data.frame(res3)
res3$Adj.Pvalue <- p.adjust(res3$causalDisturbance, method = "BY")
res4 <- res3[res3$Adj.Pvalue < 0.05 , ]
res3$Adj.ORA <- p.adjust(res3$fTestRes, method = "BY")
res5 <- res3[res3$Adj.ORA < 0.05, ]

















## Below code is for sample study


x <- read.csv("ANOVARES.csv", header = T)
X
x
top <- x
top$adj.p.val <- p.adjust(top$p.value.Malignant.vs..Normal., method = "BY")
top <- top[!is.na(top$Entrez.Gene),]
top
whic(top$Entrez.Gene == "---")
which(top$Entrez.Gene == "---")
which(top$Entrez.Gene != "---")
top <- top[which(top$Entrez.Gene != "---"),]
top
tg1 <- top[top$adj.p.val < 0.05,]
tg1
abs(tg1$Fold.Change.Malignant.vs..Normal. > 2)
which(abs(tg1$Fold.Change.Malignant.vs..Normal. > 2))
ref <- abs(tg1$Fold.Change.Malignant.vs..Normal. > 2)
ref <- as.logical(ref)
ref
tg1 <- tg1[ref,]
tg1
tail(tg1)
min(tg1$Fold.Change.Malignant.vs..Normal.)
tg1 <- top[top$adj.p.val < 0.05,]
min(tg1$Fold.Change.Malignant.vs..Normal.)
which(abs(tg1$Fold.Change.Malignant.vs..Normal.) > 2)
tg1 <- tg1[which(abs(tg1$Fold.Change.Malignant.vs..Normal.) > 2),]
tg1
pathList <- read.csv("pathlists.csv", header = F)
pathList <- as.vector(pathList)
for( i in pathList){
temp <- paste(i,".xml", sep="")
temp <- paste("0", temp, sep="")
}
pathList <- temp
deKID <- translateGeneID2KEGGID(tg1$Entrez.Gene)
deKID
top
which(top$Entrez.Gene == " ")
which(top$Entrez.Gene == "")
top <- top[!which(top$Entrez.Gene == ""), ]
top
top <- x
top$adj.p.val <- p.adjust(top$p.value.Malignant.vs..Normal., method = "BY")
top <- top[!is.na(top$Entrez.Gene),]
top
top <- top[which(top$Entrez.Gene != ""), ]
top
allKID <- translateGeneID2KEGGID(top$Entrez.Gene)
allKID
source("causalDisruption2.R")
res2 <- list()
for ( i in pathList ) {
res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,3000 )))
}
library(RBGL)
res2 <- list()
for ( i in pathList ) {
res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,3000 )))
}
keggRefGenes <- geneKEGGListExtractor(pathList)
which(pathList == "05418.xml")
pathList <- pathList[-163]
pathList
res2 <- list()
for ( i in pathList ) {
res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,3000 )))
}
res2
View(res2)
res2 <- res2[which(res2$V8 != 0),]
res2
source("causalDisruption5.R")
res2 <- list()
for ( i in pathList ) {
res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,3000 )))
}
library(igraph)
source("causalDisruption5.R")
res2 <- list()
for ( i in pathList ) {
res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,3000 )))
}
res2
res2 <- as.data.frame(res2)
res2
res2 <- res2  [which(res2$V8 != 0), ]
res2
res2$p.adj <- p.adjust(res2$V10, method = "BY")
res3 <- res2[res2$p.adj < 0.05,]
res3
res3 <- res2[res2$p.adj < 0.25,]
res3
res2$ora.adj <- p.adjust(res2$V5, method = "BY")
res4 <- res2[res2$ora.adj < 0.25,]
res2
res3
res4
res4 <- res2[res2$ora.adj < 0.1,]
res3 <- res2[res2$p.adj < 0.1,]




#some other fixin codesource("causalDisruption6.R")
res2 <- list()
for ( i in pathList ) {
res2 <- rbind( res2 ,unlist(processPathway(i, deKID, allKID, keggRefGenes,10000 )))
}
View(res2)
res2 <- as.data.frame(res2)
res2 <- res2[res2$V8 !=0, ]
res2$adj.Dist <- p.adjust(res2$V10, method = "BY")
res2$adj.ORA  <- p.adjust(res2$V5, method = "BY")





i <- 1
pathList2 <- list()
for (pathway in pathList){

    pathList2[[pathway]] <-  parseKGML2Graph(file = paste('data-raw/',pathway , sep = ""))

}

pathList2[[1]]


pathList[2]

pathNames <- list()
i <- 1
for (path in pathList){
   pathNames[i] <- KEGGgraph::getTitle(parseKGML(file = paste('data-raw/',path , sep = "")))
   i <- i+1
   }
