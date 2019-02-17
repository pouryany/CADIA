

library(KEGGgraph)
library(Rgraphviz)
if(!require("hgu133plus2.db")){
    biocLite("hgu133plus2.db")
}
library(hgu133plus2.db)


if(require("SPIA") == 0){
    biocLite("SPIA")
}
library("SPIA")
library("RBGL")
if(!require("igraph")){
  install.packages("igraph")
}
library(igraph)

if(!require("annotate")){
    biocLite("annotate")
}

library(annotate)



pathList<- read.csv("pathlists.csv", header = F)
pathList <- as.vector(pathList)
for( i in pathList){
temp <- paste(i,".xml", sep="")
temp <- paste("0", temp, sep="")
}
pathList <- temp


# This is an empty/corrupted file
pathList[163]
pathList <- pathList[-163]
pathList.names <- pathList.names[-163]
pathList

pathways.collection <- list()

for (i in 1:length(pathList)){
    pathways.collection[[i]] <- parseKGML2Graph(pathList[i])
    print(cat(i, "\n"))
}

names(pathways.collection) <- pathList

pathways.collection <- lapply(pathways.collection, removeSelfLoops)
pathways.nodeInfo   <- unlist(lapply(pathways.collection, numNodes))
pathways.edgeInfo   <- unlist(lapply(pathways.collection, numEdges))
pathways.compInfo   <- list()



for ( i in 1: length(pathways.collection)){

    pathways.compInfo [i] <- max(sapply(connectedComp(pathways.collection[[i]]) ,
                                        length))
}

pathways.compInfo <- unlist(pathways.compInfo)

pathways.comps <- lapply(pathways.collection,connectedComp)
xx <- lapply(pathways.comps, function(X) lapply(X, length))
pathways.isolated <-sapply(xx, function (X) (length(which(X ==1))))






pathways.summary <- cbind(pathways.nodeInfo,pathways.edgeInfo
                          ,pathways.compInfo,pathways.isolated)
pathways.summary <- as.matrix(pathways.summary)
a <- as.data.frame(pathways.summary)

#a <-  a[!(a$pathways.edgeInfo < 0.1 * a$pathways.nodeInfo),]


b <- rownames(a[(a$pathways.edgeInfo < 0.2 * a$pathways.nodeInfo),])
c <- rownames(a[(a$pathways.edgeInfo < 20),])
d <- rownames(a[(a$pathways.isolated > 0.7 * a$pathways.nodeInfo),])
#d <- rownames(a[(a$pathways.compInfo < 0.1  * a$pathways.nodeInfo),])
e <- rownames(a[a$pathways.compInfo < 10,])
f <-Reduce(union, list(b,c,d,e))
length(d)
#intersect(e,d)

pathways.collection.names[pathways.collection.names %in% f]
zzz <- cbind(a[f,], pathways.collection.names[f])
write.table(zzz,file = "badpathways.csv")
f <- rownames(a) %in% f
sum(f)
a <-  a[!f,]
a[rownames(a) == "04978.xml",]
pathways.collection <- pathways.collection[rownames(a)]

save(pathways.collection, file = "pathwaysData")
pathways.collection.names <- names(pathways.collection)

load(file = "pathwaysData")



pathNames <- list()
j <- 1
 for (i in pathways.collection.names){
     pathway <- KEGGgraph::parseKGML(i)
     pathNames[i] <- KEGGgraph::getTitle(pathway)
     j <- j +1
 }
pathways.collection.names <- unlist(pathNames)
save(pathways.collection.names, file = "pathwayDataNames")
length(pathways.collection.names)
pathways.collection

#

# 30 pathways were removed because of the connected component
# Commands for updating data in the package.

load(file = "pathwaysData")
load(file = "pathwayDataNames")
devtools::use_data(pathways.collection,overwrite = T)
devtools::use_data(pathways.collection.names, overwrite = T)




