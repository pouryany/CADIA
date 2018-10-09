library(KEGGgraph)
library(Rgraphviz)

if(!require("hgu133plus2.db")){
        biocLite("hgu133plus2.db")
    } else {
        library(hgu133plus2.db)
    }

if(!require("SPIA")){
    biocLite("SPIA")
    } else{
        library("SPIA")
    }

library("RBGL")

if(!require("igraph")){
      install.packages("igraph")
    } else {
        library(igraph)
    }


if(!require("annotate")){
        biocLite("annotate")
    } else {
        library(annotate)
    }





# Loading pathways to be analyzed

source("causalDisruption2.R")
pathList<- read.csv("pathlists.csv", header = F)
pathList <- as.vector(pathList)

for( i in pathList){
    temp <- paste(i,".xml", sep="")
    temp <- paste("0", temp, sep="")
    }
pathList <- temp


# This one is a corrupt xml file
pathList[163]
pathList <- pathList[-163]
#pathList.names <- pathList.names[-163]


pathways.collection <- list()
pathway.name <- list()
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

    pathways.compInfo [i] <- max(sapply(connectedComp(pathways.collection[[i]]) , length))
}

pathways.compInfo <- unlist(pathways.compInfo)

pathways.comps <- lapply(pathways.collection,connectedComp)
xx <- lapply(pathways.comps, function(X) lapply(X, length))
pathways.isolated <-sapply(xx, function (X) (length(which(X ==1))))






pathways.summary <- cbind(pathways.nodeInfo,pathways.edgeInfo
                          ,pathways.compInfo,pathways.isolated)
pathways.summary <- as.matrix(pathways.summary)
a <- as.data.frame(pathways.summary)


# Few filtering lines
b <- rownames(a[(a$pathways.edgeInfo < 0.2 * a$pathways.nodeInfo),])
c <- rownames(a[(a$pathways.edgeInfo < 20),])
d <- rownames(a[(a$pathways.isolated > 0.5 * a$pathways.nodeInfo),])
e <- rownames(a[a$pathways.compInfo < 10,])
f <- Reduce(union, list(b,c,d,e))
length(f)
#intersect(e,d)

pathways.collection       <- pathways.collection[rownames(a)]
pathways.collection.names <- names(pathways.collection)


pathNames <- list()
j <- 1
for (i in pathways.collection.names){
    pathway <- KEGGgraph::parseKGML(i)
    pathNames[i] <- KEGGgraph::getTitle(pathway)
    j <- j +1
}
pathways.collection.names <- unlist(pathNames)
length(pathways.collection.names)


zzz <- cbind(a[f,], pathways.collection.names[f])
write.table(zzz,file = "badpathways.csv")
f <- rownames(a) %in% f

pathways.collection <-  pathways.collection[!f]
pathways.collection.names<- pathways.collection.names[!f]

zzz <- cbind(a[!f,],pathways.collection.names)
write.table(zzz,file = "GoodPathways.csv")


#These commands create a file which is then moved to the package directory
save(pathways.collection, file = "pathwaysData")
save(pathways.collection.names, file = "pathwayDataNames")



length(pathways.collection)
length(pathways.collection.names)

# 30 pathways were removed because of the connected component
# Commands for updating data in the package.

#load(file = "pathwaysData")
#load(file = "pathwayDataNames")
#devtools::use_data(pathways.collection,overwrite = T)
#devtools::use_data(pathways.collection.names, overwrite = T)



library(RBGL)

mtx.collection <-sapply(pathways.collection, function(X)(as(X,"matrix")))

eigen.collection <- lapply(mtx.collection, eigen, only.value = T )




library(dplyr)

largest.eigen <- function(b) {
    c <-  b %>% unlist()  %>% Im() == 0
    c %>% subset.default(x = b) %>% Re() %>% max()
}

largest.eigen2 <- function(b) {
    c <-  b %>% unlist()  %>% Re() %>% max()
}



a.test <- unlist(eigen.collection[[148]])
largest.eigen(unlist(a.test))
a.test  %>% unlist() %>% Re() %>% max()
a.test %>% subset.default(x =  unlist(eigen.collection[[148]])) %>% Re() %>% max()

d <- sapply(eigen.collection,function(X)(largest.eigen2(unlist(X))))
length(d)

pathways.collection.names[names(d[d>=10])]





