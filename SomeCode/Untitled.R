library(PathwayDisturbance)
a.path <- PathwayDisturbance::pathways.collection[[1]]
deKID    <- translateGeneID2KEGGID(tT.de.names)
allKID   <- translateGeneID2KEGGID(tT.all.names)
iterationNo <- 5000



totPathNodes     <- nodes(a.path)
totGenes         <- length(totPathNodes)
eyeTot           <- diag(totGenes)

sizeDE       <- sum(totPathNodes %in% deKID)
eyeDE        <- diag(totGenes - sizeDE)
samplingData <- rep(0,iterationNo)
pathMat      <- as(a.path, "matrix")
#totalPaths   <- pathCounter.(pathMat,eyeTot,0.5)
deGenesInd   <- totPathNodes %in% deKID
deGenes      <- totPathNodes[deGenesInd]
eye          <- diag(nrow(pathMat))
deMatUnRef   <- pathMat[deGenesInd,deGenesInd]



newpath.centrality2 <- function(adj.matrix, alpha, beta){

    # few error throws to be implemented
    # the matrix is not square
    # the alpha is larger then the inverse of the largest eigen value

    eye <- diag(nrow(adj.matrix))
    cent.out  <- solve(eye - alpha * adj.matrix)
    cent.in   <- solve(eye - alpha * t(adj.matrix))
    cent.tot  <- (cent.out) + beta * (cent.in)
    return(cent.tot)

}


        centr.mat <- newpath.centrality2(pathMat, 0.4, beta = 1)
        paths.tot <- rowSums(centr.mat)
        cdist.tot <- rowSums(centr.mat[deGenesInd,])
        paths.log <- sum(log2(paths.tot))
        cdist.log <- sum(log2(cdist.tot))

        causalDisturbance <- cdist.log/paths.log

        for (i in 1:iterationNo) {
            randPerm     <- logical(totGenes)
            posPerm      <- sample(1:totGenes, sizeDE,replace = F)
            randPerm[posPerm] = TRUE


            cdist.tot.rand  <- rowSums(centr.mat[randPerm,])
            cdist.tot.rand  <- sum(log2(cdist.tot.rand))

            samplingData[i] <- (cdist.tot.rand /paths.log)
        }

        sampleDist      <- stats::ecdf(unlist(samplingData))
        disturbProb     <- 1 - sampleDist(causalDisturbance)
        iter2 <- iterationNo

        if(disturbProb == 0)  {
            disturbProb <- NA
        }



        return(disturbProb)
        # return(list(samplingData, causalDisturbance))


















pathSampler.test1(a.path,iterationNo = 2000, deKID,allKID,alpha = 0.4,statEval = 1)
hist(a.res[[2]])
a.res[[1]]

see.this <- 1- a.res[[2]]

hist(see.this)
plot(ecdf(see.this))





