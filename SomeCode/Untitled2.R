pathSampler.test1 <- function(pathwayRefGraphExpanded,iterationNo, deKID,
                                allKID,alpha, statEval) {


    totPathNodes     <- nodes(pathwayRefGraphExpanded)
    totGenes         <- length(totPathNodes)
    eyeTot           <- diag(totGenes)

    sizeDE       <- sum(totPathNodes %in% deKID)
    eyeDE        <- diag(totGenes - sizeDE)
    samplingData <- rep(0,iterationNo)
    pathMat      <- as(pathwayRefGraphExpanded, "matrix")
    #totalPaths   <- pathCounter.(pathMat,eyeTot,0.5)
    deGenesInd   <- totPathNodes %in% deKID
    deGenes      <- totPathNodes[deGenesInd]
    eye          <- diag(nrow(pathMat))
    deMatUnRef   <- pathMat[deGenesInd,deGenesInd]



    if (statEval == 0) {

        if (length(deMatUnRef) < 2)
        {
            causalDisturbance <- 1
            return(1)
        } else if (sizeDE == 0){
            return(1)
        }else{


            centr.mat <- newpath.centrality(pathMat, alpha, beta = 1)
            paths.tot <- sum(centr.mat)
            cdist.tot <- sum(centr.mat[deGenesInd,])
            causalDisturbance <- cdist.tot/paths.tot

            for (i in 1:iterationNo) {
                randPerm     <- logical(totGenes)
                posPerm      <- sample(1:totGenes, sizeDE,replace = F)
                randPerm[posPerm] = TRUE


                cdist.tot.rand <- sum(centr.mat[randPerm,])
                samplingData[i] <- (cdist.tot.rand / paths.tot)
            }

            sampleDist      <- stats::ecdf(unlist(samplingData))
            disturbProb     <- 1 - sampleDist(causalDisturbance)
            iter2 <- iterationNo

            if(disturbProb == 0)  {
                disturbProb <- NA
            }



            return(disturbProb)
            # return(list(samplingData, causalDisturbance))
        }
    }



    if (statEval == 1) {

        if (length(deMatUnRef) < 2)
        {
            causalDisturbance <- 1
            return(1)
        } else if (sizeDE == 0){
            return(1)
        }else{


            centr.mat <- newpath.centrality(pathMat, alpha, beta = 1)
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

                samplingData[i] <- (cdist.tot.rand / paths.log)
            }

            sampleDist      <- stats::ecdf(unlist(samplingData))
            disturbProb     <- 1 - sampleDist(causalDisturbance)
            iter2 <- iterationNo

            if(disturbProb == 0)  {
                disturbProb <- NA
            }



            return(disturbProb)
            # return(list(samplingData, causalDisturbance))
        }
    }

}
