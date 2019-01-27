#' Calculates the distrution of pathway disturbance
#'
#' @param pathwayRefGraphExpanded An KEGGGraph object of a pathway
#' @param interationNo number of rounds of sampling
#' @param deKID KEGG ID of the differentially expressed genes
#' @param allKID reference of all the genes with their respective KEGG IDs
#' @param iter the number of iterations for causal disturbance
#' @param alpha the dampening factor for Source/Sink Centrality
#' @param beta the relative Source vs Sink factor
#' @param statEval Choose 1 for product-based, 0 for summation-based
#'
#' @return a  list of two objects. First is the sampling data for
#' causal disturbance. The second one is the actual causal disturbance
#'
#'
#' @export
pathSampler.newPath <- function(pathwayRefGraphExpanded,iterationNo, deKID,
                                allKID,alpha,beta, statEval) {


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


            centr.mat <- newpath.centrality(pathMat, alpha, beta)
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
    else if (statEval == 1) {

        if (length(deMatUnRef) < 2)
        {
            causalDisturbance <- 1
            return(1)
        } else if (sizeDE == 0){
            return(1)
        }else{

            tryCatch({


            centr.mat <- newpath.centrality(pathMat, alpha, beta)
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
                disturbProb <- 1/iterationNo
            }



            return(disturbProb)
            # return(list(samplingData, causalDisturbance))
            }, error = function(e) {
                disturbProb <- NA
                return(disturbProb)
            })
            }
    }

}

