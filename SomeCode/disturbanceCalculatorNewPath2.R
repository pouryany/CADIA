#' Calculates the distrution of pathway disturbance
#'
#' @param pathwayGraph2 An KEGGGraph object of a pathway
#' @param interationNo number of rounds of sampling
#' @param deKID KEGG ID of the differentially expressed genes
#' @param allKID reference of all the genes with their respective KEGG IDs
#'
#' @return a  list of two objects. First is the sampling data for
#' causal disturbance. The second one is the actual causal disturbance
#'
#'
#'
pathSampler.newPath2 <- function(pathwayGraph,iterationNo, deKID, allKID,alpha) {

    pathwayGraph2  <- pathwayGraph
    pathwayGraph2.comps <- connectedComp(pathwayGraph2)
    subgg <- unique(unlist(pathwayGraph2.comps[which((lapply(pathwayGraph2.comps
                                            ,length)) != 1)]))
    pathwayGraph2 <- subGraph(subgg,pathwayGraph2)

    totPathNodes     <- nodes(pathwayGraph2)
    totGenes         <- length(totPathNodes)
    eyeTot           <- diag(totGenes)

    sizeDE       <- sum(totPathNodes %in% deKID)
    eyeDE        <- diag(totGenes - sizeDE)
    samplingData <- rep(0,iterationNo)
    pathMat      <- as(pathwayGraph2, "matrix")
    #totalPaths   <- pathCounter.(pathMat,eyeTot,0.5)
    deGenesInd   <- totPathNodes %in% deKID
    deGenes      <- totPathNodes[deGenesInd]
    eye          <- diag(nrow(pathMat))
    deMatUnRef   <- pathMat[deGenesInd,deGenesInd]

    if (length(deMatUnRef) < 2){
            causalDisturbance <- 1
            return(1)
        } else if (sizeDE == 0){
            return(1)
        }else{

            #cdist.out <- sum(pathMat.out[deGenesInd,])
            #cdist.in  <- sum(pathMat.in[deGenesInd,])

            centr.mat <- newpath.centrality(pathMat, alpha, beta =1)
            paths.tot <- sum(centr.mat)
            cdist.tot <- sum(centr.mat[deGenesInd,])

            #deTotalPathsUnRef <- pathCounter.katz(deMatUnRef,eyeDE,0.5)
            causalDisturbance <- 1 - (cdist.tot/paths.tot)


            for (i in 1:iterationNo) {
                randPerm <- logical(totGenes)
                posPerm <- sample(1:totGenes, sizeDE,replace = F)
                randPerm[posPerm] = TRUE


                cdist.tot.rand <- sum(centr.mat[randPerm,])
                # cdist.out.rand  <- sum(pathMat.out[randPerm,])
                # cdist.in.rand   <- sum(pathMat.in[randPerm,])
                # cdist.tot.rand  <- cdist.in.rand + cdist.out.rand
                samplingData[i] <- 1 - (cdist.tot.rand / paths.tot)
            }

            sampleDist      <- stats::ecdf(unlist(samplingData))
            disturbProb     <- 1 - sampleDist(causalDisturbance)
            iter2 <- iterationNo

            # while(disturbProb == 0 && iter2 < 6 * iterationNo){
            #
            #      for (i in iter2:(iter2 +iterationNo))
            #         {
            #          randPerm <- logical(totGenes)
            #          posPerm <- sample(1:totGenes, sizeDE,replace = F)
            #          randPerm[posPerm] = TRUE
            #
            #          cdist.tot.rand <- sum(centr.mat[randPerm,])
            #          # cdist.out.rand  <- sum(pathMat.out[randPerm,])
            #          # cdist.in.rand   <- sum(pathMat.in[randPerm,])
            #          # cdist.tot.rand  <- cdist.in.rand + cdist.out.rand
            #           samplingData[i] <- 1 - (cdist.tot.rand / paths.tot)
            #         }
            #
            #
            #
            #     iter2         <- iter2 +iterationNo
            #     sampleDist    <- stats::ecdf(unlist(samplingData))
            #     disturbProb   <- 1 - sampleDist(causalDisturbance)
            # }
            #
            if(disturbProb == 0)  {
                disturbProb <- 1/iterationNo
            }



        return(disturbProb)
       # return(list(samplingData, causalDisturbance))
    }
}

