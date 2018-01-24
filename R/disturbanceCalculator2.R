#' Calculates the distrution of pathway disturbance
#'
#' @param pathwayRefGraphExpanded An KEGGGraph object of a pathway
#' @param interationNo number of rounds of sampling
#' @param deKID KEGG ID of the differentially expressed genes
#' @param allKID reference of all the genes with their respective KEGG IDs
#'
#' @return a  list of two objects. First is the sampling data for
#' causal disturbance. The second one is the actual causal disturbance
#'
#'
#'
pathSampler <- function(pathwayRefGraphExpanded,iterationNo, deKID, allKID) {



    totPathNodes     <- nodes(pathwayRefGraphExpanded)
    totGenes         <- length(totPathNodes)
    eyeTot           <- diag(totGenes)


    sizeDE       <- sum(totPathNodes %in% deKID)
    eyeDE        <- diag(totGenes - sizeDE)
    samplingData <- rep(0,iterationNo)
    pathMat      <- as(pathwayRefGraphExpanded, "matrix")
    totalPaths   <- pathCounter(pathMat,eyeTot)
    deGenesInd   <- totPathNodes %in% deKID
    deGenes      <- totPathNodes[deGenesInd]



    deMatUnRef <- pathMat[!deGenesInd,!deGenesInd]

    if (length(deMatUnRef) < 2){
        causalDisturbance <- 1
    } else{
        deTotalPathsUnRef <- pathCounter(deMatUnRef,eyeDE)
        causalDisturbance <- 1 - (deTotalPathsUnRef/totalPaths)
    }



    for (i in 1:iterationNo) {
        randPerm <- logical(totGenes)
        posPerm <- sample(1:totGenes, sizeDE,replace = F)
        randPerm[posPerm] = TRUE


        randMatUnRefSample      <- pathMat[!randPerm,!randPerm]
        totalPathsUnRefSample   <- pathCounter(randMatUnRefSample,eyeDE)
        samplingData[i]         <- 1 - (totalPathsUnRefSample / totalPaths)
    }


    return(list(samplingData, causalDisturbance))
}

