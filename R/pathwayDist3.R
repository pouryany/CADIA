#' Calculates some information regarding the pathways
#'
#' @param pGraph a GraphNEL object containing the pathways
#' @param pathID KEGG ID of an input pathway
#' @param deKIDs list of differentially expressed genes
#' @param Name the pathway name
#' @param allKIDs list of all genes with their associated KEGG IDs
#' @param keggRef list of reference genes from KEGG
#' @param iter the number of iterations for causal disturbance
#'
#'
#'
#'

processPathway <- function(pGraph,Name, deKIDs, allKIDs, keggRef, iter) {

    tempNodes   <- nodes(pGraph) %in% allKIDs
    nodeNums    <- nodes(pGraph)[tempNodes]
    isDiffExp   <- nodes(pGraph) %in% deKIDs
    isDiff      <- sum(isDiffExp)
    deSize      <- length(deKIDs)
    allSize     <- length(allKIDs)
    totPath     <- sum(tempNodes)



    fTestRes        <- stats::phyper(isDiff - 1, totPath, allSize - totPath,
                                     deSize, lower.tail = F)
    sampledGraphs   <- pathSampler(pGraph, iter, deKIDs, allKIDs)
    sampleDist      <- stats::ecdf(unlist(sampledGraphs[1]))
    disturbProb     <- 1 - sampleDist(unlist(sampledGraphs[2]))
    causalDist      <- stats::pchisq(-2 * sum(log(fTestRes), log(disturbProb)),
                                     df = 4, lower.tail = FALSE)


    #cat("pathway done: ", Name,"\n")
    return(list(Name,numNodes(pGraph),numEdges(pGraph),fTestRes,
                isDiff, disturbProb, causalDist))

}
