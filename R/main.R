#' Calculates enrichment analysis of KEGG pathways based on causal disturbance
#' model. For more information, see publication doi>10.1145/3107411.3107488
#'
#' @param deIDs  a list of differentially expressed genes
#' @param allIDs list of all genes in the experiment
#' @param iter   Iteration round for sampling causal disturbance
#'
#'
#' @return list containing pathway names, number of pathway nodes and edges,
#'         Over-representation P-values, size of differentially expressed genes
#'         in the pathway, and causal disturbance
#' @importFrom KEGGgraph    translateGeneID2KEGGID
#' @export
#'
#'
causalDisturbance <- function(deIDs, allIDs, iter = 2000, alpha){

    deKID   <- KEGGgraph::translateGeneID2KEGGID(deIDs)
    allKID  <- KEGGgraph::translateGeneID2KEGGID(allIDs)
    len     <- length(pathways.collection.names)
    res     <- vector ("list", length = len)

    for ( i in 1:len ) {
        res[i] <- list(unlist(processPathway(pathways.collection[[i]],pathways.collection.names[[i]],
                                        deKID, allKID, keggRefGenes,iter, alpha )))
        print(cat("Pathway done \n pathway name:", pathways.collection.names[[i]]))
    }

    res <- as.data.frame(res)
    colnames(res) <- names(pathways.collection.names)
    rownames(res)  <- c("Name","nodes","edges","P_ORA",
                        "No. DE","disturbance index", "causal Disturbance")

    return(as.data.frame(t(res)))
}
