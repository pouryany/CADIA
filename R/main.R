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
causalDisturbance <- function(deIDs, allIDs, iter = 2000){

    deKID   <- KEGGgraph::translateGeneID2KEGGID(deIDs)
    allKID  <- KEGGgraph::translateGeneID2KEGGID(allIDs)
    len     <- nrow(cleanPathNames)
    res     <- vector ("list", length = len)

    for ( i in 1:len ) {
        res[i] <- list(unlist(processPathway(cleanPathList[[i]],cleanPathNames[i,1],
                                        deKID, allKID, keggRefGenes,iter )))
        print(cat("Pathway done \n pathway name:", cleanPathNames[i,1]))
    }

    res <- as.data.frame(res)
    colnames(res) <- cleanPathNames[,1]
    rownames(res)  <- c("Name","nodes","edges","P_ORA",
                        "No. DE","disturbance index", "causal Disturbance")

    return(as.data.frame(t(res)))
}
