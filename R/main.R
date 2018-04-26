#' Calculates enrichment analysis of KEGG pathways based on causal disturbance
#' model. For more information, see publication doi>10.1145/3107411.3107488
#'
#' @param deIDs  a list of differentially expressed genes
#' @param allIDs list of all genes in the experiment
#' @param iter   Iteration round for sampling causal disturbance
#' @param alpha  Decay factor for Source/Sink centrality, default 0.1.
#' @param beta   Relative importance of sink versus source
#' @param statEval 0 for summation of source/sink centrality 1 for product
#' @param fdrMethod choice of false disovery rate method. inherited from p.adjust
#'
#' @return list containing pathway names, number of pathway nodes and edges,
#'         Over-representation P-values, size of differentially expressed genes
#'         in the pathway, and causal disturbance
#'
#' @export
#'
#'
causalDisturbance <- function(deIDs, allIDs, iter = 2000, alpha = 0.1,
                              beta =1,statEval = 1 , fdrMethod = "BH"){


    len     <- length(pathways.collection.names)
    res     <- vector ("list", length = len)

    for ( i in 1:len ) {
        res[i] <- list(unlist(processPathway(pathways.collection[[i]],
                                             pathways.collection.names[[i]],
                                             deIDs, allIDs,iter,
                                             alpha,statEval )))
        print(cat("Pathway done \n pathway name:", pathways.collection.names[[i]]))
    }

    res <- as.data.frame(res)
    colnames(res) <- names(pathways.collection.names)
    rownames(res)  <- c("Name","nodes","edges","P_ORA",
                        "No. DE","disturbance index", "causal Disturbance")

    res <- as.data.frame(t(res))

    res.clean <- res[!is.na(res$`disturbance index`),]
    res.clean$cadia  <- p.adjust(as.numeric(as.character(
        res.clean$`causal Disturbance`)),method = fdrMethod)

    res.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                (res.clean$P_ORA)),method = fdrMethod)

    res.clean$KEGGID <- str_sub(rownames(res.clean), end = -5)
    rownames(res.clean) <- NULL

    res.clean <- res.clean[order(res.clean$cadia),]

    res.clean[,c(4,6,7,8,9)] <- mapply(as.character,res.clean[,c(4,6,7,8,9)])
    res.clean[,c(4,6,7,8,9)] <- mapply(as.numeric,res.clean[,c(4,6,7,8,9)])
    res.clean[,c(4,6,7,8,9)] <- mapply(formatC,res.clean[,c(4,6,7,8,9)],
                                    MoreArgs = list(format = "e", digits = 3))


    return(res.clean)
}
