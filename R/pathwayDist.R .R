#'  Calculates some information regarding the pathways
#'
#'  @param pathID KEGG ID of an input pathway
#'  @param deKIDs lsit of differentially expressed genes
#'  @param allKIDs list of all genes with their associated KEGG IDs
#'  @param keggRefGenes list of reference genes from KEGG
#'  @param iternationNo the number of iterations for causal disturbance
#'
#'
#'
#'

processPathway <- function(pathID, deKIDs, allKIDs , keggRefGenes, iterationNo){



    pathGraphExp <- parseKGML2Graph(pathID, expandGenes = T)
    pathwayFile  <- parseKGML(pathID)




    postDAG  <- DAGprocessor(pathGraphExp)
    pathName <- getTitle(pathwayFile)
    edgesPre <- numEdges(pathGraphExp)
    edgesPos <- numEdges(postDAG)






    allKID <- inputKEGGunifier(allKIDs , keggRefGenes)
    deKID  <- inputKEGGunifier(deKIDs  , keggRefGenes)


    tempNodes <- nodes(pathGraphExp) %in% allKIDs
    nodeNums <- nodes(pathGraphExp)[tempNodes]
    isDiffExp <- nodes(pathGraphExp) %in% deKIDs



    isDiff  <- sum(isDiffExp)
    deSize  <- length(deKIDs)
    allSize <- length(allKIDs)
    totPath <- sum(tempNodes)



    # fTestMat <-  matrix(c(isDiff, deSize - isDiff, totPath - isDiff,
    #                        allSize - deSize - totPath + isDiff),
    #                        nrow = 2,
    #                         dimnames =
    #                             list(c("Pathway", "NotPathway"),
    #                                    c("DE", "NotDE")))
    #


    #fTestRes <- fisher.test(fTestMat, alternative = "greater")

    fTestRes <- phyper(isDiff -1 , totPath, allSize - totPath,..
                       deSize, lower.tail =  F )


    sampledGraphs <- pathSampler(postDAG, iterationNo, deKIDs, allKIDs)

    sampleDist <- ecdf(unlist(sampledGraphs[1]))
    disturbProb <- 1 - sampleDist(unlist(sampledGraphs[2]))


    causalDisturbance <- pchisq(-2*sum(log(fTestRes) , log(disturbProb)),..
                                df = 4 , lower.tail = FALSE)

    cat("pathway done: ", pathName)
    return(list(pathName,edgesPre,edgesPos,pathID,fTestRes, isDiff ,totPath,..
                disturbProb , disturbProb*fTestRes, causalDisturbance))

}
