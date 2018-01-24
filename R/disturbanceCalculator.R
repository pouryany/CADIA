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

    # Pathway Reference to Graph and Matrix
    gene2pathNodeRef <- list()

    # for(i in nodes(pathwayRefGraphUnexpanded)){
    #   tempNode <- getKEGGnodeData(pathwayRefGraphUnexpanded, i)
    #
    #   for (j in getName(tempNode)){
    #     gene2pathNodeRef[[j]]  <- i
    #   }
    # }
    #
    #Should put in graph derefrences

    totGenes <- length(nodes(pathwayRefGraphExpanded))
    # totMat = diag(totGenes)


    sizeDE <- sum(nodes(pathwayRefGraphExpanded) %in% deKID)
    # the above number is the sizeDE from



    samplingData <- rep(0,iterationNo)

    pathwayRefDAG <- DAGprocessor(pathwayRefGraphExpanded)
    pathwayRefDAG <- removeSelfLoops(pathwayRefDAG)
    pathMat <- as(pathwayRefDAG, "matrix")
    totalPaths <- pathCounter(pathMat)


    deGenesInd <- nodes(pathwayRefGraphExpanded) %in% deKID
    deGenes <- nodes(pathwayRefGraphExpanded)[deGenesInd]



    # deInPathNames<- unique(unlist(  deInPath, use.names = FALSE ))
    # deInPathRevise <- nodes(pathwayRefGraphUnexpanded) %in% deInPathNames
    #
    deMatUnRef <- pathMat[!deGenesInd,!deGenesInd]

    if (length(deMatUnRef) < 2){

        #### This will mess up the whole analysis in case CD = 1 because
        #of Statistical assessment
        causalDisturbance <- 1
    } else{



        deTotalPathsUnRef <- pathCounter(deMatUnRef)


        causalDisturbance <- 1 - (deTotalPathsUnRef/totalPaths)

    }



    for (i in 1:iterationNo) {
        randPerm <- logical(totGenes)
        posPerm <- sample(1:totGenes, sizeDE,replace = F)
        randPerm[posPerm] = TRUE

        # inPathSample <- nodes(pathwayRefGraphExpanded)[randPerm]


        randMatUnRefSample <- pathMat[!randPerm,!randPerm]

        if(length(randMatUnRefSample) < 2){
            samplingData[i] <- 1
        } else{
            totalPathsUnRefSample <- pathCounter(randMatUnRefSample)


            samplingData[i] <- 1 - (totalPathsUnRefSample / totalPaths)
        }}


    return(list(samplingData, causalDisturbance))

}
