pathCounter <- function(adjMatDAG) {
  
  
  #Should add a DAG Checker function
  
  if(length(adjMatDAG) == 1 | is.na(adjMatDAG)){
    
    return(0)
  } else{
    
    
    
    
    rank <- nrow(adjMatDAG)
    totMat = diag(rank)
    
    initial <- diag(rank)
    tot2<- initial - adjMatDAG
    totalPaths <- sum(solve(tot2))
    return(totalPaths)
  }
}


#
#
#
#
#
#
#
#
#
#
#
#




geneKEGGListExtractor <- function(listOfPathwayFiles){
  
  keggGeneList <- list()
  
  for (i in listOfPathwayFiles){
    
   tryCatch({
              pathGraphExp <- parseKGML2Graph(i, expandGenes = T)
              keggGeneList <- c(keggGeneList, nodes(pathGraphExp))
            }, error = function(e) {print(paste("Pathway Not Found", i))})
  }
  
  keggGeneList <- unique(unlist(keggGeneList))
  
  return(keggGeneList)
}





inputKEGGunifier <- function(someGeneList,keggReferencePathwayGenome){
  
  keggUnified <- someGeneList
  keggUnified  <- keggUnified[keggUnified %in% keggReferencePathwayGenome]
  
  return(keggUnified)
  
}






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
    
    #### This will mess up the whole analysis in case CD = 1 because of Statistical assessment 
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








DAGprocessor <- function(inputGraphDAG)
{
  
  if(length(tsort(inputGraphDAG) == 0))
  {
    return(inputGraphDAG)
  } else {
    
    
    cyclesDetector <- RBGL::strongComp(inputGraphDAG)
    compsLen <- lapply(cyclesDetector, length)
    t1 <- which(compsLen > 1)
    
    if(length(t1) == 0){
      # print("Graph is DAG")
      return(inputGraphDAG)
    } else{
      
      for (i in cyclesDetector[t1])
      {
        someLoop <- subGraph(unlist(i), inputGraphDAG)
        nodeHead <- nodes(someLoop)[1]
        nodeTail <- unlist(adj(someLoop,nodeHead))
        inputGraphDAG <- removeEdge(nodeHead , nodeTail, inputGraphDAG)
      }
      return(DAGprocessor(inputGraphDAG))
      
    }
  } 
}







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
  
  fTestRes <- phyper(isDiff -1 , totPath, allSize - totPath, deSize, lower.tail =  F )
  
  
  sampledGraphs <- pathSampler(postDAG, iterationNo, deKIDs, allKIDs)
  
  sampleDist <- ecdf(unlist(sampledGraphs[1]))
  disturbProb <- 1 - sampleDist(unlist(sampledGraphs[2]))
  
  
  causalDisturbance <- pchisq(-2*sum(log(fTestRes) , log(disturbProb)), df = 4 , lower.tail = FALSE)
  
  cat("pathway done: ", pathName)
  return(list(pathName,edgesPre,edgesPos,pathID,fTestRes, isDiff ,totPath,  disturbProb , disturbProb*fTestRes, causalDisturbance))
  
}


processPathwayOnlyORA <- function(pathID, deKIDs, allKIDs , keggRefGenes){
  
  
  
  pathGraphExp <- parseKGML2Graph(pathID, expandGenes = T)
  pathwayFile  <- parseKGML(pathID)
  pathName <- getTitle(pathwayFile)
  
  
  
  # allKID <- inputKEGGunifier(allKIDs , keggRefGenes)
# deKID  <- inputKEGGunifier(deKIDs  , keggRefGenes)
  
  
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
  
  fTestRes <- phyper(isDiff -1 , totPath, allSize - totPath, deSize, lower.tail =  F )
 # fTestRes <- phyper(isDiff -1 , deSize, allSize - totPath, deSize, lower.tail =  F )
  
  
  
  
  print("pathway done")
  return(list(pathName,pathID,fTestRes, isDiff ,totPath))
  
}



CADISA <- function(somethingForNow){
  return(FALSE)
}




