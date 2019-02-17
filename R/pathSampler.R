#' Calculates the distribution of pathway disturbance
#'
#' @param inputGraph A graphNEL object of a pathway
#' @param interationNo number of rounds of sampling
#' @param deKID The name of the nodes in the subset of interest. e.g. DEG
#' @param iter the number of iterations for causal disturbance
#' @param alpha the dampening factor for Source/Sink Centrality
#' @param beta the relative Source vs Sink factor
#' @param statEval Choose 1 for product-based, 0 for summation-based
#'
#' @return the probability of observing a subset of
#'
#'
#' @export
pathSampler <-  function(inputGraph,deKID,iterationNo,alpha,
                         beta, statEval =1 ) {


                totPathNodes <- nodes(inputGraph)
                totGenes     <- length(totPathNodes)
                sizeDE       <- sum(totPathNodes %in% deKID)
                samplingData <- rep(0,iterationNo)
                pathMat      <- as(inputGraph, "matrix")
                deGenesInd   <- totPathNodes %in% deKID
                deGenes      <- totPathNodes[deGenesInd]
                deMatUnRef   <- pathMat[deGenesInd,deGenesInd]



                if (statEval == 0) {

                    if (length(deMatUnRef) < 2)
                        {
                        causalDisturbance <- 1
                        return(1)
                    } else if (sizeDE == 0){
                        return(1)
                    }else{


                          centr.mat <- source.sink.centrality(pathMat, alpha,
                                                              beta)
                          paths.tot <- sum(centr.mat)
                          cdist.tot <- sum(centr.mat[deGenesInd,])
                          causalDisturbance <- cdist.tot/paths.tot

                        for (i in 1:iterationNo) {
                            randPerm <- logical(totGenes)
                            posPerm  <- sample(1:totGenes, sizeDE,replace = F)
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


                        centr.mat <- source.sink.centrality(pathMat, alpha,
                                                            beta)
                        paths.tot <- rowSums(centr.mat)
                        cdist.tot <- rowSums(centr.mat[deGenesInd,])
                        paths.log <- sum(log2(paths.tot))
                        cdist.log <- sum(log2(cdist.tot))

                        causalDisturbance <- cdist.log/paths.log

                        for (i in 1:iterationNo) {
                            randPerm <- logical(totGenes)
                            posPerm  <- sample(1:totGenes, sizeDE,replace = F)
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

