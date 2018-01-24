#' Find the total number of paths in a DAG
#'
#' Takes an adjacency matrix input
#' @param adjMatDAG An adjacency matrix of a DAG
#' @return the total number of paths in the graph


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
