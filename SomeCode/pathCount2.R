#' Find the total number of paths in a DAG
#'
#' Takes an adjacency matrix input
#' @param adjMatDAG An adjacency matrix of a DAG
#' @param eye  is the identity function
#' @return the total number of paths in the graph


pathCounter <- function(adjMatDAG, eye) {

    if (length(adjMatDAG) <= 1) {
        return(0)
    } else {
        totalPaths <- sum(solve(eye - adjMatDAG))
        return(totalPaths)
    }
}
