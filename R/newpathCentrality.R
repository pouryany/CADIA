#' Calculates the Source/Sink Centrality of a pathway
#'
#' @param adj.matrix  the adjacency matrix of a graph to be analyzed
#' @param alpha is the decay factor, similar to the one in Katz centrality
#' @param beta is the relative effect of the sinks vs sources.

#' @return A matrix in which the row i hold the vector of influence of the
#' node i on all the other nodes.
#'
#'
#' @export
newpath.centrality <- function(adj.matrix, alpha, beta){

    # few error throws to be implemented
    # the matrix is not square
    # the alpha is larger then the inverse of the largest eigen value

    eye <- diag(nrow(adj.matrix))
    cent.out  <- solve(eye - alpha * adj.matrix)
    cent.in   <- solve(eye - alpha * t(adj.matrix))
    cent.tot  <- cent.out + beta * cent.in
    return(cent.tot)

}
