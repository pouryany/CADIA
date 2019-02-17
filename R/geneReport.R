#' Returns the list of DE genes in a set of pathways
#'
#' @param deIDs  a list of differentially expressed genes
#' @param pathwayIDs a list KEGG IDs of pathways for which we want the DEGs
#'
#' @importFrom magrittr %>%
#' @importFrom KEGGgraph    translateGeneID2KEGGID
#' @importFrom KEGGgraph    translateKEGGID2GeneID
#' @importFrom graph    nodes
#'
#' @export
#'
#'
geneReport <- function(deIDs, pathwayList = NULL){


    if(is.null(pathwayList)) {
        stop("No Pathway List provided")
    }else{
          tryCatch({
            pathList <- paste0(pathwayList,".xml")

            a        <- CADIA::pathways.collection[pathList]
            deKIDs   <- KEGGgraph::translateGeneID2KEGGID(deIDs)

            genes <- sapply(a, function(X){
                isDiffExp   <- nodes(X) %in% deKIDs
                isDiffExp   <- nodes(X)[isDiffExp]
                x <- KEGGgraph::translateKEGGID2GeneID(isDiffExp)
                q <- stringi::stri_join_list(list(x), sep = "/",
                                             collapse = NULL)
                return(q)
            })

            path.names <- CADIA::pathways.collection.names[pathList]
            return(data.frame(path.names,genes))
            }, error = function(e){
            stop("The pathway list may contain invalid IDs or pathways without
                     differentially expressed nodes")})
          }
}
