#' Extracts all the genes in a given list of pathways
#'
#' @param  listOfPathwayFiles An input containing a list of pathways where the
#' original files are in the format of KGML.
#' @return A list of genes using associated KEGG ID
#' @export
geneKEGGListExtractor <- function(listOfPathwayFiles) {

    keggGeneList <- list()

    for (i in listOfPathwayFiles) {

        tryCatch({
            pathGraphExp <- parseKGML2Graph(i, expandGenes = T)
            keggGeneList <- c(keggGeneList, nodes(pathGraphExp))
        }, error = function(e) {
            print(paste("Pathway Not Found", i))
        })
    }

    keggGeneList <- unique(unlist(keggGeneList))

    return(keggGeneList)
}


