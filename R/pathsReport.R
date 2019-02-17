#' Generates a summary of the pathways used in CADIA
#'
#'
#' @return The list of pathways used in CADIA and their KEGG IDs
#'
#'
#' @export
cadia.paths <- function() {

    KEGG_IDs <- stringr::str_sub(names(pathways.collection.names), end = -5)
    df <- data.frame(KEGG_IDs,pathways.collection.names,
                     stringsAsFactors = F)
    row.names(df) <- NULL
    return(df)

}

