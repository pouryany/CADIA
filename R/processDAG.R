#'  Finds an underlying DAG of a given graph pathway
#'
#'  @param inputGraphDAG is the graph of the pathways
#'  @return processed graph without anyloops or self loops
#'
#'
#'
#'
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

