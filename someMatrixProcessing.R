library(RBGL)

mtx.collection <-sapply(pathways.collection, function(X)(as(X,"matrix")))

eigen.collection <- lapply(mtx.collection, eigen, only.value = T )

sapply(eigen.collection,function (X)(max(X["values"])))
a <- eigen.collection[[1]]

b <- unlist(a[[148]])

1/ max(Re(b[Im(unlist(b)) == 0]))


library(dplyr)

largest.eigen <- function(b) {
c <-  b %>% unlist()  %>% Im() == 0
c %>% subset.default(x = b,) %>% Re() %>% max()
}

max(sapply(eigen.collection,function(X)(largest.eigen(unlist(X)))))

bugbug <- unlist(pathways.collection[["04070.xml"]])
bugbug
library(RBGL)
bugmat <- as(bugbug,"matrix")

newpath.centrality(bugmat,alpha = 100, beta = 1)

tryCatch({
    solve(bugmat)
}, error = function(e) {
   "singular"
})


 sapply(mtx.collection, function(X) (tryCatch({
    solve(X)
}, error = function(e) {
    "singular"
})
))

