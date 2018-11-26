#Random Graph testing

install.packages("GGally")
install.packages("network")
install.packages("sna")
install.packages("Hmisc")


library(GGally)
library(network)
library(sna)
library(ggplot2)
library(RBGL)
library(KEGGgraph)
library(org.Hs.eg.db)
library(annotate)
library(Hmisc)
library(CADIA)


## Plotting some test cases for CADIA

apop.Graph   <- pathways.collection[["04012.xml"]]
node.names   <- nodes(apop.Graph)
node.ENTREZ  <- KEGGgraph::translateKEGGID2GeneID(node.names)
node.genes   <- getSYMBOL(node.ENTREZ,data='org.Hs.eg')
apop.mat     <- as(apop.Graph,"matrix")
cent.mat     <- PathwayDisturbance::newpath.centrality(apop.mat,alpha = 0.1,
                                                       beta = 1)
cent.vec     <- rowSums(cent.mat)



# Testing for beta parameter
#
#
test.mat   <- matrix(NA, length(cent.vec),10000)
beta.array <- seq(0,10,by = 0.001)



for (i in 1:10000) {
    mat           <- PathwayDisturbance::newpath.centrality(apop.mat,
                                                      alpha = 0.1,
                                                      beta = beta.array[i])
    vec           <- rowSums(mat)
    test.mat[,i]  <- vec/(1+beta.array[i])

}


sink.mat     <- PathwayDisturbance::newpath.centrality(t(apop.mat),
                                                  alpha = 0.1,
                                                  beta = 0)
sink.vec     <- rowSums(sink.mat)
sink.rank    <- rank(sink.vec, ties.method = "min")
ssc.mat      <- PathwayDisturbance::newpath.centrality(t(apop.mat),
                                                       alpha = 0.1,
                                                       beta = 1)
ssc.vec      <- rowSums(ssc.mat)
ssc.rank     <- rank(ssc.vec, ties.method = "min")


row.names(test.mat) <- names(cent.vec)



dif.array4 <- rep(NA,10000)

for (i in 2:10000){
   dif.array4[i] <-  sqrt(sum((test.mat[,i] - sink.vec) **2))  *
                     sqrt(sum((test.mat[,i] - test.mat[,1]) **2))
   }


plot(beta.array[4:10000],dif.array4[4:10000],type = "l", xlab = "Beta",
     ylab = "product of L-2 norms of Source and Sink",cex.lab=1.5)
abline(v = 1, col="red", lwd=3, lty=2)





# Plotting the ErbB graph

ii      <- cut2(cent.vec, m =1 , g = 25)
ii.sort <- cut2(sort(cent.vec), m =1 , g = 25)
colors  <- colorRampPalette(c("grey80", "red4"))(25)[ii]


names(node.genes) <- NULL
apop.network <- network::as.network.matrix(apop.mat)


set.seed(1)
ggnet2(apop.network, node.label = node.genes,node.size = 10,   node.color =  colors,
       arrow.size = 4,label.size = 5,label.trim = T,arrow.gap = 0.015,
       mode = "fruchtermanreingold") +
       theme(legend.title=element_blank())




# Comparing centrality scores across different methods


# SSC with beta =1
names(cent.vec) <- node.genes
ssc.centrality  <- sort(cent.vec,decreasing = T)
ssc.rank        <- rank(ssc.centrality, ties.method = "min")
ssc.tab         <- as.data.frame(ssc.rank)


# Betweenness Centrality
beet.vec        <- sna::betweenness(apop.network, cmode = "directed")
names(beet.vec) <- node.genes
bet.centrality  <- sort(beet.vec , decreasing = T)
bet.rank        <- rank(bet.centrality, ties.method = "min")
bet.tab         <- as.data.frame(bet.rank)


# Degree Centrality
degr.vec        <- sna::degree(apop.network,cmode = "outdegree")
names(degr.vec) <- node.genes
deg.centrality  <- sort(degr.vec , decreasing = T)
deg.rank        <- rank(deg.centrality, ties.method = "min")
deg.tab         <- as.data.frame(deg.rank)


# Katz Centrality
katz.mat        <- PathwayDisturbance::newpath.centrality(apop.mat,alpha =0.1,
                                                          beta =0)
katz.vec        <- rowSums(katz.mat)
names(katz.vec) <- node.genes
ktz.centrality  <-sort(katz.vec , decreasing = T)
ktz.rank        <- rank(ktz.centrality, ties.method = "min")
ktz.tab         <- as.data.frame(ktz.rank)


# Creating a merged list
a <- merge(ssc.tab,bet.tab, by = 0)
b <- merge(deg.tab,ktz.tab, by = 0)

c <- merge(a,b, by =1)
colnames(c) <- c("Gene", "SSC", "Bet", "Deg","Katz")

d <- c[order(c$SSC,decreasing = T),]
d1 <- d[1:30,]
d2 <- d[31:60,]
d3 <- d[61:88,]
d3 <- rbind(d3, c("","","",""))
d3 <- rbind(d3, c("","","",""))
e  <- cbind(d1,d2,d3)


library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

print(xtable(e), include.rownames = FALSE)



# Plotting some demonstration graphs

set.seed(2)
net      <-  network(rgraph(14, tprob = 0.1), directed = TRUE,
                     vertex.attrnames = 1:14)
mat.net  <- as.matrix.network(net)
cent.mat <- PathwayDisturbance::newpath.centrality(mat.net,alpha = 0.4,
                                                   beta = 1)
cent.vec <- rowSums(cent.mat)

ii       <- cut(cent.vec, breaks = seq(min(cent.vec), max(cent.vec), len =14),
                include.lowest = TRUE)
ii.sort  <- cut(sort(cent.vec), breaks = seq(min(cent.vec), max(cent.vec),
                                             len =14),include.lowest = TRUE)
colors   <- colorRampPalette(c("grey80", "red4"))(14)[ii]


colors.sort  <- colorRampPalette(c("grey80", "red4"))(14)[ii.sort]
colors.sort2 <- colorRampPalette(c("grey80", "red4"))(14)[1:14]

set.seed(5)
newpath.color <- ggnet2(net, node.label = 1:14, node.size = 14, arrow.size = 12,
                        arrow.gap = 0.025, node.color =  colors,
       mode = "fruchtermanreingold")

