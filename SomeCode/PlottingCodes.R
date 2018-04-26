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

pathways.collection.names


## Plotting some test cases for CADIA

apop.Graph   <- pathways.collection[["04012.xml"]]
node.names   <- nodes(apop.Graph)
node.ENTREZ  <- KEGGgraph::translateKEGGID2GeneID(node.names)
node.genes   <- getSYMBOL(node.ENTREZ,data='org.Hs.eg')
apop.mat     <- as(apop.Graph,"matrix")
cent.mat     <- PathwayDisturbance::newpath.centrality(apop.mat,alpha = 0.1, beta = 1)
cent.vec     <- rowSums(cent.mat)


names(cent.vec) <- node.genes
sort(cent.vec)


#ii      <- cut(cent.vec, breaks = seq(min(cent.vec), max(cent.vec), len =40),
#          include.lowest = TRUE)
ii      <- cut2(cent.vec, m =1 , g = 25)
#ii.sort <- cut(sort(cent.vec), breaks = seq(min(cent.vec), max(cent.vec), len =14),
#               include.lowest = TRUE)
ii.sort <- cut2(sort(cent.vec), m =1 , g = 25)
colors  <- colorRampPalette(c("grey80", "red4"))(25)[ii]
colors2  <- colorRampPalette(c("grey80", "red4"))(25)[ii.sort2]


names(node.genes) <- NULL
apop.network <- network::as.network.matrix(apop.mat)


set.seed(1)
ggnet2(apop.network, node.label = node.genes,node.size = 10,   node.color =  colors,
       arrow.size = 4,label.size = 5,label.trim = T,arrow.gap = 0.015, mode = "fruchtermanreingold") +
    theme(legend.title=element_blank())

?ggnet2


barplot(rep(1,88), col=colors2,yaxt='n', ann=FALSE)

# Comparing centrality scores

ssc.centrality  <- sort(cent.vec,decreasing = T)
ssc.rank        <- rank(ssc.centrality, ties.method = "min")
ssc.tab         <- as.data.frame(ssc.rank)


beet.vec        <- sna::betweenness(apop.network, cmode = "directed")
names(beet.vec) <- node.genes
bet.centrality  <- sort(beet.vec , decreasing = T)
bet.rank        <- rank(bet.centrality, ties.method = "min")
bet.tab         <- as.data.frame(bet.rank)


degr.vec        <- sna::degree(apop.network,cmode = "outdegree")
names(degr.vec) <- node.genes
deg.centrality  <- sort(degr.vec , decreasing = T)
deg.rank        <- rank(deg.centrality, ties.method = "min")
deg.tab         <- as.data.frame(deg.rank)

katz.mat        <- PathwayDisturbance::newpath.centrality(apop.mat,alpha =0.1, beta =0)
katz.vec        <- rowSums(katz.mat)
names(katz.vec) <- node.genes
ktz.centrality  <-sort(katz.vec , decreasing = T)
ktz.rank        <- rank(ktz.centrality, ties.method = "min")
ktz.tab         <- as.data.frame(ktz.rank)


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






#set.seed(13)
#net <-  network(rgraph(10, tprob = 0.15), directed = TRUE)
set.seed(2)
net <-  network(rgraph(14, tprob = 0.1), directed = TRUE,vertex.attrnames = 1:14)
mat.net  <- as.matrix.network(net)
cent.mat <- PathwayDisturbance::newpath.centrality(mat.net,alpha = 0.4, beta = 1)
cent.vec <- rowSums(cent.mat)

ii <- cut(cent.vec, breaks = seq(min(cent.vec), max(cent.vec), len =14),
          include.lowest = TRUE)
ii.sort <- cut(sort(cent.vec), breaks = seq(min(cent.vec), max(cent.vec), len =14),
          include.lowest = TRUE)

colors <- colorRampPalette(c("grey80", "red4"))(14)[ii]

colors.sort <- colorRampPalette(c("grey80", "red4"))(14)[ii.sort]
colors.sort2 <- colorRampPalette(c("grey80", "red4"))(14)[1:14]
#a <- ggnet2(net, arrow.size = 12, arrow.gap = 0.015, node.color =  colors)

set.seed(5)
newpath.color <- ggnet2(net, node.label = 1:14, node.size = 14, arrow.size = 12, arrow.gap = 0.025, node.color =  colors,
       mode = "fruchtermanreingold")


katz.color
katz.rev.color
newpath.color

k  <- katz.color + ggtitle("Katz") +
    theme(plot.title = element_text(size = 30, face = "bold"))

k1 <- katz.rev.color + ggtitle("Katz Reveresed") +
    theme(plot.title = element_text(size = 30, face = "bold"))

k2 <- newpath.color + ggtitle("Source/Sink") +
    theme(plot.title = element_text(size = 30, face = "bold"))


z <- {barplot(1:10, col = colors.sort, axes = F)
mtext("Color ~ Centaliry Score", side=3, adj=0, line=1.2, cex=2, font=2); }
z <- 1:10

library(cowplot)
z <- ggplot(data = as.data.frame(1:88), aes(88:1)) + geom_bar(aes(fill = colors2)) +
    scale_fill_manual(values=colors2) +  ggtitle("Color ~ Centaliry Score")+
    theme(legend.position="none", panel.grid.major = element_blank()
          ,panel.grid.minor = element_blank(), axis.line = element_blank(),
          axis.text = element_blank(),panel.background = element_blank(),
          axis.title  = element_blank(), axis.ticks = element_blank(),
          plot.title = element_text(size = 10, face = "bold"))


library(RColorBrewer)
display.pal(colors.sort)

library(gridExtra)

grid.arrange(k,k1,k2, z, nrow =2)
