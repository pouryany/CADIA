# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Feb 20 21:18:14 EST 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO
# 30 ovaria LMP vs 60 ovarian cancer
gset <- getGEO("GSE12172", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("11111010110111001100111111110101001010101111101111",
               "1101111010011110011110000011111010010111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)



tT.filter  <- tT[!is.na(tT$Gene.ID),]
tT.filter  <- tT.filter[!duplicated(tT.filter$Gene.ID),]
tT.deGenes <- tT.filter[tT.filter$adj.P.Val < 0.05, ]
tT.deGenes <- tT.deGenes[abs(tT.deGenes$logFC) >1,]
tT.deGenes

tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)
deKID    <- translateGeneID2KEGGID(tT.de.names)
allKID   <- translateGeneID2KEGGID(tT.all.names)

library(stringr)
tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 2000,
                                 alpha = 0.1 , statEval = 1)

names(tT.pathways)
tT.pathways.clean <- tT.pathways[!is.na(tT.pathways$`disturbance index`),]
tT.pathways.clean$CDIST  <- p.adjust(as.numeric(as.character(
    tT.pathways.clean$`causal Disturbance`))
    ,method = "BH")
tT.pathways.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                        (tT.pathways.clean$P_ORA)),method = "fdr")


tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,]
tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]

head(tT.pathways.clean[order(tT.pathways.clean$CDIST),],20)


library(PathwayDisturbance)




#exporting results

tT.pathways.clean$KEGGID <- str_sub(rownames(tT.pathways.clean), end = -5)

rownames(tT.pathways.clean) <- NULL

Hodgkins.cdist  <- tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,]
Hodgkins.ora    <- tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]
sapply(Hodgkins.cdist, mode)


Hodgkins.cdist <- Hodgkins.cdist[order(Hodgkins.cdist$CDIST),c(1,10,4,6,8,9)]
Hodgkins.ora <- Hodgkins.ora[order(Hodgkins.ora$ORAFDR),c(1,10,8,9)]
Hodgkins.cdist[,c(3,4,5,6)] <- mapply(as.character,Hodgkins.cdist[,c(3,4,5,6)])
Hodgkins.cdist[,c(3,4,5,6)] <- mapply(as.numeric,Hodgkins.cdist[,c(3,4,5,6)])

Hodgkins.cdist[,c(3,4,5,6)] <- mapply(formatC,Hodgkins.cdist[,c(3,4,5,6)],
                                      MoreArgs = list(format = "e", digits = 2))
Hodgkins.ora[,c(3,4)]   <- mapply(formatC,Hodgkins.ora[,c(3,4)],
                                  MoreArgs = list(format = "e", digits = 2))

Hodgkins.ora
Hodgkins.cdist


library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

print(xtable(Hodgkins.cdist), include.rownames = FALSE)
print(xtable(Hodgkins.ora), include.rownames = FALSE)




#SPIA RESULTs


library(SPIA)
deSPIA <- tT.deGenes$logFC
names(deSPIA)<- tT.deGenes$Gene.ID
allSPIA <- tT.filter$Gene.ID


resSPIA <- spia(de=deSPIA,all=allSPIA,organism="hsa",nB=2000
                ,plots=FALSE,beta=NULL,combine="fisher",verbose=FALSE)


head(resSPIA)
resSPIA.report <- resSPIA[order(resSPIA$pGFdr),c(1,2,9)]
resSPIA.report[,3] <- mapply(formatC,resSPIA.report[,3],
                             MoreArgs = list(format = "e", digits = 2))

resSPIA.report <- resSPIA.report[as.numeric(resSPIA.report$pGFdr) <0.05,]
print(xtable(resSPIA.report), include.rownames = FALSE)









library(stringr)
library(gage)
data("kegg.gs")
kg.hsa=kegg.gsets("hsa")
kegg.gs=kg.hsa$kg.sets[kg.hsa$sigmet.idx]
av.paths <- names(kegg.gs)
av.paths <- str_sub(av.paths, start =10)
pathways.collection.names[!(pathways.collection.names %in% av.paths)]
p.kegg.gsets <- lapply(pathways.collection,nodes)
p.kegg.gsets <- mapply(str_sub,p.kegg.gsets, MoreArgs = list(start = 5) )

names(p.kegg.gsets) <- pathways.collection.names
p.kegg.gsets
expdata <- exprs(gset)
head(tT.filter)
expdata.clean <- expdata[row.names(expdata) %in% rownames(tT.filter),]
head(expdata.clean)

tmp1 <- tT.filter[rownames(expdata.clean),]
rownames(expdata.clean) <- tmp1$Gene.ID
rownames(expdata.clean) <-
    nsample <- grep("G0", fl)
csample <- grep("G1", fl)


a <- gage(expdata.clean, gsets = p.kegg.gsets, ref = nsample, sample = csample,
          compare = "unpaired",saaTest = gs.KSTest)

head(a$greater[,1:5],20)







tT.filter$Gene.ID

log.fc <- tT.filter$logFC
class(as.vector(log.fc))
names(log.fc) <- tT.filter$Gene.ID
class(log.fc)

a <- gage(log.fc, gsets = p.kegg.gsets, ref = NULL, sample = NULL)
tT.filter

head(a$greater[,1:5],40)















#Random Testing

??data_frame

library(dplyr)
err.samples.ora <- data_frame()
err.samples.cdist <- data_frame()
for(i in 1:50){


    for(j in 1:10){
        tT.de.names   <- sample(tT.all.names,i*100)


        tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 2000, alpha = 0.1,statEval = 1)
        tT.pathways.clean<- tT.pathways #[tT.pathways$`disturbance index` !=0,]
        tT.pathways.clean$CDIST  <- p.adjust(as.numeric(as.character(
            tT.pathways.clean$`causal Disturbance`))
            ,method = "fdr")
        tT.pathways.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                                        (tT.pathways.clean$P_ORA)),method = "fdr")

        #hist(as.numeric(as.character(tT.pathways.clean$`causal Disturbance`)))

        err.samples.cdist[j,i] <- nrow(tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,])
        err.samples.ora[j,i]   <- nrow(tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,])
    }

}

err.cdist <- t(err.samples.cdist)
err.cdist <- as.data.frame(err.cdist)
err.cdist$av.err <- rowMeans(err.cdist) /148
tail(t(err.samples.ora))

err.ora <- t(err.samples.ora)
err.ora <- as.data.frame(err.ora)
err.ora$av.err <- rowMeans(err.ora) / 148

plot(err.ora$av.err, col ="#e41a1c",cex =1.9,pch = "o",
     xlab="Size of sample set (x100)" , ylab="ratio of false positives",cex.lab=1.8)
#par(new=TRUE)
points(err.cdist$av.err,cex =1.5,col ="#00441b",pch = 15)
legend(x = 0, y = 0.002, legend = c("ORA", "CADIA"),
       col= c("#e41a1c","#00441b"), cex = 2, pch = c(15,15))

d <-cbind(err.cdist$av.err,err.ora$av.err)
d <- as.data.frame(d)
names(d) <- c("CADIA","ORA")
d$iteration  <- 1:50 *100
library(ggplot2)
library(reshape)
dm <- melt(d,id.vars = 3)
ggplot(data = dm, aes(x = iteration, y = value, shape = variable,color = variable))+
    geom_point(size = 5, stroke = 1.5)+
    scale_shape( solid = TRUE)+
    scale_color_manual(values = c("CADIA" = 'red4','ORA' = '#018571')) +
    scale_shape_manual(values = c("CADIA" = 5,'ORA' = 4)) +
    geom_hline(yintercept=0.05, size = 2, linetype="dashed", color = "red")+
    theme_bw()+
    labs(y = "Average false positives")+
    theme(
        axis.title = element_text(size = 30),
        legend.text = element_text(size = 25),
        axis.text=element_text(size= 30),
        axis.text.y= element_text(size = 30),
        legend.title=element_blank())

write.csv(err.samples.cdist, file = "err_cdist")
write.csv(err.samples.cdist, file = "err_ora")

head(tT.pathways.clean[order(tT.pathways.clean$ORAFDR),],20)

gg  <- pathways.collection[["04144.xml"]]
cgg <- connectedComp(gg)
subgg <- unique(unlist(cgg[which((lapply(cgg, length)) != 1)]))
sgg <- subGraph(subgg,gg)
# Testing the stability of  newpath Centrality

testGraph <- pathways.collection[["05219.xml"]]
testGraph.adj <- as(testGraph, "matrix")




cent.matrix <- PathwayDisturbance::newpath.centrality(testGraph.adj,alpha = 0.5, beta = 0.5)








causalDisturbance2 <- function(deIDs, allIDs, iter = 2000, alpha = 0.1,
                              beta =1,statEval = 1 , fdrMethod = "BH"){


    len     <- length(pathways.collection.names)
    res     <- vector ("list", length = len)

    for ( i in 1:len ) {
        res[i] <- list(unlist(processPathway(pathways.collection[[i]],
                                             pathways.collection.names[[i]],
                                             deIDs, allIDs,iter,
                                             alpha,statEval )))
        print(cat("Pathway done \n pathway name:", pathways.collection.names[[i]]))
    }

    res <- as.data.frame(res)
    colnames(res) <- names(pathways.collection.names)
    rownames(res)  <- c("Name","nodes","edges","P_ORA",
                        "No. DE","disturbance index", "causal Disturbance")

    res <- as.data.frame(t(res))

    res.clean <- res[!is.na(res$`disturbance index`),]
    res.clean$cadia  <- p.adjust(as.numeric(as.character(
        res.clean$`causal Disturbance`)),method = fdrMethod)

    res.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                            (res.clean$P_ORA)),method = fdrMethod)

    res.clean$KEGGID <- str_sub(rownames(res.clean), end = -5)
    rownames(res.clean) <- NULL

    res.clean <- res.clean[order(res.clean$cadia),]

    res.clean[,c(4,6,7,8,9)] <- mapply(as.character,res.clean[,c(4,6,7,8,9)])
    res.clean[,c(4,6,7,8,9)] <- mapply(as.numeric,res.clean[,c(4,6,7,8,9)])
    res.clean[,c(4,6,7,8,9)] <- mapply(formatC,res.clean[,c(4,6,7,8,9)],
                                       MoreArgs = list(format = "e", digits = 3))


    return(res.clean)
}




