
# The lines 4- 50 are generated from NCBI GEO portal
# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Feb 20 21:18:14 EST 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)
library(KEGGgraph)
# load series and platform data from GEO
# 30 ovarian LMP vs 60 ovarian cancer
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
sml     <- paste("G", sml, sep="")    # set group names
fl      <- as.factor(sml)
gset$description <- fl
design  <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit     <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2    <- contrasts.fit(fit, cont.matrix)
fit2    <- eBayes(fit2, 0.01)
tT      <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)



tT.filter  <- tT[!is.na(tT$Gene.ID),]
tT.filter  <- tT.filter[!duplicated(tT.filter$Gene.ID),]
tT.deGenes <- tT.filter[tT.filter$adj.P.Val < 0.05, ]
tT.deGenes <- tT.deGenes[abs(tT.deGenes$logFC) >1,]

tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)

deKID    <- translateGeneID2KEGGID(tT.de.names)
allKID   <- translateGeneID2KEGGID(tT.all.names)

pathSampler.newPath2 <- function(pathwayRefGraphExpanded,iterationNo, deKID,
                                allKID,alpha, statEval) {


    totPathNodes     <- nodes(pathwayRefGraphExpanded)
    totGenes         <- length(totPathNodes)
    eyeTot           <- diag(totGenes)

    sizeDE       <- sum(totPathNodes %in% deKID)
    eyeDE        <- diag(totGenes - sizeDE)
    samplingData <- rep(0,iterationNo)
    pathMat      <- as(pathwayRefGraphExpanded, "matrix")
    #totalPaths   <- pathCounter.(pathMat,eyeTot,0.5)
    deGenesInd   <- totPathNodes %in% deKID
    deGenes      <- totPathNodes[deGenesInd]
    eye          <- diag(nrow(pathMat))
    deMatUnRef   <- pathMat[deGenesInd,deGenesInd]



    if (statEval == 0) {

        if (length(deMatUnRef) < 2)
        {
            causalDisturbance <- 1
            return(1)
        } else if (sizeDE == 0){
            return(1)
        }else{


            centr.mat <- newpath.centrality(pathMat, alpha, beta = 1)
            paths.tot <- sum(centr.mat)
            cdist.tot <- sum(centr.mat[deGenesInd,])
            causalDisturbance <- cdist.tot/paths.tot

            for (i in 1:iterationNo) {
                randPerm     <- logical(totGenes)
                posPerm      <- sample(1:totGenes, sizeDE,replace = F)
                randPerm[posPerm] = TRUE


                cdist.tot.rand <- sum(centr.mat[randPerm,])
                samplingData[i] <- (cdist.tot.rand / paths.tot)
            }

            sampleDist      <- stats::ecdf(unlist(samplingData))
            disturbProb     <- 1 - sampleDist(causalDisturbance)
            iter2 <- iterationNo

            if(disturbProb == 0)  {
                disturbProb <- NA
            }



            return(disturbProb)
            # return(list(samplingData, causalDisturbance))
        }
    }
    else if (statEval == 1) {

        if (length(deMatUnRef) < 2)
        {
            causalDisturbance <- 1
            return(1)
        } else if (sizeDE == 0){
            return(1)
        }else{

            tryCatch({


                centr.mat <- newpath.centrality(pathMat, alpha, beta = 1)
                paths.tot <- rowSums(centr.mat)
                cdist.tot <- rowSums(centr.mat[deGenesInd,])
                paths.log <- sum(log2(paths.tot))
                cdist.log <- sum(log2(cdist.tot))

                causalDisturbance <- cdist.log/paths.log

                for (i in 1:iterationNo) {
                    randPerm     <- logical(totGenes)
                    posPerm      <- sample(1:totGenes, sizeDE,replace = F)
                    randPerm[posPerm] = TRUE


                    cdist.tot.rand  <- rowSums(centr.mat[randPerm,])
                    cdist.tot.rand  <- sum(log2(cdist.tot.rand))

                    samplingData[i] <- (cdist.tot.rand /paths.log)
                }

                sampleDist      <- stats::ecdf(unlist(samplingData))
                disturbProb     <- 1 - sampleDist(causalDisturbance)
                iter2 <- iterationNo

                if(disturbProb == 0)  {
                    disturbProb <- 1/iterationNo
                }



                #return(disturbProb)
                return(list(samplingData, causalDisturbance))
            }, error = function(e) {
                disturbProb <- NA
                return(disturbProb)
            })
        }
    }

}



ras.path <- pathways.collection[["04510.xml"]]


a <- pathSampler.newPath2(ras.path,500000,deKID,allKID,alpha = 0.1,statEval = 1)

samples <- a[[1]]


library(ggplot2)
library(MASS)

fit  <- fitdistr(samples,"normal")
para <- fit$estimate


ggplot(data=as.data.frame(samples), aes(samples)) +
    geom_histogram(bins =500,aes( y = ..density..)) +
    stat_function(linetype="dashed",size =1.5 , fun = dnorm,
                  args = list(mean = mean(samples), sd = sd(samples)),
                  colour = "coral3")+
    scale_x_continuous(name = "Aggregated source/sink score") +
    scale_y_continuous(name = "Probability density")    +
    theme_bw()+
    geom_vline(xintercept = a[[2]], size = 1, colour = "dodgerblue3",
               linetype = "dashed")+
    theme(
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 9),
          legend.title=element_text(face = "bold", size = 9),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())


