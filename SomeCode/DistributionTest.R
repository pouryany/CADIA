
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


a <- pathSampler.newPath2(ras.path,500000,deKID,allKID,alpha = 0.1, statEval = 1)

samples <- a[[1]]


library(ggplot2)
library(MASS)

fit <- fitdistr(samples,"normal")
?fitdistr
plot(fit)
para <- fit$estimate

hist(samples,prob =T)
curve(dnorm(x, para[1],para[2]), add=TRUE)
median(samples)


install.packages("psych")
library(psych)
describe(samples)

ggplot(data=as.data.frame(samples), aes(samples)) + geom_histogram(bins =500,aes( y = ..density..)) +
    stat_function(linetype="dashed",size =1.5 , fun = dnorm,
                  args = list(mean = mean(samples), sd = sd(samples)),colour = "coral3")+
    scale_x_continuous(name = "Aggregated source/sink score") +
    scale_y_continuous(name = "Probability density")    +
    theme_bw()+
    geom_vline(xintercept = a[[2]], size = 1, colour = "dodgerblue3",linetype = "dashed")+
    theme(
          plot.title = element_text(size = 20),
          axis.title = element_text(size = 20),
          legend.text = element_text(size = 9),
          legend.title=element_text(face = "bold", size = 9),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())

a[[2]]





# Checking the law of  large numbers

i <- 1
ras.path <- pathways.collection[["04014.xml"]]

stat.samples <- list()
while(i < 201){
    i <-100
    deKID2 <- sample(nodes(ras.path),i)
    a <- pathSampler.newPath2(ras.path,10000,deKID2,allKID,alpha = 0.1, statEval = 1)
    samples <- a[[1]]
    hist(samples)
    stat.samples <- rbind(stat.samples, c(i,mean(samples)/i,sd(samples)/i))
    i <- i+1
    print(i)
}

stat.samples <- as.data.frame(stat.samples)


ras.mat <- as(ras.path,"matrix")

centr.mat <- newpath.centrality(ras.mat, 0.1, beta = 1)

cent.val  <- log(rowSums(centr.mat),base = 2)
mean(cent.val/sum(cent.val))
sd(cent.val/sum(cent.val))
hist(cent.val/sum(cent.val),prob=T)


stat.samples[1,] <- c(1,mean(cent.val/sum(cent.val)),sd(cent.val/sum(cent.val))
)
