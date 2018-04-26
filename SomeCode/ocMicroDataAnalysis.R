# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Tue Feb 13 12:31:25 EST 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE32867", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL6884", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("10101010101010101010101010101010101010101010101010",
               "10101010101010101010101010101010101010101010100101",
               "0101010110101010")
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

#boxplot(exprs(gset))

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
gset$batch <- gset$`batch:ch1`
gset$sid <- gset$`subject:ch1`
gset$batch

design <- model.matrix(~ description +gset$sid + gset$batch , gset)
design
fit <- lmFit(gset, design)
fit1 <- eBayes(fit)
summary(decideTests(fit1))

tT <- topTable(fit1, coef = "descriptionG1", adjust.method = "fdr"
               , number = Inf,sort.by = "B")
tT.filter  <- tT[!is.na(tT$Gene.ID),]
tT.filter  <- tT.filter[!duplicated(tT.filter$Gene.ID),]
tT.deGenes <- tT.filter[tT.filter$adj.P.Val < 0.001, ]
tT.deGenes <- tT.deGenes[abs(tT.deGenes$logFC) >1,]
tT.deGenes



tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)
deKID    <- translateGeneID2KEGGID(tT.de.names)
allKID   <- translateGeneID2KEGGID(tT.all.names)
allKID

tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 5000,0.7)
tT.pathways.clean<- tT.pathways[tT.pathways$`disturbance index` !=0,]
tT.pathways.clean$CDIST  <- p.adjust(as.numeric(as.character(
                                    tT.pathways.clean$`causal Disturbance`))
                                    ,method = "fdr")
tT.pathways.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                    (tT.pathways.clean$P_ORA)),method = "fdr")


tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,]
tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]
head(tT.pathways.clean[order(tT.pathways.clean$CDIST),],10)

