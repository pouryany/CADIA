library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE4122", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL201", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "2222221202X21120210202101112010101201201011222X022022222120X2212122"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

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
cont.matrix <- makeContrasts(MvsN = "G2-G0", G1-G0,MvsB = "G2-G1",
                             MvsAll ="G2-(G1+G0)/2", levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, coef="MvsAll", adjust="BH", sort.by="B", number=Inf)

vennDiagram(decideTests(fit2))


tT.filter  <- tT[!is.na(tT$Gene.ID),]
tT.filter  <- tT.filter[!duplicated(tT.filter$Gene.ID),]
tT.deGenes <- tT.filter[tT.filter$adj.P.Val < 0.05, ]
tT.deGenes <- tT.deGenes[abs(tT.deGenes$logFC) >1,]
tT.deGenes


tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)
deKID    <- translateGeneID2KEGGID(tT.de.names)
allKID   <- translateGeneID2KEGGID(tT.all.names)

tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 20000, 0.4)
tT.pathways.clean<- tT.pathways #[tT.pathways$`disturbance index` !=0,]
tT.pathways.clean$CDIST  <- p.adjust(as.numeric(as.character(
    tT.pathways.clean$`causal Disturbance`))
    ,method = "fdr")
tT.pathways.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                                (tT.pathways.clean$P_ORA)),method = "fdr")



tT.pathways.clean[tT.pathways.clean$CDIST <0.05,]
tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]




#Exporting results

library(stringr)

tT.pathways.clean$KEGGID <- str_sub(rownames(tT.pathways.clean), end = -5)

rownames(tT.pathways.clean) <- NULL

microOC.cdist  <- tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,]
microOC.ora    <- tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]
sapply(microOC.cdist, mode)


microOC.cdist <- microOC.cdist[order(microOC.cdist$CDIST),c(1,10,4,6,8,9)]
microOC.ora <- microOC.ora[order(microOC.ora$ORAFDR),c(1,10,8,9)]
microOC.cdist[,c(3,4,5,6)] <- mapply(as.character,microOC.cdist[,c(3,4,5,6)])
microOC.cdist[,c(3,4,5,6)] <- mapply(as.numeric,microOC.cdist[,c(3,4,5,6)])

microOC.cdist[,c(3,4,5,6)] <- mapply(formatC,microOC.cdist[,c(3,4,5,6)],
                                      MoreArgs = list(format = "e", digits = 2))
microOC.ora[,c(3,4)]   <- mapply(formatC,microOC.ora[,c(3,4)],
                                  MoreArgs = list(format = "e", digits = 2))

microOC.ora
microOC.cdist


library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

print(xtable(microOC.cdist), include.rownames = FALSE)
print(xtable(microOC.ora), include.rownames = FALSE)




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



