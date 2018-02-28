# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Thu Feb 22 13:16:21 EST 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE6631", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL8300", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "10101010101010101010101010101010101010101010"
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

tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 20000, 0.4)
tT.pathways.clean<- tT.pathways[tT.pathways$`disturbance index` !=0,]
tT.pathways.clean$CDIST  <- p.adjust(as.numeric(as.character(
    tT.pathways.clean$`causal Disturbance`))
    ,method = "fdr")
tT.pathways.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                                (tT.pathways.clean$P_ORA)),method = "fdr")


tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,]
tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]

head(tT.pathways.clean[order(tT.pathways.clean$CDIST),],10)


## Exporting Results

tT.pathways.clean$KEGGID <- str_sub(rownames(tT.pathways.clean), end = -5)

rownames(tT.pathways.clean) <- NULL

neckC.cdist  <- tT.pathways.clean[tT.pathways.clean$CDIST < 0.1,]
neckC.ora    <- tT.pathways.clean[tT.pathways.clean$ORAFDR <0.1,]
sapply(neckC.cdist, mode)


neckC.cdist <- neckC.cdist[order(neckC.cdist$CDIST),c(1,10,4,6,8,9)]
neckC.ora <- neckC.ora[order(neckC.ora$ORAFDR),c(1,10,8,9)]
neckC.cdist[,c(3,4,5,6)] <- mapply(as.character,neckC.cdist[,c(3,4,5,6)])
neckC.cdist[,c(3,4,5,6)] <- mapply(as.numeric,neckC.cdist[,c(3,4,5,6)])

neckC.cdist[,c(3,4,5,6)] <- mapply(formatC,neckC.cdist[,c(3,4,5,6)],
                                      MoreArgs = list(format = "e", digits = 2))
neckC.ora[,c(3,4)]   <- mapply(formatC,neckC.ora[,c(3,4)],
                                  MoreArgs = list(format = "e", digits = 2))

neckC.ora
neckC.cdist


library(xtable)
options(xtable.floating = FALSE)
options(xtable.timestamp = "")

print(xtable(neckC.cdist), include.rownames = FALSE)
print(xtable(neckC.ora), include.rownames = FALSE)









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

nsample <- grep("G0", fl)
csample <- grep("G1", fl)


a <- gage(expdata.clean, gsets = p.kegg.gsets, ref = nsample, sample = csample,
          compare = "unpaired",same.dir = F, saaTest = gs.KSTest)

head(a$greater[,1:5],10)


