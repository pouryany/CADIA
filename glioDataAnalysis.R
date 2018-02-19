# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Sat Feb 17 14:41:42 EST 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE4290", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("121222X02100020201002210022023XX00002X3200X01X20XX",
               "12302000220X02X02X2X2X3X3X222XXXXX2X0XX23X22222122",
               "22223213XX232232X1X2X221222112X2321XX22XXX2212X212",
               "1X2232X2X2222XX1212X2X22222X2X")
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




head(gset)


# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts("OGvsN" = "G3-G0","GliovsN"= "G2-G0", "AstrovsN" = "G1-G0",
                             "AllGliovsN" = (G3 +G2 + G1)/3 -G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2,coef ="AllGliovsN" , adjust="BH", sort.by="B", number=Inf)


tT.filter  <- tT[!is.na(tT$Gene.ID),]
tT.filter  <- tT.filter[!duplicated(tT.filter$Gene.ID),]
tT.deGenes <- tT.filter[tT.filter$adj.P.Val < 0.0001, ]
tT.deGenes <- tT.deGenes[abs(tT.deGenes$logFC) >2,]
tT.deGenes

tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)
deKID    <- translateGeneID2KEGGID(tT.de.names)
allKID   <- translateGeneID2KEGGID(tT.all.names)

tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 10000)
tT.pathways.clean<- tT.pathways[tT.pathways$`disturbance index` !=0,]
tT.pathways.clean$CDIST  <- p.adjust(as.numeric(as.character(
                                    tT.pathways.clean$`causal Disturbance`))
                                    ,method = "fdr")
tT.pathways.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                    (tT.pathways.clean$P_ORA)),method = "fdr")


tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,]
tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]

head(tT.pathways.clean[order(tT.pathways.clean$`causal Disturbance`),],20)
