# Version info: R 3.2.3, Biobase 2.30.0, GEOquery 2.40.0, limma 3.26.8
# R scripts generated  Mon Mar 19 09:24:25 EDT 2018

################################################################
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE21510", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- paste0("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",
               "XXXX00000000000000000000000001111111111111111111")
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
tT.deGenes <- tT.filter[tT.filter$adj.P.Val < 0.005, ]
tT.deGenes <- tT.deGenes[abs(tT.deGenes$logFC) >1,]
tT.deGenes


# Saving the list of Differentially expressed genes
write.table(unlist(tT.deGenes$Gene.ID), "Colorectad_tTDEG", col.names = F,
            row.names = F)


tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)
deKID    <- translateGeneID2KEGGID(tT.de.names)
allKID   <- translateGeneID2KEGGID(tT.all.names)


set.seed(1)
tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 10000,
                                 alpha = 0.1,statEval = 1)
tT.pathways.clean<- tT.pathways #[tT.pathways$`disturbance index` ==0,]
tT.pathways[is.na(tT.pathways$`disturbance index`),]
tT.pathways.clean <- tT.pathways[!is.na(tT.pathways$`disturbance index`),]

tT.pathways.clean$CDIST  <- p.adjust(as.numeric(as.character(
                                     tT.pathways.clean$`causal Disturbance`))
                                     ,method = "fdr")
tT.pathways.clean$ORAFDR <- p.adjust(as.numeric(as.character
                                                (tT.pathways.clean$P_ORA)),method = "fdr")


tT.pathways.clean[tT.pathways.clean$CDIST < 0.05,]
tT.pathways.clean[tT.pathways.clean$ORAFDR <0.05,]





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



### Processing for Gene Set Enrichment analysis.


library(stringr)
library(gage)
library(CADIA)


# Using Gage library, we define a customized gene set
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
rownames(expdata.clean)
nsample <- grep("G0", fl)
csample <- grep("G1", fl)


gsea.Res   <- gage(expdata.clean, gsets = p.kegg.gsets, ref = nsample,
                   sample = csample, compare = "as.group",same.dir = F,
                   rank.test = F)



gsea.great <- data.frame(gsea.Res$greater)
gsea.great <- tibble::rownames_to_column(gsea.great,"Name")
gsea.great <- select(gsea.great, c("Name","p.val","q.val"))
gsea.great <- dplyr::inner_join(gsea.great,
                                tT.pathways.clean[,c("Name","KEGGID")],
                                by = "Name")

gsea.great %<>% mutate(.,ID = KEGGID) %>%
    select(.,c("Name","ID","p.val","q.val" )) %>%
    filter(.,q.val < 0.1)



gsea.great.rep          <- gsea.great
gsea.great.rep[,c(3,4)] <- mapply(formatC,gsea.great.rep[,c(3,4)],
                                  MoreArgs = list(format = "e", digits = 2))


print(xtable(gsea.great.rep), include.rownames = FALSE)

#library(devtools)
#install_github("ctlab/fgsea")



require(fgsea)
library(ggplot2)

deg.values <- tT.filter$logFC
names(deg.values) <-  tT.filter$Gene.ID

fgseaRes <- fgsea(pathways = p.kegg.gsets,
                  stats = deg.values,
                  minSize=15,
                  maxSize=500,
                  nperm=10000)

fgseaRes[padj < 0.1,]

?fgsea

fgseaRes   <- as_data_frame(fgseaRes)
fgseaRes   <- select(fgseaRes, c("pathway","pval","padj"))
fgseaRes   <- dplyr::inner_join(fgseaRes,
                                tT.pathways.clean[,c("Name","KEGGID")],
                                by = c("pathway" = "Name"))

fgseaRes   %<>% mutate(.,ID = KEGGID) %>%
    select(.,c("pathway","ID","pval","padj" )) %>%
    filter(.,padj < 0.05)



fgseaRes.rep          <- fgseaRes
fgseaRes.rep[,c(3,4)] <- mapply(formatC,fgseaRes[,c(3,4)],
                                MoreArgs = list(format = "e", digits = 2))

nrow(fgseaRes.rep)
print(xtable(fgseaRes.rep), include.rownames = FALSE)















