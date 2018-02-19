library(annotate)
library(GEOquery)
library(limma)
require(org.Hs.eg.db)
library(affy)
library(KEGGgraph)
biocLite("hgu133b.db")
library(hgu133b.db)






GSE8397         <- getGEO('GSE8397', destdir = "../../GEOData/", GSEMatrix = T)
park.data       <- GSE8397$`GSE8397-GPL97_series_matrix.txt.gz`
park.pheno      <- phenoData(park.data)
park.ph.Tab     <- pData(park.pheno)
park.sample.id  <- park.ph.Tab$title
park.sample.id  <- strsplit2(as.character(park.sample.id) , split = " ")
park.sample.id  <- park.sample.id[,4]
park.type       <- factor(park.sample.id)
park.exp        <- exprs(park.data)


levels(park.type) <- c("control", "dControl","treatment")


park.design <- model.matrix(~0+park.type)
colnames(park.design) <- c ("control","dControl", "treatment")
park.contrast <- makeContrasts( treatment - control , levels = park.design)
park.fit <- lmFit(park.exp, park.design)
park.TvsC <- contrasts.fit(park.fit, park.contrast)
park.TvsC <- eBayes(park.TvsC)


park.results <- topTable(park.TvsC, sort.by = "none", number = Inf)
park.results$ID <- rownames(park.results)
hgu133.ent <- hgu133bENTREZID
park.results$ENTREZ <- unlist(as.list(hgu133.ent[park.results$ID]))

park.filter  <- park.results[!is.na(park.results$ENTREZ),]
park.filter  <- park.filter[!duplicated(park.filter$ENTREZ),]
park.deGenes <- park.filter[park.filter$adj.P.Val < 0.05, ]
park.deGenes <- park.deGenes[abs(park.deGenes$logFC) >1,]


park.all.names <- as.vector(park.filter$ENTREZ)
park.de.names  <- as.vector(park.deGenes$ENTREZ)
deKID    <- translateGeneID2KEGGID(park.de.names)
allKID   <- translateGeneID2KEGGID(park.all.names)

park.pathways <- causalDisturbance(park.de.names,park.all.names,iter = 5000)
park.pathways.clean<- park.pathways[park.pathways$`disturbance index` !=0,]
park.pathways.clean


park.pathways.clean$CDIST <- p.adjust(as.numeric(
    as.character(park.pathways.clean$`causal Disturbance`)), method = "BH")
as.numeric(as.character(park.pathways.clean$P_ORA))

park.pathways.clean[as.numeric(as.character(park.pathways.clean$CDIST)) < 0.5,]




