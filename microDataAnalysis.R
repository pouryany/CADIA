library(annotate)
library(GEOquery)
library(limma)
biocLite("affy")
biocLite("hgu95av2.db")

require(org.Hs.eg.db)
library(affy)
library(KEGGgraph)


GSE13601 <- getGEO('GSE13601', destdir="../../GEOData/", GSEMatrix = T)


# These commands apply when GSEMatrix is F
Meta(GSE13601)
names(GSMList(GSE13601))
##

GSE.data <- GSE13601$GSE13601_series_matrix.txt.gz
GSE.pheno <- phenoData(GSE.data)

GSE.pheno.table <- pData(GSE.pheno)
GSE.pheno.table$source_name_ch1
GSE.pheno.table$source_name_ch1
GSE.exp<-exprs(GSE13601$GSE13601_series_matrix.txt.gz)
GSE.log <- log2(GSE.exp)



GSE.sample.id <- GSE.pheno.table$title
GSE.sample.id <- strsplit2(as.character(GSE.sample.id),split = " ")
GSE.sample.id <- GSE.sample.id[,2]
GSE.patient   <- factor(GSE.sample.id)
GSE.tissue    <- factor(GSE.pheno.table$source_name_ch1)









##Paired analysis for the samples

        GSE.sample.id <- GSE.pheno.table$title
        GSE.sample.id <- strsplit2(as.character(GSE.sample.id),split = " ")
        GSE.sample.id <- GSE.sample.id[,2]
        GSE.patient   <- factor(GSE.sample.id)
        GSE.tissue    <- factor(GSE.pheno.table$source_name_ch1)


        GSE.design2      <- model.matrix(~ GSE.tissue+ GSE.patient)
        contrast.paired  <- makeContrasts( GSE.tissueTumor - GSE.tissueNormal
                                           , levels = GSE.design2)
        colnames(GSE.design2)[1] <- "Intercept"
        GSE.fit2 <- lmFit(GSE.log,GSE.design2)
        GSE.paired <- contrasts.fit(GSE.fit2, contrast.paired)
        GSE.paired <- eBayes(GSE.paired)
        GSE.paired.results <- topTable(GSE.paired, sort.by = "none", number = Inf)
        GSE.paired.results$ID <- row.names(GSE.paired.results)
        hgu.ent <-hgu95av2ENTREZID
        list.hgu.ent <- as.list(hgu.ent)
        GSE.paired.results$ENTREZ   <- unlist(list.hgu.ent[GSE.paired.results$ID])
        GSE.paired.filter           <- GSE.paired.results[!is.na(GSE.paired.results$ENTREZ),]
        GSE.paired.filter  <- GSE.paired.filter[!duplicated(GSE.paired.filter$ENTREZ),]
        GSE.paired.deGenes <- GSE.paired.filter[GSE.paired.filter$adj.P.Val < 0.05, ]
        GSE.paired.deGenes <- GSE.paired.filter[abs(GSE.paired.filter$logFC) > 2.0,]
        GSE.allGene.names  <- as.vector(GSE.paired.filter$ENTREZ)
        GSE.deGenes.names  <- as.vector(GSE.paired.deGenes$ENTREZ)
        deKID              <- translateGeneID2KEGGID(GSE.deGenes.names)
        allKID             <- translateGeneID2KEGGID(GSE.allGene.names)

        GSE.paired.Pathways <- causalDisturbance(GSE.deGenes.names,GSE.allGene.names,iter = 500)
        GSE.paired.Pathways.Clean <- GSE.paired.Pathways[GSE.paired.Pathways$`No. DE` !=0,]
        GSE.paired.Pathways.Clean$CDIST <- p.adjust(as.numeric(
            as.character(GSE.paired.Pathways.Clean$`causal Disturbance`)), method = "BH")
        as.numeric(as.character(GSE.paired.Pathways.Clean$P_ORA))

        GSE.paired.Pathways.Clean[as.numeric(as.character(GSE.paired.Pathways.Clean$CDIST)) < 0.5,]






 # just comparing normal and tumors.


        GSE.design <- model.matrix(~0+GSE.pheno.table$source_name_ch1)
        head(GSE.design)
        colnames(GSE.design) <- c("LN", "Normal", "Tumor")
        contrast.TvsN  <- makeContrasts( Tumor - Normal , levels = GSE.design)
        GSE.fit <- lmFit(GSE.log,GSE.design)

        GSE.TvsN <- contrasts.fit(GSE.fit, contrast.TvsN)
        GSE.TvsN <- eBayes(GSE.TvsN)
        hist(GSE.TvsN$p.value)
        GSE.results <- topTable(GSE.TvsN,sort.by = "none",
                                adjust.method = "BH", number = Inf)

        GSE.results$ID <- row.names(GSE.results)

        hgu.ent <-hgu95av2ENTREZID
        list.hgu.ent <- as.list(hgu.ent)
        GSE.results$ENTREZ<-unlist(as.list(list.hgu.ent[GSE.results$ID]))

        GSE.filter  <- GSE.results[!is.na(GSE.results$ENTREZ),]
        GSE.filter  <- GSE.filter[!duplicated(GSE.filter$ENTREZ),]
        GSE.deGenes <- GSE.filter[GSE.filter$adj.P.Val < 0.05, ]
        GSE.deGenes <- GSE.filter[abs(GSE.filter$logFC) > 2.0,]
        GSE.deGenes

        GSE.allGene.names <- as.vector(GSE.filter$ENTREZ)
        GSE.deGenes.names <- as.vector(GSE.deGenes$ENTREZ)
        deKID    <- translateGeneID2KEGGID(GSE.deGenes.names)
        allKID   <- translateGeneID2KEGGID(GSE.allGene.names)




        GSE.Pathways <- causalDisturbance(GSE.deGenes.names,GSE.allGene.names,iter = 1000)


        GSE.Pathways

        GSE.Pathways.Clean <- GSE.Pathways[GSE.Pathways$`No. DE` !=0,]

        GSE.Pathways.Clean$CDIST <- p.adjust(as.numeric(
            as.character(GSE.Pathways.Clean$`causal Disturbance`)), method = "BH")
        as.numeric(as.character(GSE.Pathways.Clean$P_ORA))

        GSE.Pathways.Clean[as.numeric(as.character(GSE.Pathways.Clean$CDIST)) < 0.5,]






# Have to be careful for the high-high expressions. Lets try now

e<-exprs(GSE13601$GSE13601_series_matrix.txt.gz)
boxplot(e)
test.e <- log2(e)

boxplot(test2.e, col=as.numeric(GSE.pheno.table$source_name_ch1)+1)


# MA plot

Index <- as.numeric(GSE.pheno.table$source_name_ch1)
 d <- rowMeans(test.e[,Index==2]) - rowMeans(test.e[, Index==3])
 a <- rowMeans(test.e)
 smoothScatter(a, d, main="MA plot", xlab="A", ylab="M")
 abline(h=c(-1,1), col="red")




 #library(preprocessCore)

test2.e <- normalize.quantiles(test.e)

d <- rowMeans(test2.e[,Index==2]) - rowMeans(test2.e[, Index==3])
a <- rowMeans(test2.e)
smoothScatter(a, d, main="MA plot", xlab="A", ylab="M")
abline(h=c(-1,1), col="red")




f <- rma(e)
hist(e[,1], xlab = "log")
















#### removal of samples is the same

# {
#    GSE.data <- GSE13601$GSE13601_series_matrix.txt.gz
#    GSE.data <- GSE.data[,-c(1,10,17,18,23,28,29,32,33,37,44,45,46,47,48,53,58)]
#    GSE.pheno <- phenoData(GSE.data)
#
#    GSE.pheno.table <- pData(GSE.pheno)
#    GSE.pheno.table$source_name_ch1
#    GSE.pheno.table$source_name_ch1
#    GSE.exp<-exprs(GSE.data)
#    GSE.log <- log2(GSE.exp)
#
#
#
#    GSE.sample.id <- GSE.pheno.table$title
#    GSE.sample.id <- strsplit2(as.character(GSE.sample.id),split = " ")
#    GSE.sample.id <- GSE.sample.id[,2]
#    GSE.patient   <- factor(GSE.sample.id)
#    GSE.tissue    <- factor(GSE.pheno.table$source_name_ch1)
#    GSE.tissue
#    GSE.design2   <- model.matrix(~ GSE.tissue+ GSE.patient)
#
#
#     GSE.fit2 <- lmFit(GSE.log,GSE.design2)
#
#     GSE.fit2 <- eBayes(GSE.fit2)
#    topTable(GSE.fit2,coef= "GSE.tissueTumor", sort.by = "none")
#     # GSE.paired
#    # GSE.paired.results <- topTable(GSE.paired, sort.by = "none", number = Inf)
#    # GSE.paired.results[(GSE.paired.results$adj.P.Val <0.001)
#    #                    & abs(GSE.paired.results$logFC) > 2,]
#
#
#
#
#
#
#
#
#
#
#
#
#    }
