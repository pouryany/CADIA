CADIA Accompanying Manual and Walk Through
================

Package Installation and Preparation
------------------------------------

Provide an package description

Provide a session info

Example Data Analysis
---------------------

If you have a prepared list of differential expression data skip to the next part.

Here, we will use the standard NCBI automatically generated code to import and calculate differential expressions. Make sure the libraries are installed and updated.

``` r
# The lines are generated from NCBI GEO portal
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
```

Selecting differentially expressed genes and proper filtering.

``` r
tT.filter    <- tT[!is.na(tT$Gene.ID),]
tT.filter    <- tT.filter[!duplicated(tT.filter$Gene.ID),]
tT.deGenes   <- tT.filter[tT.filter$adj.P.Val < 0.05, ]
tT.deGenes   <- tT.deGenes[abs(tT.deGenes$logFC) >1,]
tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)
```

Pathway Enrichment Analysis with CADIA
--------------------------------------

Just run the following code.

``` r
library(CADIA)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following object is masked from 'package:Biobase':
    ## 
    ##     combine

    ## The following objects are masked from 'package:BiocGenerics':
    ## 
    ##     combine, intersect, setdiff, union

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
set.seed(1)
tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 10000)
head(tT.pathways)
```

    ## # A tibble: 6 x 10
    ##   Name  nodes edges   P_ORA `No. DE`   P_SSC `causal Disturb…   cadia
    ##   <chr> <dbl> <dbl>   <dbl>    <dbl>   <dbl>            <dbl>   <dbl>
    ## 1 Micr…    54    98 3.66e-8       16 2.65e-1      0.000000189 2.70e-5
    ## 2 Oocy…    17    82 3.13e-4       10 3.00e-4      0.00000162  1.09e-4
    ## 3 p53 …    80   120 2.83e-7       10 4.80e-1      0.00000228  1.09e-4
    ## 4 PI3K…    56    55 7.34e-3       17 2.50e-3      0.000219    7.82e-3
    ## 5 Ras …    47    13 3.65e-1        8 2.00e-4      0.000768    2.19e-2
    ## 6 Foca…    43    28 1.01e-2       12 9.80e-3      0.00101     2.41e-2
    ## # ... with 2 more variables: ORAFDR <dbl>, KEGGID <chr>

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
