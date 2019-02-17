CADIA Accompanying Manual and Walk Through
================
Pourya Naderi Yeganeh
2019-02-16

Overview
========

This document contains the walk-through/tutorial for the Causal Disturbance Analysis (CADIA) as described in the corresponding publication by Naderi and Mostafavi (References to be provided). CADIA is an enrichment analysis tool for interpreting gene perturbations by contrasting them with underlying graphs of annotated pathways. This program takes an input list of differentially expressed gene IDs and a gene universe and produces p-values that indicate pathway enrichments.

Dependencies and Installation guide
===================================

CADIA depends on the following packages: `KEGGgraph`, `RBGL`, `graph`, `stringr`, `dplyr`, `magrittr`. Make sure the packages are installed before using CADIA. To install the packages you can use the following code chunk.

``` r
dependencies  <- c("KEGGgraph", "RBGL", "graph", "stringr", "dplyr", "magrittr")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Please note the parameters of BiocManager::install() and modify accordingly

for (i in dependencies) {
    if(!require(i))
        BiocManager::install(i,update = F)
}
```

Package Installation and Preparation
------------------------------------

You can install CADIA using devtools library. To install devtool run the following commnad:

``` r
install.packages("devtools")
```

After having `devtools` installed, you can install `CADIA` from its github repository

``` r
library(devtools)
devtools::install_github("pouryany/CADIA")
```

File guides
-----------

The github page of CADIA packages contains associated the code and documentation. The folder data-raw and its subfolders contains relevant data and code accompanying the package and the original manuscripts. In particular, KEGGPATHS contains the raw unprocessed XML files of the original pathways. The TestCodes folder contains the test cases and codes for the accompaniying publications. Additional instructions and guides regarding the specific test cases and codes are provided in the GuideMe.R document inthe data-raw folder.

Enrichment Analysis with CADIA
==============================

The current version of CADIA works by taking two inputs of differentially expressed genes (DEG) and the gene universe. The output of CADIA is a list of KEGG pathways along with a group of p-values that describe the association of the pathway with a the DEG.

This document contains a test case of CADIA using an ovarian cancer dataset by Bowtell and colleagues, with 60 High-grade serous ovarian cancer and 30 Low malignant potential tumors. This data is available from NCBI-GEO portal through accession code GSE12172 (Anglesio et al. 2008). The datasets platform us affymetrix HG-U133b.

Finding differentially expressed genes
--------------------------------------

NCBI-GEO portal provides and R-code generation tool for differential expression analysis of its datasets. The below code shows the automatically generated script for differential expression analysis of the GSE12172 dataset. In general, it uses the `limma` package to contrast the normalized expressions of two predefined subsets of the samples. For running the next code chunk you would need to have these packages installed: `Biobase`, `GEOquery`, and `limma`. If you already have a `limma` topTable object you may skip to the next section where we use a predefined table of this analysis.

``` r
# The following lines are based on the code generated from NCBI GEO portal
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

# group names for all samples, 1 For cancer and 0 for non-cancer
gsms <- paste0("11111010110111001100111111110101001010101111101111",
               "1101111010011110011110000011111010010111")
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# log2 transform
ex   <- exprs(gset)
qx   <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
    (qx[6]-qx[1] > 50 && qx[2] > 0) ||
    (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml              <- paste("G", sml, sep="")    # set group names
fl               <- as.factor(sml)
gset$description <- fl
design           <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit              <- lmFit(gset, design)
cont.matrix      <- makeContrasts(G1-G0, levels=design)
fit2             <- contrasts.fit(fit, cont.matrix)
fit2             <- eBayes(fit2, 0.01)
tT               <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
```

Preparing differentially expressed genes
----------------------------------------

Here, we will an topTable object that was produced from an ovarian cancer data set from NCBI GEO (GSE12172) (Anglesio et al. 2008). The instructions for producing the topTable object are provided in the previous sections. The next step for enrichment analysis is to identify a subset of differentially expressed genes. The specific platform of this data has multiple probes for each annotated genes. Also, some of the probes do not correspond to any annotated gene. After appropriate filtering for such instances, we use an adjusted p-value and fold-change criteria to identify the subset of DEG. Note that the end results of this analysis step are two lists `tT.all.names` and `tT.de.names`.

``` r
library(CADIA)
```

    ## 
    ## Attaching package: 'CADIA'

    ## The following object is masked _by_ '.GlobalEnv':
    ## 
    ##     tT

``` r
data(tT)
# Filter out unannotated and duplicated genes
tT.filter    <- tT[!is.na(tT$Gene.ID),]
tT.filter    <- tT.filter[!duplicated(tT.filter$Gene.ID),]

# Select a p-value and fold-change cut-off for differential expressions
tT.deGenes   <- tT.filter[tT.filter$adj.P.Val < 0.05, ]
tT.deGenes   <- tT.deGenes[abs(tT.deGenes$logFC) >1,]

# Parameters needed for CADIA. The name of all genes and DE genes.
tT.all.names <- as.vector(tT.filter$Gene.ID)
tT.de.names  <- as.vector(tT.deGenes$Gene.ID)
```

Pathway Enrichment Analysis with CADIA
--------------------------------------

The current version of CADIA only works with ENTREZ IDs. To run the enrichment analysis, simply provide the function `causalDisturbance()` with minimum of two inputs The first argument is the ENTREZ ID of the DEG. The second argument is the ENTREZ ID of all genes from the experiment. If your data has an output other than ENTREZ IDs refer to the subsection at the end of this sections.

Note that, as described in the manuscript, CADIA uses random sampling for calculating the outputs. To ensure the consistency and reproducibility of your results, always the random seed before running the method.

``` r
set.seed(1)
cadia.res <- causalDisturbance(tT.de.names,tT.all.names,iter = 10000)
```

The following are the additional arguments for `causalDisturbance()`. `iter` denotes the number of iterations for random bootstrap sampling (preferrably larger than 2000). `alpha` is the dampening factor of Source-Sink centrality (described in the manuscript, defaul = 0.1), `beta` is the relative importance of sink centrality compared to the source centrality (default = 1). `statEval` denotes whether to use a product or summation for calculating the aggeregate centrality of the perturbations (default = 1, product-based aggregation). `fdrMethod` is the choice multiple hypothesis correction method (default = "B", see `p.adjust` function documentaion from stats package for options). `verbose` denotes whether to generate progress report during calculations.

``` r
head(cadia.res)
```

    ## # A tibble: 6 x 10
    ##   Name  nodes edges   P_ORA `No. DE`   P_SSC `causal Disturb…   cadia
    ##   <chr> <fct> <fct>   <dbl> <fct>      <dbl>            <dbl>   <dbl>
    ## 1 Micr… 299   518   3.66e-8 30       2.65e-1      0.000000189 2.70e-5
    ## 2 Oocy… 124   424   3.13e-4 18       3.00e-4      0.00000162  1.09e-4
    ## 3 p53 … 69    84    2.83e-7 18       4.80e-1      0.00000228  1.09e-4
    ## 4 PI3K… 341   2677  7.34e-3 34       2.50e-3      0.000219    7.82e-3
    ## 5 Ras … 229   1396  3.65e-1 16       2.00e-4      0.000768    2.19e-2
    ## 6 Foca… 201   1854  1.01e-2 22       9.80e-3      0.00101     2.41e-2
    ## # … with 2 more variables: ORAFDR <dbl>, KEGGID <chr>

`cadia.res` is output table. The column `cadia` is the FDR corrected causal disturbance, which indicates the statistical significance of a pathway enrichment. The output table is sorted based on this value and can be taken as the enrichment resutls. Additional outputs of the method denote different aspects of the analysis. `P_ORA` is the p-value of over-representation analysis through hypergeometric test. `P_SSC` is the topological evidence, p-value of Source-Sink Centrality as described in the original manuscript. The column `Causal Disturbance` is the combined `P_SSC` and `P_ORA` using Fisher's method. The column `ORA` is the adjusted `P_ORA`. `KEGGID` is the associated ID of a pathway in the KEGG Database.

`cadia.res` is a tibble and, thus, can be manipulated using tibble related methods. For example, to create a table with only pathways with `cadia < 0.05`.

``` r
library(dplyr)
top.res <- cadia.res %>% dplyr::filter(.,cadia < 0.05) %>% select(., Name, KEGGID)
head(top.res)
```

    ## # A tibble: 6 x 2
    ##   Name                       KEGGID
    ##   <chr>                      <chr> 
    ## 1 MicroRNAs in cancer        05206 
    ## 2 Oocyte meiosis             04114 
    ## 3 p53 signaling pathway      04115 
    ## 4 PI3K-Akt signaling pathway 04151 
    ## 5 Ras signaling pathway      04014 
    ## 6 Focal adhesion             04510

Translating and preparing inputs for CADIA
------------------------------------------

As mentioned, the current version of CADIA only works with ENTREZ IDs. If the results of your analysis is in some other format, e.g. symbol and ENSEMBEL, we suggest using the functions of `clusterProfiler` package for appropriate translations (Yu et al. 2012). To do this, you would need appropriate library installations. The code below depicts a procedure for translating between different formats in the hypothetical case where your experimental results are in gene symbol notation.

``` r
library(clusterProfiler)
library(org.Hs.eg.db)
library(CADIA)
# the object geneList contains a list of all genes in the universe

# the object deGenes contains a list of differentially expressed genes.


gene.df    <- bitr(geneList, fromType = "SYMBOL",
                   toType = c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)
deGenes.df <- bitr(deGenes, fromType = "SYMBOL",
                   toType = c("ENTREZID","ENSEMBL"),
                   OrgDb = org.Hs.eg.db)

set.seed(1)
cadia.res  <- CADIA::causalDisturbance(deGenes.df$ENTREZID,gene.df$ENTREZID,
                                       iter = 5000)
```

Additional functions in CADIA package
=====================================

The package CADIA provides additional functionality for pathway enrichment analysis and other applications. This section provides a brief overview of the functions.

Internal data
-------------

The graph objects of KEGG pathways that are used in CADIA can be accessed using the internal data of package. The data `pathways.collection` contains the processed pathway graphs in the graphNEL format.

``` r
library(CADIA)
data("pathways.collection")
pathways.collection[[1]]
```

    ## A graphNEL graph with directed edges
    ## Number of Nodes = 229 
    ## Number of Edges = 1396

If you are interested in getting a list of pathways that are analyzed by CADIA you can use the built-in function `cadia.paths()`.

``` r
head(CADIA::cadia.paths())
```

    ##           KEGG_IDs pathways.collection.names
    ## 04014.xml    04014     Ras signaling pathway
    ## 04015.xml    04015    Rap1 signaling pathway
    ## 04010.xml    04010    MAPK signaling pathway
    ## 04012.xml    04012    ErbB signaling pathway
    ## 04310.xml    04310     Wnt signaling pathway
    ## 04330.xml    04330   Notch signaling pathway

One might be intresed to retrieve the differentially expressed genes associated with the significantly enriched pathways, or a subset of them. The function `geneReport()` faciliates this operation by returning a list containing the pathways and the list of their DEG (in ENTREZ format) concatenated in the rows. See the following.

``` r
reports <- geneReport(tT.de.names,top.res$KEGGID)
rownames(reports) <- NULL
head(reports)
```

    ##                   path.names
    ## 1        MicroRNAs in cancer
    ## 2             Oocyte meiosis
    ## 3      p53 signaling pathway
    ## 4 PI3K-Akt signaling pathway
    ## 5      Ras signaling pathway
    ## 6             Focal adhesion
    ##                                                                                                                                                                      genes
    ## 1                 993/995/7422/4170/1545/3925/407040/4233/4318/595/898/9134/27086/7430/7078/10253/4893/3667/8660/2146/113130/9493/1026/23414/7329/5743/2261/1029/4193/5156
    ## 2                                                                                    983/4085/9133/29945/991/9232/891/898/9134/995/5347/6790/91860/9700/5241/3480/3479/699
    ## 3                                                                                1111/898/9134/6241/1643/3479/1029/1026/64393/51512/8795/983/4193/595/64065/27244/891/9133
    ## 4 3667/1969/2261/3480/3643/4233/5156/10161/4193/595/5105/2247/2259/3479/7422/4893/5521/5563/1026/3673/3675/3691/3693/1298/1311/3908/3914/3918/5747/3574/4170/4602/898/9134
    ## 5                                                                                        4893/2247/2259/3479/7422/1969/2261/3480/3643/4233/5156/51196/5602/5924/91860/5063
    ## 6                                                             3691/3693/3673/3675/55742/5747/824/1298/1311/3908/3914/3918/858/595/3480/4233/5156/3479/7422/5602/5063/10451

The graph analysis functionality of CADIA is implemented separately for those who wish to utilize them in different lines of research. The function `source.sink.centrality()` implements the concept of Source/Sink Centrality as describe in the original manuscript. This function returns a matrix whose individual elements corresponding to row i and column j can be interpreted as the influence of node i on node j. One can alternatively calculate the centrality of individual nodes by using the function `rowSums()`

``` r
library(CADIA)
library(graph)
data("pathways.collection")
test.graph    <- pathways.collection[[1]]
test.matrix   <- as(test.graph,"matrix")
ssc.influence <- CADIA::source.sink.centrality(test.matrix,alpha = 0.1, beta =1)
head(rowSums(ssc.influence))
```

    ##  hsa:5594  hsa:5595  hsa:5604  hsa:5605  hsa:5894 hsa:22800 
    ##  4.664992  4.664992  3.004958  3.004958  4.585584  7.524828

The notion of the aggregate score in the original paper can be used in applications where one is interested in computing a centrality score for a subset of nodes. The following provides a showcase of how to calculate this notion of subgraph centrality using the built-in function `pathSampler()`.

``` r
test.nodes <- graph::nodes(test.graph)
set.seed(1)
subs.nodes <- sample(test.nodes,10, replace = F)

set.seed(1)
subs.prob  <- pathSampler(inputGraph = test.graph ,deKID = subs.nodes, 
                          iterationNo = 10000, alpha = 0.1, beta = 1,
                          statEval = 1  )
subs.prob
```

    ## [1] 0.8295

``` r
# sub.prob is a probability, one can turn it into a score as following

-log(subs.prob,base = 10)
```

    ## [1] 0.08118361

References
==========

Anglesio, Michael S, Jeremy M Arnold, Joshy George, Anna V Tinker, Richard Tothill, Nic Waddell, Lisa Simms, et al. 2008. “Mutation of Erbb2 Provides a Novel Alternative Mechanism for the Ubiquitous Activation of Ras-Mapk in Ovarian Serous Low Malignant Potential Tumors.” *Molecular Cancer Research* 6 (11). AACR: 1678–90.

Yu, Guangchuang, Li-Gen Wang, Yanyan Han, and Qing-Yu He. 2012. “ClusterProfiler: An R Package for Comparing Biological Themes Among Gene Clusters.” *Omics: A Journal of Integrative Biology* 16 (5). Mary Ann Liebert, Inc. 140 Huguenot Street, 3rd Floor New Rochelle, NY 10801 USA: 284–87.
