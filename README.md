CADIA Accompanying Manual and Walk Through
================

File guides
-----------

This page contains the code and documentation of CADIA package.

In the folder KEGGPATHS, you can find the XML files of the original pathways. In the folder TestCodes, you can find the test cases and codes for the accompaniying publications. The rest of this document showcases how to use CADIA.

Package Installation and Preparation
------------------------------------

You can install CADIA using devtool library. To install devtool run the following commnad:

``` r
install.packages("devtool")
```

After having `devtools` installed, you can install `CADIA` from its github repository

``` r
library(devtools)
devtools::install_github("pouryany/CADIA")
```

Example Data Analysis
---------------------

If you have a prepared list of differential expression data skip to the next part.

Here, we will use an example data from an ovarian cancer data set from NCBI GEO (GSE12172). Make sure the libraries are installed and updated. Selecting differentially expressed genes and proper filtering.

``` r
library(CADIA)

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

Just run the following code. The first argument is the ENTREZ ID of the differentially expressed genes. The second argument is the ENTREZ ID of all genes from the experiment.

``` r
## CADIA uses random sampling, make sure to set the random seed for reproducing
set.seed(1)
tT.pathways <- causalDisturbance(tT.de.names,tT.all.names,iter = 10000)
head(tT.pathways)
```

    ## # A tibble: 6 x 10
    ##   Name  nodes edges   P_ORA `No. DE`   P_SSC `causal Disturb…   cadia
    ##   <chr> <fct> <fct>   <dbl>    <dbl>   <dbl>            <dbl>   <dbl>
    ## 1 Micr… 299   518   3.66e-8       16 2.65e-1      0.000000189 2.70e-5
    ## 2 Oocy… 124   424   3.13e-4       10 3.00e-4      0.00000162  1.09e-4
    ## 3 p53 … 69    84    2.83e-7       10 4.80e-1      0.00000228  1.09e-4
    ## 4 PI3K… 341   2677  7.34e-3       17 2.50e-3      0.000219    7.82e-3
    ## 5 Ras … 229   1396  3.65e-1        8 2.00e-4      0.000768    2.19e-2
    ## 6 Foca… 201   1854  1.01e-2       12 9.80e-3      0.00101     2.41e-2
    ## # ... with 2 more variables: ORAFDR <dbl>, KEGGID <chr>

In the above results, P\_ORA is the p-value of hypergeometric test. P\_SSC is the p-value of Source-Sink Centrality, Causal Disturbance is the combined p-value, and CADIA is the FDR corrected causal disturbance.
