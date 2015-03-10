# Response to Buettner, Natarajan "Computational analysis of cell-to-cell heterogeneity in single-cell RNA-sequencing data reveals hidden subpopulations of cells"
### Andrew McDavid, Raphael Gottardo, Greg Finak
#### Fred Hutchinson Cancer Research Center
<!-- using knitrBootstrap::knit_bootstrap -->


## Loading mESC data
We first load the expression data that the authors provide in the paper, although note that the final gene in the mESC data set contained corrupted data in the file included in the Nature supplement.   John Marioni kindly provided us with a version that fixed that corruption, which is what is being used here.

We also load the list of cell cycle genes, and an estimate of the Panama (scLVM) latent variable.  In the remainder of this document we refer to this variable as just the panama latent variable, since that is what BNC are doing in the python code provided.

```r
## Supplementary table 1, sheet 1
ccgenes <- read.xls('data/cell_cycle_genes.xlsx', as.is=TRUE, header=F)$V1
## converted from excel .xlsx, Supplementary Data 2, fixing a corrupted character in final column.
##These are Log10(ercc_normalized_counts+1).
uncorr <- fread('data/mesc_cell_cycle_uncorrected.csv')
## Generated from above using authors' code from https://github.com/PMBio/scLVM
Xpanama <- scan('data/ES_X_ARD_1.txt')
mesc <- as.matrix(uncorr[,-1,with=FALSE])
## using totals of each phase listed in paper
mescCdat <- data.frame(cycle_class=factor(c(rep("G1",59),  rep("S", 58),  rep("G2M",65)), levels=c('G1', 'S', 'G2M')), idx=seq_len(nrow(mesc)), panama=Xpanama, geomsizeUncorr=rowSums(mesc))
```

We'll also need a crosswalk between the gene symbols and ensembl ids provided in `ccgenes`:

```r
eg <- as.data.table(as.data.frame(org.Mm.egSYMBOL2EG))
setkey(eg, gene_id)
ensembleg <- as.data.table(as.data.frame(org.Mm.egENSEMBL))
setkey(ensembleg, gene_id)
fdat <- ensembleg[eg,,nomatch=0]
setkey(fdat, symbol)
```

Now join with symbols provided.  Some symbols map to multiple ensembl ids--we'll just take the first one listed.

```r
mescFdat <- fdat[colnames(mesc),,mult='first']
mescFdat <- mescFdat[,ccRanked:=factor((ensembl_id %in% ccgenes), labels=c('FALSE'='Unranked', 'TRUE'='Ranked'))]
mescFdatDup <- fdat[colnames(mesc)]
```
This affects 77 records.

## Loading mouse T-cell data
We repeat the process as above.

```r
## Supplementary Data 1, sheet 1
T_cell_raw <- fread("data/T_cell_uncorrected.csv")
## sheet 2
T_cell_corrected <- fread("data/T_cell_corrected.csv") 
setnames(T_cell_raw, "V1", "cell_id")
setnames(T_cell_corrected, "V1", "cell_id")
## Supplementary Data 1, sheet 3
cluster <- fread("data/T_cell_cluster.csv")
setnames(cluster, 'Gata3HighCLuster', 'clusterid')
cluster <- cluster[,clusterid:=factor(clusterid)]
## Again, we generated using authors' code from https://github.com/PMBio/scLVM
lv_factor <- scan("data/X_ARD_1.txt")
```

Two gene symbols are duplicated, we'll just take the first one listed.

```r
T_cell_matrix <- as.matrix(T_cell_raw[,-1,with=FALSE])
T_cell_matrix_corrected <- as.matrix(T_cell_corrected[,-1,with=FALSE])
# Remove duplicated gene names (only 2)
T_cell_matrix <- T_cell_matrix[,unique(colnames(T_cell_matrix))]
T_cell_matrix_corrected <- T_cell_matrix_corrected[,unique(colnames(T_cell_matrix_corrected))]
stopifnot(all(colnames(T_cell_matrix_corrected)==colnames(T_cell_matrix)))

TcellCdat <- data.frame(cluster, panama=lv_factor, geomsizeUncorr=rowSums(T_cell_matrix))
```
Again, we'll need to map symbols to ensembl.

```r
TcellFdat <- fdat[colnames(T_cell_matrix),,mult='first']
TcellFdat <- TcellFdat[,ccRanked:=factor((ensembl_id %in% ccgenes), labels=c('FALSE'='Unranked', 'TRUE'='Ranked'))]
TcellFdatDup <- fdat[colnames(T_cell_matrix)]
```

There are 182 cells and 9571 genes in the mESC data and
81 cells and 7071 genes in the T-cell data.


## $R^2$ attributable to Cell Cycle in mESC
First, we explore how much of the variance in the mESC data set can be explained by cell cycle.

Below, we define the $R^2$ as a ratio of linear model deviances (sums of squares) in a full and null model.

```r
lrTestLM <- function(Fit, drop){
    nFit <- update(Fit, as.formula(paste0('. ~ . -', drop)))
    lambda <- (deviance(nFit)-deviance(Fit))
    nDF <- (Fit$rank-nFit$rank)
    data.frame(lambda=lambda, P=1-pchisq(lambda, nDF), pctDev=pmin(lambda/deviance(nFit)*100))
    }
```
This function just compares the residual squared error with and without the factor `drop` (in our case, this will always be cell cycle), and calculates this change (`lambda`) and its ratio over the null deviance as a percentage.

BNC offer an adjustment in their variance decompositions for a factor that reflects technical variability.  We didn't attempt to replicate their calculation of this factor.  Below we compare the unadjusted values to some benchmark values from their paper and find reasonable agreement nonetheless, suggesting that these findings are robust to adjustment for putative sources of technical variability.

We consider regressions of `cycle_class` on the uncorrected log expression values.

```r
afit <- lm(mesc~cycle_class, data=mescCdat)
atest <- lrTestLM(afit, 'cycle_class')
P <- atest$P 
pctDev <- atest$pctDev
lambda <- atest$lambda
devModels <- data.table(model='cycle_class', symbol=colnames(coef(afit)), P, lambda, pctDev)
setkey(devModels, symbol)
setkey(mescFdat, symbol) 
devModels <- devModels[mescFdat,,allow.cartesian=TRUE]
```
We have joined the linear model results to the gene annotation data here.

Next, we report the variance estimates in median/90% gene in ranked/unranked sets, which gives the following result:

```r
devSum <- devModels[,{ 
    Q <- quantile(pctDev, prob=c(.5, .9), na.rm=TRUE)
    list('50%'=Q[1], '90%'=Q[2], Ngenes=.N)
}, keyby=list(ccRanked)]

kable(devSum)
```



|ccRanked |      50%|      90%| Ngenes|
|:--------|--------:|--------:|------:|
|Unranked | 3.479381| 15.15893|   8994|
|Ranked   | 7.968757| 26.93517|    609|

### Comparing to Supplementary Figure 4b
We sought to compare our estimate to a description of these quantities in the paper.  We located an estimate of the effect of cell cycle by reading numbers off of figures S4b/S4d.

From this figure, fewer than 22% (1-7638/9571) of genes have >10% of variance attributable to cell cycle via "gold standard" including both ranked and unranked genes.
We reproduce the rest of the numbers from reading the gene counts annotated above the box plots.
Although this ignores the difference between ranked/unranked genes, we see that  BNC's result lies somewhere between the results for ranked/unranked genes.

```r
## right endpoints of R^2 bins
R2 <- c(0, 10, 20, 30, 40, 50, 90)
## number of genes in bin
counts <- c(0, 7638, 1079, 483, 225, 97, 49)
## number of genes in bin and all bins below
countsBelow <- cumsum(counts)
figureData <- data.frame(R2, countsBelow, pctBelow=countsBelow/max(countsBelow), source='From Figure S4b')
ggplot()+stat_ecdf(data=devModels, mapping=aes(x=pctDev, col=ccRanked)) +geom_line(data=figureData, mapping=aes(x=R2, y=pctBelow, col=source)) + xlab('% R^2 due to cell cycle') + ylab('Cumulative Distribution') + scale_color_discrete("Source") + theme(legend.position=c(.7, .3)) + geom_abline(aes(intercept=c(.5, .9), slope=0), lty=2)
```

![plot of chunk varianceEstimatesFigureS4](figure/varianceEstimatesFigureS4-1.png) 
#### Remark
Since there is no manifest variable that explains the "technical" variance that BNC subtract in their calculations, we think that it makes more sense to consider the undeflated estimates, since these represent the gene expression values that experimenters must use in practice, until a factor can be identified as a manifest variable for the technical variability.

## PCA
### PCA On Uncorrected mESC
PCA on all genes shows that the geometric size `geomsize` is principal component 1, explaining 9%-29% of the variance.

```r
getPCA <- function(exprMat, cData, scale.=FALSE){
    pcOut <- prcomp(exprMat, scale.=scale.)
    geomsize <- rowSums(exprMat)
    pcall <- data.frame(pcOut$x[,1:3], geomsize, cData)
    varExp <- pcOut$sdev^2/sum(pcOut$sdev^2)
    list(pcPairs=pcall, varianceExplained=varExp)
}

mescPca <- getPCA(mesc, mescCdat[,c('cycle_class', 'panama')])
```

We see that PC1 one is highly correlated with geometric size (`geomsize`) and the `panama` factor that we estimated.
We also observe that the cell cycle phase of a cell is correlated with `geomsize`/PC1 and discuss that point at greater length below.

```r
ggpairs(mescPca$pcPairs, color='cycle_class')
```

![plot of chunk plotmescpca](figure/plotmescpca-1.png) 

### PCA on Uncorrected Tcells
In the Tcells, we again see that PC1 one is highly correlated with geometric size (`geomsize`) and the `panama`  factor.

```r
tcpca <- getPCA(T_cell_matrix, TcellCdat[,c('panama', 'clusterid')])
ggpairs(tcpca$pcPairs, columns=c(1:5), color='clusterid')
```

![plot of chunk plotTcpca](figure/plotTcpca-1.png) 


## Percentage variance due to PC1

```r
varExplained <- data.frame('data set'=c('mESC', 'T-cell'),
                           'PC1 Variance Explained'=c(mescPca$varianceExplained[1],
                               tcpca$varianceExplained[1]), check.names=F)
kable(varExplained) 
```



|data set | PC1 Variance Explained|
|:--------|----------------------:|
|mESC     |              0.0892713|
|T-cell   |              0.2895370|

### Remark on T-cell clustering
The T-cell clustering reported in the corrected T-cell data appears to depend substantially on the dimension reduction method used.  We are unable to replicate the clustering in the corrected or uncorrected data using linear PCA, or non-linear t-stochastic neighborhood embedding.  As the exact software version and parameters BNC used for the 'non-linear PCA' were not provided, we didn't attempt to apply that algorithm.

Using ordinary PCA, we found instead that the GATA3 `cluster` variable corresponds roughly to a half-space defined by PC2, but since multiple modes are not evident, we would not detect distinct clusters, but rather a gradient of cells.

```r
tcCorrectedPCA <- getPCA(T_cell_matrix_corrected, TcellCdat[,c('panama', 'clusterid')])
qp <- qplot(data=tcCorrectedPCA$pcPairs, x=PC1, y=PC2, col=clusterid, size=geomsize)+geom_density2d(aes(col=NA))
qp1 <- qp + ggtitle('PCA Corrected')
qp2 <- qp %+% tcpca$pcPairs + ggtitle('PCA Uncorrected')
grid.arrange(qp1, qp2)
```

![plot of chunk PCAforClusters](figure/PCAforClusters-1.png) 

We also try to do the same using `t-SNE` as implemented in the `Rtsne` package, which is a non-linear dimension reduction technique, which has been successfully applied to single-cell data.

The figures above (PCA) and below (t-SNE) do now show any obvious clustering of the cells either before or after correction according to the clusters reported in BNC. An obvious grouping of the cells by geomsize is present.

```r
tsne <- data.frame(tsne=Rtsne(T_cell_matrix, dims = 2, initial_dims = 50, perplexity = 20, theta = 0.5, check_duplicates = TRUE, pca = TRUE)$Y, geomsize=rowSums(T_cell_matrix))
tsneCorr <- data.frame(tsne=Rtsne(T_cell_matrix_corrected, dims = 2, initial_dims = 50, perplexity = 20, theta = 0.5, check_duplicates = TRUE, pca = TRUE)$Y, geomsize=rowSums(T_cell_matrix_corrected)) 
qptsne <- qplot(data=cbind(tsne, TcellCdat), x=tsne.1, y=tsne.2, col=clusterid, size=geomsize) +geom_density2d(aes(col=NA))
qp3 <- qptsne %+% cbind(tsneCorr, TcellCdat) + ggtitle('TSNE Corrected')
qp4 <- qptsne + ggtitle('TSNE Uncorrected')
grid.arrange(qp3, qp4)
```

![plot of chunk TSNEforClusters](figure/TSNEforClusters-1.png) 


## Relationship between `geomsize`, `panama` and cell cycle
From the plots of the PCA in mESC, we observe that average geometric size and average panama factor both depend on the cell cycle.  Here we formally consider how much cell cycle explains in these factors by considering a regression of cell cycle on `geomsizeUncorr` (geometric size in uncorrected data) and cell cycle on `panama`.

```r
models <- c('geomsizeUncorr ~ cycle_class', 'panama ~ cycle_class', 'panama ~ geomsizeUncorr +cycle_class')
R2duetoCC <- sapply(models, function(m){
    l <- lm(as.formula(m), data=mescCdat)
    dr <- anova(l)
    totalss <- sum(dr[,'Sum Sq'])
    ccss <- sum(dr['cycle_class', 'Sum Sq'])
    ccss/totalss
})

kable(data.frame(Model=models, 'R^2 due to cycle class' = R2duetoCC, check.names=FALSE), row.names=FALSE)
```



|Model                                | R^2 due to cycle class|
|:------------------------------------|----------------------:|
|geomsizeUncorr ~ cycle_class         |              0.5313693|
|panama ~ cycle_class                 |              0.6479790|
|panama ~ geomsizeUncorr +cycle_class |              0.0333929|


## Gene Set Enrichment Analysis using Corrected Data
To conclude our analysis, we investigated if the scLVM factor would be successful in eliminating cell cycle signal when testing for differences between the T-cell clusters.

We considered gene set enrichment analysis on the corrected T-cell data.

```r
eSet <- ExpressionSet(t(T_cell_matrix_corrected))
pData(eSet) <- TcellCdat
```

We considered gene set enrichment analysis on the corrected T-cell data.



We used the Broad Institute's "Reactome" module available here: http://www.broadinstitute.org/gsea/msigdb/collections.jsp.

```r
c2_set <- getGmt("data/c2.cp.reactome.v4.0.symbols.gmt")
gene_ids <- geneIds(c2_set)
## # Camera requires gene-indices
design <- model.matrix(~clusterid, eSet)
sets_indices <- ids2indices(gene_ids, toupper(rownames(eSet))) 
res <- camera(eSet, sets_indices, contrast=2, design=design)
```




```r
res$Direction <- ifelse(res$Direction=='Up', 1, -1)
res$set <- row.names(res)
setDT(res)
setkey(res, FDR) 
res <- res[,':='(totalrank=rank(PValue))]
resFDR <- res[FDR<.05,]
nFDR <- nrow(resFDR)
```
In the 34 modules with a FDR-q value 5% or less, a substantial number relate to cell cycle.

Among the top 20 modules, there is little suggestion of enrichment for immue-related modules.  Many modules relate to cell cycling.

```r
kable(resFDR[1:20,list(Direction, FDR, totalrank, set)]) 
```



| Direction|       FDR| totalrank|set                                                                                      |
|---------:|---------:|---------:|:----------------------------------------------------------------------------------------|
|        -1| 0.0000036|         1|REACTOME_MITOTIC_PROMETAPHASE                                                            |
|        -1| 0.0000153|         2|REACTOME_E2F_MEDIATED_REGULATION_OF_DNA_REPLICATION                                      |
|        -1| 0.0000346|         3|REACTOME_DEPOSITION_OF_NEW_CENPA_CONTAINING_NUCLEOSOMES_AT_THE_CENTROMERE                |
|        -1| 0.0000346|         4|REACTOME_G2_M_CHECKPOINTS                                                                |
|        -1| 0.0000951|         5|REACTOME_G1_S_SPECIFIC_TRANSCRIPTION                                                     |
|        -1| 0.0002130|         6|REACTOME_CELL_CYCLE                                                                      |
|        -1| 0.0002834|         7|REACTOME_CHROMOSOME_MAINTENANCE                                                          |
|        -1| 0.0015164|         8|REACTOME_ACTIVATION_OF_ATR_IN_RESPONSE_TO_REPLICATION_STRESS                             |
|        -1| 0.0015164|         9|REACTOME_DNA_REPLICATION                                                                 |
|        -1| 0.0015164|        10|REACTOME_CELL_CYCLE_MITOTIC                                                              |
|        -1| 0.0015164|        11|REACTOME_KINESINS                                                                        |
|        -1| 0.0015164|        12|REACTOME_MITOTIC_M_M_G1_PHASES                                                           |
|        -1| 0.0018450|        13|REACTOME_DOUBLE_STRAND_BREAK_REPAIR                                                      |
|        -1| 0.0020784|        14|REACTOME_CYCLIN_A_B1_ASSOCIATED_EVENTS_DURING_G2_M_TRANSITION                            |
|        -1| 0.0032273|        15|REACTOME_G2_M_DNA_DAMAGE_CHECKPOINT                                                      |
|        -1| 0.0037052|        16|REACTOME_FANCONI_ANEMIA_PATHWAY                                                          |
|        -1| 0.0037108|        17|REACTOME_ACTIVATION_OF_THE_PRE_REPLICATIVE_COMPLEX                                       |
|        -1| 0.0048701|        18|REACTOME_HOMOLOGOUS_RECOMBINATION_REPAIR_OF_REPLICATION_INDEPENDENT_DOUBLE_STRAND_BREAKS |
|        -1| 0.0061528|        19|REACTOME_APC_CDC20_MEDIATED_DEGRADATION_OF_NEK2A                                         |
|        -1| 0.0089812|        20|REACTOME_E2F_ENABLED_INHIBITION_OF_PRE_REPLICATION_COMPLEX_FORMATION                     |


## Gene Set Enrichment Analysis using Counts Per Million and GATA3 level
We explored GSEA on the uncorrected data, as well.
To relax the assumption of faithful normalization via the ERCC spike-ins we calculated the counts per million for each gene.


```r
counts <- 10^(T_cell_matrix)-1
scounts <- rowSums(counts)
logcpm <- log10(counts/scounts*1e6+1)
geomsizecpm <- rowSums(logcpm)
```
Based on the observation from PCA that GATA3 appears to exist more a gradient than in discrete clusters, we considered an alternate GSEA analysis. We used GATA3 expression levels in hopes they might partially define the T cell differentiation state, then removed gata3 from the expression matrix.


```r
eSet_adj <- ExpressionSet(t(logcpm[, setdiff(colnames(logcpm), 'Gata3')]))
gata3 <- logcpm[,'Gata3']
pData(eSet_adj) <- cbind(TcellCdat, geomsizecpm=geomsize, gata3=gata3)
```

Based on observation that `geomsizecpm` still is the leading principal component,  we adjust for it as a nuisance source of variation, and run GSEA, testing gata3:

```r
design <- model.matrix(~gata3+geomsizecpm, eSet_adj)
sets_indices <- ids2indices(gene_ids, toupper(rownames(eSet))) 
rescpm <- camera(eSet_adj, sets_indices, contrast = 2, design=design)
rescpm$Direction <- ifelse(rescpm$Direction=='Up', 1, -1)
rescpm$set <- row.names(rescpm)
setDT(rescpm)
setkey(rescpm, FDR)
rescpm[,totalrank:=seq_len(.N)]
```
The results do not change too much if `geomsizecpm` is not included.

Although no modules are FDR significant, the top few modules relate to immue signaling

```r
kable(rescpm[1:20,list(Direction, FDR, totalrank, set)]) 
```



| Direction|       FDR| totalrank|set                                                                       |
|---------:|---------:|---------:|:-------------------------------------------------------------------------|
|        -1| 0.7352039|         1|REACTOME_DOWNREGULATION_OF_ERBB2_ERBB3_SIGNALING                          |
|        -1| 0.7352039|         2|REACTOME_TRANSFERRIN_ENDOCYTOSIS_AND_RECYCLING                            |
|        -1| 0.7352039|         3|REACTOME_LATENT_INFECTION_OF_HOMO_SAPIENS_WITH_MYCOBACTERIUM_TUBERCULOSIS |
|         1| 0.7352039|         4|REACTOME_BETA_DEFENSINS                                                   |
|        -1| 0.7352039|         5|REACTOME_ENDOGENOUS_STEROLS                                               |
|        -1| 0.7352039|         6|REACTOME_IRON_UPTAKE_AND_TRANSPORT                                        |
|        -1| 0.7352039|         7|REACTOME_INSULIN_RECEPTOR_RECYCLING                                       |
|        -1| 0.7352039|         8|REACTOME_ACTIVATION_OF_CHAPERONES_BY_ATF6_ALPHA                           |
|         1| 0.7352039|         9|REACTOME_PHOSPHOLIPASE_C_MEDIATED_CASCADE                                 |
|         1| 0.7352039|        10|REACTOME_THE_NLRP3_INFLAMMASOME                                           |
|         1| 0.7352039|        11|REACTOME_INFLAMMASOMES                                                    |
|         1| 0.7352039|        12|REACTOME_TIE2_SIGNALING                                                   |
|         1| 0.7352039|        13|REACTOME_SIGNALING_BY_FGFR_MUTANTS                                        |
|        -1| 0.7352039|        14|REACTOME_REGULATION_OF_BETA_CELL_DEVELOPMENT                              |
|        -1| 0.7352039|        15|REACTOME_REGULATION_OF_GENE_EXPRESSION_IN_BETA_CELLS                      |
|         1| 0.7352039|        16|REACTOME_ION_TRANSPORT_BY_P_TYPE_ATPASES                                  |
|         1| 0.7352039|        17|REACTOME_ION_CHANNEL_TRANSPORT                                            |
|         1| 0.7352039|        18|REACTOME_CA_DEPENDENT_EVENTS                                              |
|         1| 0.7352039|        19|REACTOME_PLC_BETA_MEDIATED_EVENTS                                         |
|         1| 0.7352039|        20|REACTOME_DAG_AND_IP3_SIGNALING                                            |

