### R code from vignette source 'MIGSA.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: MIGSA.Rnw:133-136 (eval = FALSE)
###################################################
## ## try http:// if https:// URLs are not supported
## source("https://bioconductor.org/biocLite.R");
## biocLite("MIGSA");


###################################################
### code chunk number 2: MIGSA.Rnw:182-198
###################################################
library(GSEABase);
gs1 <- GeneSet(c("10", "1544", "1548", "1549", "1553", "7498", "9"), 
    setName="hsa00232", 
    setIdentifier="Caffeine metabolism");
gs1;

gs2 <- GeneSet(c("10229", "27235", "3242", "51004", "51805", "6898", "84274"),
    setName="hsa00130", 
    setIdentifier="Ubiquinone and other terpenoid-quinone biosynthesis");
gs3 <- GeneSet(c("11019", "387787", "51601"),
    setName="hsa00785", 
    setIdentifier="Lipoic acid metabolism");

## And now construct the GeneSetCollection object.
gsetsColl <- GeneSetCollection(list(gs1, gs2, gs3));
gsetsColl;


###################################################
### code chunk number 3: MIGSA.Rnw:211-221 (eval = FALSE)
###################################################
## ## Not run:
## 
## ## Load cellular component gene sets (another possibility would be "MF" or "BP")
## ccGsets <- loadGo("CC"); # It is a GeneSetCollection object
## 
## ## Load KEGG and Reactome gene sets
## keggReact <- downloadEnrichrGeneSets(c("KEGG_2015", "Reactome_2015"));
## ## It is a list object containing two GeneSetCollection objects
## 
## ## End(Not run)


###################################################
### code chunk number 4: Speedup1
###################################################
library(BiocParallel);
library(mGSZ);
library(MIGSA);
library(MIGSAdata);

data(tcgaMAdata);
subtypes <- tcgaMAdata$subtypes;
geneExpr <- tcgaMAdata$geneExpr;

## MA data: filter genes with less than 30% of genes read per condition
dim(geneExpr);
geneExpr <- geneExpr[
    rowSums(is.na(geneExpr[, subtypes == "Basal" ])) <
        .3*sum(subtypes == "Basal") &
    rowSums(is.na(geneExpr[, subtypes == "LumA" ])) <
        .3*sum(subtypes == "LumA")
    , ];
dim(geneExpr);


###################################################
### code chunk number 5: Speedup2 (eval = FALSE)
###################################################
## ## Not run:
## 
## ## Download GO and KEGG gene sets using MIGSA
## gSets <- list(
##             KEGG=downloadEnrichrGeneSets("KEGG_2015")[[1]],
##             BP=loadGo("BP"),
##             CC=loadGo("CC"),
##             MF=loadGo("MF"));
## gSetsList <- do.call(c, lapply(gSets, MIGSA:::asList));
## rm(gSets);
## 
## nCores <- c(1,2,4,8,10,12,14);
## allRes <- lapply(nCores, function(actCores) {
##     # setting in how many cores to run
##     bp_param <- MulticoreParam(workers=actCores, threshold="DEBUG",
##         progressbar=TRUE);
##     
##     set.seed(8818);
##     newtimeSpent <- Sys.time();
##     MIGSAmGSZres <- MIGSAmGSZ(geneExpr, gSetsList, subtypes, 
##         bp.param=bp_param);
##     newtimeSpent <- Sys.time()-newtimeSpent;
##     
##     res <- list(timeSpent=newtimeSpent, res=MIGSAmGSZres);
##     
##     return(res);
## })
## 
## set.seed(8818);
## timeSpent <- Sys.time();
## mGSZres <- mGSZ(geneExpr, gSetsList, subtypes);
## timeSpent <- Sys.time()-timeSpent;
## 
## mGSZres <- mGSZres$mGSZ;
## ## this tests that the returned values are equal, must give all TRUE
## lapply(allRes, function(actRes) {
##     actRes <- actRes$res;
##     actRes <- actRes[,1:4];
##     mergedRes <- merge(mGSZres, actRes, by="gene.sets",
##         suffixes=c("mGSZ", "MIGSAmGSZ"));
##     
##     all(unlist(lapply(2:4, function(x) {
##         all.equal(mergedRes[,x], mergedRes[,x+3])
##     })));
## })
## ## End(Not run)


###################################################
### code chunk number 6: Speedup3
###################################################
## As last chunk of code was not executed, we load that data:
library(MIGSAdata);

data(mGSZspeedup);
nCores <- mGSZspeedup$nCores;
allRes <- mGSZspeedup$allRes;
timeSpent <- mGSZspeedup$timeSpent;
## End(Loading data)

newtimeSpent <- lapply(allRes, function(actRes) {
    actRes$timeSpent;
})
names(newtimeSpent) <- nCores;

speeduptable <- c(timeSpent, unlist(newtimeSpent));
names(speeduptable) <- c(1, nCores);

## Let's put all times in the same unit in order to measure speedup
newtimeSpent <- lapply(newtimeSpent, function(acttime) {
    units(acttime) <- "secs";
    return(acttime);
});
units(timeSpent) <- "secs";

speedup <- do.call(c, lapply(newtimeSpent, function(acttime) 
    as.numeric(timeSpent)/as.numeric(acttime)));
speeduptable <- rbind(speeduptable, c(1, speedup));

## calculate efficiency
speeduptable <- rbind(speeduptable,
    speeduptable[2,] / as.numeric(colnames(speeduptable)));

rownames(speeduptable) <- c("Runtime", "Speedup", "Efficiency");
round(speeduptable, 2);


###################################################
### code chunk number 7: MIGSAmGSZ example
###################################################
library(MIGSA);

## Let's create our gene expression matrix with 200 genes and 8 subjects
nSamples <- 8; # 8 subjects
nGenes <- 200; # 200 genes
geneNames <- paste("g", 1:nGenes, sep = ""); # with names g1 ... g200

## Create random gene expression data matrix.
set.seed(8818);
exprMatrix <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);

## It must have rownames, as they will be treated as the gene names!
rownames(exprMatrix) <- geneNames;

## There will be 10 differentially expressed genes.
nDeGenes <- 10;
## Let's generate the offsets to sum to the differentially expressed genes.
deOffsets <- matrix(2*abs(rnorm(nDeGenes*nSamples/2)), ncol=nSamples/2);

## Randomly select which are the DE genes.
deIndexes <- sample(1:nGenes, nDeGenes, replace=FALSE);
exprMatrix[deIndexes, 1:(nSamples/2)] <-
    exprMatrix[deIndexes, 1:(nSamples/2)] + deOffsets;

## 4 subjects with condition C1 and 4 with C2.
conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));

nGSets <- 50; # 50 gene sets
## Let's create randomly 50 gene sets, of 10 genes each
gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
names(gSets) <- paste("set", as.character(1:nGSets), sep="");
## with names set1 ... set50

## And simply execute MIGSAmGSZ
MIGSAmGSZres <- MIGSAmGSZ(exprMatrix, gSets, conditions);

## It is just a simple data.frame
head(MIGSAmGSZres);


###################################################
### code chunk number 8: MIGSA example
###################################################
library(MIGSA);

## Let's simulate two expression matrices of 300 genes and 16 subjects.
nGenes <- 300; # 300 genes
nSamples <- 16; # 16 subjects
geneNames <- paste("g", 1:nGenes, sep = ""); # with names g1 ... g300

## Create the random gene expression data matrices.
set.seed(8818);
exprData1 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
rownames(exprData1) <- geneNames;
exprData2 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
rownames(exprData2) <- geneNames;

## There will be 30 differentially expressed genes.
nDeGenes <- nGenes/10;
## Let's generate the offsets to sum to the differentially expressed genes.
deOffsets <- matrix(2*abs(rnorm(nDeGenes*nSamples/2)), ncol=nSamples/2);

## Randomly select which are the DE genes.
deIndexes1 <- sample(1:nGenes, nDeGenes, replace=FALSE);
exprData1[deIndexes1, 1:(nSamples/2)] <-
    exprData1[deIndexes1, 1:(nSamples/2)] + deOffsets;

deIndexes2 <- sample(1:nGenes, nDeGenes, replace=FALSE);
exprData2[deIndexes2, 1:(nSamples/2)] <-
    exprData2[deIndexes2, 1:(nSamples/2)] + deOffsets;

exprData1 <- new("MAList",list(M=exprData1));
exprData2 <- new("MAList",list(M=exprData2));

## 8 subjects with condition C1 and 8 with C2.
conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
fitOpts <- FitOptions(conditions);

nGSets <- 30; # 30 gene sets
## Let's create randomly 30 gene sets, of 10 genes each

gSets1 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
names(gSets1) <- paste("set", as.character(1:nGSets), sep="");
myGSs1 <- as.Genesets(gSets1);

gSets2 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
names(gSets2) <- paste("set", as.character((nGSets+1):(2*nGSets)), sep="");
myGSs2 <- as.Genesets(gSets2);

igsaInput1 <- IGSAinput(name="igsaInput1", expr_data=exprData1, 
    fit_options=fitOpts);
igsaInput2 <- IGSAinput(name="igsaInput2", expr_data=exprData2, 
    fit_options=fitOpts);

experiments <- list(igsaInput1, igsaInput2);
## As we did not set gene sets for each IGSAinput, then we will have to 
## provide them in MIGSA function

## another way of generating the same MIGSA input would be setting the 
## gene sets individually to each IGSAinput:
igsaInput1 <- IGSAinput(name="igsaInput1", expr_data=exprData1, 
    fit_options=fitOpts, 
    gene_sets_list=list(myGeneSets1=myGSs1, myGeneSets2=myGSs2));
igsaInput2 <- IGSAinput(name="igsaInput2", expr_data=exprData2, 
    fit_options=fitOpts, 
    gene_sets_list=list(myGeneSets1=myGSs1, myGeneSets2=myGSs2));

experiments <- list(igsaInput1, igsaInput2);

## And then simply run MIGSA
migsaRes <- MIGSA(experiments);

## migsaRes contains the p-values obtained in each experiment for each gene set
head(migsaRes);


###################################################
### code chunk number 9: MIGSA.Rnw:494-516
###################################################
## Other possible analyses:
## If we want some gene sets to be evaluated in just one IGSAinput we 
## can do this:

## If we want to test myGSs1 in exprData1 and myGSs2 in exprData2:
igsaInput1 <- IGSAinput(name="igsaInput1", expr_data=exprData1, 
    fit_options=fitOpts, gene_sets_list=list(myGeneSets1=myGSs1));
igsaInput2 <- IGSAinput(name="igsaInput2", expr_data=exprData2, 
    fit_options=fitOpts, gene_sets_list=list(myGeneSets2=myGSs2));

experiments <- list(igsaInput1, igsaInput2);

## If we want to test myGSs1 in exprData1 and both in exprData2:
igsaInput1 <- IGSAinput(name="igsaInput1", expr_data=exprData1, 
    fit_options=fitOpts, gene_sets_list=list(myGeneSets1=myGSs1));
igsaInput2 <- IGSAinput(name="igsaInput2", expr_data=exprData2, 
    fit_options=fitOpts, 
    gene_sets_list=list(myGeneSets1=myGSs1, myGeneSets2=myGSs2));

experiments <- list(igsaInput1, igsaInput2);

## And this way, all possible combinations.


###################################################
### code chunk number 10: tcga MIGSA1
###################################################
library(edgeR);
library(limma);
library(MIGSA);
library(MIGSAdata);

data(tcgaMAdata);
data(tcgaRNAseqData);

geneExpr <- tcgaMAdata$geneExpr;
rnaSeq <- tcgaRNAseqData$rnaSeq;

subtypes <- tcgaMAdata$subtypes; # or tcgaRNAseqData$subtypes; are the same
fitOpts <- FitOptions(subtypes);

## MA data: filter genes with less than 30% of genes read per condition
dim(geneExpr);
geneExpr <- geneExpr[
    rowSums(is.na(geneExpr[, subtypes == "Basal" ])) <
        .3*sum(subtypes == "Basal") &
    rowSums(is.na(geneExpr[, subtypes == "LumA" ])) <
        .3*sum(subtypes == "LumA")
    , ];
dim(geneExpr);

## create our IGSAinput object
geneExpr <- new("MAList", list(M=geneExpr));
geneExprIgsaInput <- IGSAinput(
    name="tcgaMA",
    expr_data=geneExpr,
    fit_options=fitOpts,
    # with this treat we will get around 5% differentially expressed genes
    sea_params=SEAparams(treat_lfc=1.05));
summary(geneExprIgsaInput);


## RNAseq data: filter genes with less than 30% of genes read per 
## condition and (below)
dim(rnaSeq);
rnaSeq <- rnaSeq[
    rowSums(is.na(rnaSeq[, subtypes == "Basal" ])) <
        .3*sum(subtypes == "Basal") &
    rowSums(is.na(rnaSeq[, subtypes == "LumA" ])) <
        .3*sum(subtypes == "LumA")
    , ];
dim(rnaSeq);

## a mean less than 15 counts per condition.
rnaSeq <- rnaSeq[
    rowMeans(rnaSeq[, subtypes == "Basal" ], na.rm=TRUE) >= 15 &
    rowMeans(rnaSeq[, subtypes == "LumA"  ], na.rm=TRUE) >= 15
    , ];
dim(rnaSeq);

## create our IGSAinput object
rnaSeq <- DGEList(counts=rnaSeq);

rnaSeqIgsaInput <- IGSAinput(
    name="tcgaRNA",
    expr_data=rnaSeq,
    fit_options=fitOpts,
    # with this treat we will get around 5% differentially expressed genes
    sea_params=SEAparams(treat_lfc=1.45));
summary(rnaSeqIgsaInput);

experiments <- list(geneExprIgsaInput, rnaSeqIgsaInput);


###################################################
### code chunk number 11: tcga MIGSA2 (eval = FALSE)
###################################################
## ## Not run:
## 
## gSets <- list(
##             KEGG=downloadEnrichrGeneSets("KEGG_2015")[[1]],
##             BP=loadGo("BP"),
##             CC=loadGo("CC"),
##             MF=loadGo("MF"));
## 
## set.seed(8818);
## tcgaMigsaRes <- MIGSA(experiments, geneSets=gSets);
## 
## ## Time difference of 29.83318 mins in 10 cores
## ## End(Not run)


###################################################
### code chunk number 12: pbcmc MIGSA1
###################################################
library(limma);
library(MIGSA);
library(MIGSAdata);

data(pbcmcData);

## with these treat log fold change values we will get around 5% of 
## differentially expressed genes for each experiment
treatLfcs <- c(0.7, 0.2, 0.6, 0.25, 0.4, 0.75);
names(treatLfcs) <- c("mainz", "nki", "transbig", "unt", "upp", "vdx");

experiments <- lapply(names(treatLfcs), function(actName) {
    actData <- pbcmcData[[actName]];
    actExprs <- actData$geneExpr;
    actSubtypes <- actData$subtypes;
    
    # filtrate genes with less than 30% per condition
    actExprs <- actExprs[
        rowSums(is.na(actExprs[, actSubtypes == "Basal" ])) <
            .3*sum(actSubtypes == "Basal") &
        rowSums(is.na(actExprs[, actSubtypes == "LumA" ])) <
            .3*sum(actSubtypes == "LumA")
    , ]
    
    # create our IGSAinput object
    actExprData <- new("MAList", list(M=actExprs));
    actFitOpts <- FitOptions(actSubtypes);
    actIgsaInput <- IGSAinput(
        name=actName,
        expr_data=actExprData,
        fit_options=actFitOpts,
        sea_params=SEAparams(treat_lfc=treatLfcs[[actName]]));
    return(actIgsaInput);
})



###################################################
### code chunk number 13: pbcmc MIGSA2 (eval = FALSE)
###################################################
## ## Not run:
## 
## gSets <- list(
##             KEGG=downloadEnrichrGeneSets("KEGG_2015")[[1]],
##             BP=loadGo("BP"),
##             CC=loadGo("CC"),
##             MF=loadGo("MF"));
## 
## set.seed(8818);
## pbcmcMigsaRes <- MIGSA(experiments, geneSets=gSets);
## 
## ## Time difference of 1.26684 hours in 10 cores
## ## End(Not run)


###################################################
### code chunk number 14: MIGSAres merge1 (eval = FALSE)
###################################################
## ## Not run:
## 
## dim(pbcmcMigsaRes);
## # [1] 20425     9
## dim(tcgaMigsaRes);
## # [1] 20425     5
## 
## ## Let's merge both results in one big MIGSAres object
## bcMigsaRes <- merge(pbcmcMigsaRes, tcgaMigsaRes);
## dim(bcMigsaRes);
## # [1] 20425     11
## ## End(Not run)


###################################################
### code chunk number 15: MIGSA.Rnw:735-755
###################################################
## As last chunk of code was not executed, we load that data:
library(MIGSA);
library(MIGSAdata);
data(bcMigsaResAsList);
bcMigsaRes <- MIGSA:::MIGSAres.data.table(bcMigsaResAsList$dframe, 
bcMigsaResAsList$genesRank);
rm(bcMigsaResAsList);
## End(Loading data)

## Let's see a summary of enriched gene sets at different cutoff values
summary(bcMigsaRes);

## We will set a cutoff of 0.01 (recommended)
## A gene set will be considered enriched if its p-value is < 0.01 on 
## SEA or GSEA.
bcMigsaRes <- setEnrCutoff(bcMigsaRes, 0.01);

## The bcMigsaRes data object that is included in MIGSA package is the 
## following:
# bcMigsaRes <- bcMigsaRes[1:200,];


###################################################
### code chunk number 16: MIGSAres exploring1
###################################################
colnames(bcMigsaRes);
dim(bcMigsaRes);

summary(bcMigsaRes);

## We can see that 18,191 gene sets were not enriched, while 242 were 
## enriched in every dataset.
## Moreover, there is a high consensus between datasets, with a maximum of 679 
## enriched gene sets in common between upp and unt.
##
## Let's keep only gene sets enriched in at least one data set
bcMigsaRes <- bcMigsaRes[ rowSums(bcMigsaRes[,-(1:3)], na.rm=TRUE) > 0, ];
dim(bcMigsaRes);


###################################################
### code chunk number 17: MIGSA.Rnw:776-780
###################################################
# it is the same code as below, but saving the pdf to correctly load it in tex
pdf("bcMigsaResMigsaHeatmap.pdf");
aux <- migsaHeatmap(bcMigsaRes);
dev.off();


###################################################
### code chunk number 18: MIGSAres exploring2 (eval = FALSE)
###################################################
## ## Let's see enrichment heat map
## ## i.e. a heat map of binary data (enriched/not enriched)
## aux <- migsaHeatmap(bcMigsaRes);


###################################################
### code chunk number 19: MIGSAres exploring3
###################################################
## In this heat map we can see a high number of gene sets that are being 
## enriched in consensus by most of the datasets. Let's explore them.
## We can obtain them (enriched in at least 80% of datasets) by doing
consensusGsets <- bcMigsaRes[ rowSums(bcMigsaRes[, -(1:3)], na.rm=TRUE)
    > 6.4,];
dim(consensusGsets);

## And let's see from which sets are them
table(consensusGsets$GS_Name);


###################################################
### code chunk number 20: MIGSA.Rnw:803-808
###################################################
# it is the same code as below, but saving the pdf to correctly load it in tex
pdf("bcMigsaResGenesHeatmap.pdf");
aux <- genesHeatmap(bcMigsaRes, enrFilter=6.4, gsFilter=70,
    dendrogram="col");
dev.off();


###################################################
### code chunk number 21: MIGSAres exploring4 (eval = FALSE)
###################################################
## ## Moreover, let's see which are the genes that are mostly contributing 
## ## to gene set enrichment (genes contributing in at least 70 gene sets)
## ## i.e. a heat map showing the number of datasets in which each gene (columns) 
## ## contributed to enrich each gene set (rows).
## aux <- genesHeatmap(bcMigsaRes, enrFilter=6.4, gsFilter=70,
##     dendrogram="col");


###################################################
### code chunk number 22: MIGSAres exploring5
###################################################
## Well, we could continue exploring them, however, at the first heat map we
## can see that TCGA datasets are defining a separate cluster, this is caused 
## by a big group of gene sets that seem to be enriched mainly by TCGA.
## Let's explore them:
## (gene sets enriched by both TCGA datasets and in less than 20% of the other)
tcgaExclusive <- bcMigsaRes[
    rowSums(bcMigsaRes[, c("tcgaMA", "tcgaRNA")], na.rm=TRUE) == 2 & 
    rowSums(bcMigsaRes[, c("mainz","nki","transbig","unt","upp","vdx")],
        na.rm=TRUE) <  1.2
,];
dim(tcgaExclusive);

table(tcgaExclusive$GS_Name);

## Let's see which is this KEGG enriched gene set
tcgaExclusive[ tcgaExclusive$GS_Name == "KEGG_2015", "id" ];

## Let's see in which depths of the GO tree are these gene sets
table(getHeights(
    tcgaExclusive[ tcgaExclusive$GS_Name != "KEGG_2015", "id", drop=TRUE]));
## We can see that the most of the gene sets are between depths three and five


###################################################
### code chunk number 23: MIGSA.Rnw:846-850
###################################################
# it is the same code as below, but saving the pdf to correctly load it in tex
pdf("tcgaExclusiveMigsaGoTreeMF.pdf");
aux <- migsaGoTree(tcgaExclusive, ont="MF");
dev.off();


###################################################
### code chunk number 24: MIGSAres exploring6_1 (eval = FALSE)
###################################################
## ## And plot the GO tree of the other gene sets (except of CC, as it 
## ## has only three gene sets, and it will look bad)
## aux <- migsaGoTree(tcgaExclusive, ont="MF");


###################################################
### code chunk number 25: MIGSA.Rnw:861-865
###################################################
# it is the same code as below, but saving the pdf to correctly load it in tex
pdf("tcgaExclusiveMigsaGoTreeBP.pdf");
aux <- migsaGoTree(tcgaExclusive, ont="BP");
dev.off();


###################################################
### code chunk number 26: MIGSAres exploring6_2 (eval = FALSE)
###################################################
## aux <- migsaGoTree(tcgaExclusive, ont="BP");


###################################################
### code chunk number 27: MIGSA.Rnw:874-878
###################################################
# it is the same code as below, but saving the pdf to correctly load it in tex
pdf("tcgaExclusiveGenesBarplot.pdf");
mostEnrichedGenes <- genesBarplot(tcgaExclusive, gsFilter=12.45);
dev.off();


###################################################
### code chunk number 28: MIGSAres exploring6_3 (eval = FALSE)
###################################################
## ## Let's explore which are the genes that repeat the most in these 
## ## gene sets (that are present in at least 15% of the gene sets)
## ## i.e. a bar plot of the number of gene sets in which each gene contributed to 
## ## enrich.
## mostEnrichedGenes <- genesBarplot(tcgaExclusive, gsFilter=12.45);


###################################################
### code chunk number 29: MIGSAres exploring7
###################################################
mostEnrichedGenes$data;

## Gene 652 is contributing to enrichment in 15 gene sets. And in total 
## there are 6 genes that are being really active in TCGA enriched 
## gene sets
tcgaImportantGenes <- as.character(mostEnrichedGenes$data$id);


###################################################
### code chunk number 30: MIGSA.Rnw:900-904
###################################################
# it is the same code as below, but saving the pdf to correctly load it in tex
pdf("consensusGsetsGenesBarplot.pdf");
consMostEnrichedGenes <- genesBarplot(consensusGsets, gsFilter=53.25);
dev.off();


###################################################
### code chunk number 31: MIGSAres exploring8 (eval = FALSE)
###################################################
## ## Let's do the same analysis for the rest of the datasets, so we can filtrate 
## ## which genes are acting exclusively in TCGA datasets
## consMostEnrichedGenes <- genesBarplot(consensusGsets, gsFilter=53.25);


###################################################
### code chunk number 32: MIGSAres exploring9
###################################################
consImportantGenes <- as.character(consMostEnrichedGenes$data$id);

## Let's see which genes they share
intersect(tcgaImportantGenes, consImportantGenes);

## And get the really tcga exclusive genes (5 genes)
tcgaExclGenes <- setdiff(tcgaImportantGenes, consImportantGenes);


###################################################
### code chunk number 33: MIGSAres exploring10
###################################################
## Let's sample 4 genes from consImportantGenes (as if they are our interest 
## genes)
set.seed(8818);
myInterestGenes <- sample(consImportantGenes, 4);

## So we can get the filtered MIGSAres object by doing:
intGenesMigsa <- filterByGenes(bcMigsaRes, myInterestGenes);

dim(intGenesMigsa);

head(intGenesMigsa);


###################################################
### code chunk number 34: Session Info
###################################################
sessionInfo()


