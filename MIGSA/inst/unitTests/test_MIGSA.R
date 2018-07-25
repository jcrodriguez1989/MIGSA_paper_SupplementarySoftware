# library(covr);
library(RUnit);
library(MIGSA);

# MIGSA_coverage <- package_coverage("MIGSA"); MIGSA_coverage
# shine(MIGSA_coverage);
# BiocGenerics:::testPackage("MIGSA")

rOpts <- getOption("RUnit");
rOpts$silent <- !FALSE;
options("RUnit"=rOpts);

require(Biobase);
hasInternet <- testBioCConnection();

testAll <- FALSE;

require(BiocParallel);
if (.Platform$OS.type == "unix") {
    bp_param <- MulticoreParam(workers=1);
} else if (.Platform$OS.type == "windows") {
    bp_param <- SnowParam(workers=1);
}
register(bp_param);

######## FitOptions tests

###### FitOptions-class tests

# It doesnt have any function to test

###### FitOptions-constructor tests

#### Correct ones

test_FitOptions.default_ok <- function() {
    conditions <- c(rep("C1", 4), rep("C2", 7));
    
    checkTrue(validObject(FitOptions(conditions)));
}

test_FitOptions.data.frame_ok <- function() {
    myData <- data.frame(cond=c(rep("C1", 4), rep("C2", 7)));
    myFormula <- ~cond - 1;
    myContrast <- c(-1, 1);
    
    checkTrue(validObject(FitOptions(myData, myFormula, myContrast)));
}

test_FitOptions_bothConstructorEqual <- function() {
    conditions <- c(rep("C1", 4), rep("C2", 7));
    myData <- data.frame(cond=c(rep("C1", 4), rep("C2", 7)));
    myFormula <- ~cond - 1;
    myContrast <- c(-1, 1);
    
    fitDefault <- FitOptions(conditions);
    fitDataFrame <- FitOptions(myData, myFormula, myContrast);
    checkEquals(fitDefault, fitDataFrame);
}

#### Incorrect ones

test_FitOptions.default_wrong_oneCond <- function() {
    conditions <- rep("C1", 1);
    
    checkException(FitOptions(conditions));
}

test_FitOptions.default_wrong_threeCond <- function() {
    conditions <- c(rep("C1", 4), rep("C2", 7), rep("C3", 1));
    
    checkException(FitOptions(conditions));
}

###### FitOptions tests

# It doesnt have any exported function to test

###### Genesets-geneSetsFromFile tests

#### Correct ones

test_Genesets_geneSetsFromFile_ok_simpleFile <- function() {
    # create gene sets file
    genes <- paste("gene", 1:(10*20));
    gsets <- data.frame(
        IDs=paste("set", 1:10), 
        Names=rep("", 10), 
        matrix(genes, nrow=10, byrow=TRUE));
    
    geneSetsFile <- paste(tempdir(), "/fakeGsets.tsv", sep="");
    write.table(gsets, file=geneSetsFile, sep="\t",
        col.names=FALSE, row.names=FALSE, quote=FALSE);
    
    # Now lets load this tsv file as a Genesets object.
    myGsets <- geneSetsFromFile(geneSetsFile);
    myGsetsGO <- geneSetsFromFile(geneSetsFile, is_GO=TRUE);
    
    # Lets delete this tsv file
    unlink(geneSetsFile);
    
    require(GSEABase);
    
    checkTrue(all(unlist(lapply(myGsetsGO, function(x) 
        is(collectionType(x), "GOCollection")))));
    checkTrue(all(unlist(lapply(myGsets, function(x) 
        is(collectionType(x), "NullCollection")))));
    checkEquals(length(myGsets), 10);
    checkEquals(names(MIGSA:::asList(myGsets)), as.character(gsets$IDs));
    checkTrue(all(unlist(MIGSA:::asList(myGsets)) == genes));
}

#### Incorrect ones

test_Genesets_geneSetsFromFile_wrong_emptyPath <- function() {
    filePath <- "";
    name <- "myGenesets";
    
    checkException(
        geneSetsFromFile(filePath, name));
}

test_Genesets_geneSetsFromFile_wrong_dirPath <- function() {
    filePath <- list.dirs("..")[[2]];
    name <- "myGenesets";
    
    checkException(
        geneSetsFromFile(filePath, name));
}

test_Genesets_geneSetsFromFile_wrong_wrongPath <- function() {
    filePath <- "./notExistingPath/notExistingFile.csv";
    name <- "myGenesets";
    
    checkException(
        geneSetsFromFile(filePath, name));
}

###### Genesets-asGenesets tests

#### Correct ones

test_Genesets_asGenesets_ok_complete <- function() {
    require(GSEABase);
    
    myGs1 <- GeneSet(as.character(1:10), setName="fakeId1", setIdentifier="");
    myGs2 <- GeneSet(as.character(7:15), setName="fakeId2", setIdentifier="");
    myGs3 <- GeneSet(as.character(20:28), setName="fakeId3", setIdentifier="");
    
    constGsets <- GeneSetCollection(list(myGs1, myGs2, myGs3));
    
    listGsets <- list(geneIds(myGs1), geneIds(myGs2), geneIds(myGs3));
    names(listGsets) <- c(setName(myGs1), setName(myGs2), setName(myGs3));
    
    asGsets <- as.Genesets(listGsets);
    
    checkEquals(constGsets, asGsets);
}

test_Genesets_asGenesets_wrong_emptyGsetRemoved <- function() {
    require(GSEABase);
    myGs1 <- as.character(1:10);
    myGs2 <- "";
    myGs3 <- as.character(20:28);
    
    listGsets <- list(myGs1, myGs2, myGs3);
    names(listGsets) <- c("myGs1", "myGs2", "myGs3");
    
    asGsets <- as.Genesets(listGsets);
    checkEqualsNumeric(length(asGsets), 2);
    
    fstGset <- asGsets[[1]];
    sndGset <- asGsets[[2]];
    
    checkEquals(setName(fstGset), "myGs1");
    checkEquals(setName(sndGset), "myGs3"); # myGs2 must be deleted
    
    checkEquals(geneIds(fstGset), myGs1);
    checkEquals(geneIds(sndGset), myGs3); # myGs2 must be deleted
}

#### Incorrect ones

test_Genesets_asGenesets_wrong_noIds <- function() {
    myGs1 <- as.character(1:10);
    myGs2 <- as.character(7:15);
    
    listGsets <- list(myGs1, myGs2);
    
    checkException(
        as.Genesets(listGsets));
}

test_Genesets_asGenesets_wrong_repeatedIds <- function() {
    myGs1 <- as.character(1:10);
    myGs2 <- as.character(7:15);
    
    listGsets <- list(myGs1, myGs2);
    names(listGsets) <- c("myGs1", "myGs1");
    
    checkException(
        as.Genesets(listGsets));
}

test_Genesets_asGenesets_wrong_noGenes <- function() {
    myGs1 <- "";
    myGs2 <- "";
    
    listGsets <- list(myGs1, myGs2);
    names(listGsets) <- c("myGs1", "myGs2");
    
    checkException(
        as.Genesets(listGsets));
}

###### Genesets-enrichrGeneSets tests

test_Genesets_enrichrGeneSets_ok_wellListed <- function() {
    # if we dont have internet then dont fail the test
    require(Biobase);
    if (!hasInternet) {
        checkTrue(TRUE);
    }
    
    enrichrList <- enrichrGeneSets();
    
    checkEquals(class(enrichrList), "data.frame");
    checkTrue(ncol(enrichrList) > 0);
    checkTrue(nrow(enrichrList) > 0);
}

###### Genesets-downloadEnrichrGeneSets tests

#### Correct ones

if (testAll) {
    test_Genesets_downloadEnrichrGeneSets_ok_goCC <- function() {
        # if we dont have internet then dont fail the test
        require(Biobase);
        if (!hasInternet) {
            checkTrue(TRUE);
        }
        
        goCcName <- "GO_Cellular_Component_2013";
        goCc <- downloadEnrichrGeneSets(goCcName);
        
        checkEquals(length(goCc), 1);
        checkEquals(names(goCc), goCcName);
        
        checkEquals(names(goCc)[[1]], goCcName);
        goCc <- goCc[[1]];
#         checkTrue(goCc@is_GO);
        checkTrue(length(goCc) > 0);
    }
}

test_Genesets_downloadEnrichrGeneSets_ok_kegg <- function() {
    # if we dont have internet then dont fail the test
    require(Biobase);
    if (!hasInternet) {
        checkTrue(TRUE);
    }
    
    bioCartaName <- "BioCarta_2015";
    bioCarta <- downloadEnrichrGeneSets(bioCartaName);
    
    checkEquals(length(bioCarta), 1);
    checkEquals(names(bioCarta), bioCartaName);
    
    bioCarta <- bioCarta[[1]];
#         checkTrue(!kegg@is_GO);
    checkTrue(length(bioCarta) > 0);
}

if (testAll) {
    test_Genesets_downloadEnrichrGeneSets_ok_keggGoCC <- function() {
        # if we dont have internet then dont fail the test
        require(Biobase);
        if (!hasInternet) {
            checkTrue(TRUE);
        }
        
        libNames <- c("KEGG_2016", "GO_Cellular_Component_2013");
        libs <- downloadEnrichrGeneSets(libNames);
        
        checkEquals(length(libs), 2);
        checkEquals(names(libs), libNames);
        kegg <- libs[[1]];
        goCc <- libs[[2]];
        
        checkEquals(names(libs)[[1]], libNames[[1]]);
#         checkTrue(!kegg@is_GO);
        checkTrue(length(kegg) > 0);
        
        checkEquals(names(libs)[[2]], libNames[[2]]);
#         checkTrue(goCc@is_GO);
        checkTrue(length(goCc) > 0);
    }
}

if (testAll) {
    test_Genesets_downloadEnrichrGeneSets_ok_keggOneFakeLib <- function() {
        # if we dont have internet then dont fail the test
        require(Biobase);
        if (!hasInternet) {
            checkTrue(TRUE);
        }
        
        libNames <- c("KEGG_2013", "fakeLib");
        libs <- downloadEnrichrGeneSets(libNames);
        
        checkEquals(length(libs), 1);
        checkEquals(names(libs), libNames[[1]]);
        
        checkEquals(names(libs)[[1]], libNames[[1]]);
        kegg <- libs[[1]];
#         checkTrue(!kegg@is_GO);
        checkTrue(length(kegg) > 0);
    }
}

#### Incorrect ones

if (testAll) {
    test_Genesets_downloadEnrichrGeneSets_wrong_noLibs <- function() {
        # if we dont have internet then dont fail the test
        require(Biobase);
        if (!hasInternet) {
            checkTrue(TRUE);
        }
        
        libNames <- "";
        checkTrue(
            length(downloadEnrichrGeneSets(libNames))==0);
    }
}

if (testAll) {
    test_Genesets_downloadEnrichrGeneSets_wrong_fakeLibs <- function() {
        # if we dont have internet then dont fail the test
        require(Biobase);
        if (!hasInternet) {
            checkTrue(TRUE);
        }
        
        libNames <- c("fakeLib1", "fakeLib2");
        checkTrue(
            length(downloadEnrichrGeneSets(libNames))==0);
    }
}

###### Genesets-loadGo tests

#### Correct ones

if (testAll) {
    test_Genesets_loadGo_ok_cc <- function() {
        ccName <- "CC";
        goCc <- loadGo(ccName);
        
#         checkEquals(goCc@name, ccName);
#         checkTrue(goCc@is_GO);
        checkTrue(length(goCc) > 0);
    }
}

#### Incorrect ones

test_Genesets_loadGo_wrong_kegg <- function() {
    checkException(
        loadGo("KEGG"));
}

###### Genesets tests

# It doesnt have any exported function to test

######## GSEAparams tests

###### GSEAparams-class tests

#### Correct ones

test_GSEAparams_ok_default <- function() {
    checkTrue(validObject(GSEAparams()));
}

#### Incorrect ones

test_GSEAparams_wrong_1perms <- function() {
    checkException(
        GSEAparams(perm_number=1));
}

###### GSEAparams tests

# It doesnt have any exported function to test

######## SEAparams tests

###### SEAparams-class tests

#### Correct ones

test_SEAparams_ok_default <- function() {
    checkTrue(validObject(SEAparams()));
}

test_SEAparams_ok_bri <- function() {
    checkTrue(validObject(SEAparams(br="bri")));
}

test_SEAparams_ok_briii <- function() {
    checkTrue(validObject(SEAparams(br="briii")));
}

test_SEAparams_ok_ownbr <- function() {
    checkTrue(validObject(SEAparams(
                                br=as.character(1:10))));
}

test_SEAparams_ok_over1Lfc <- function() {
    checkTrue(validObject(SEAparams(
                                treat_lfc=1.1)));
}

#### Incorrect ones

test_SEAparams_wrong_negLfc <- function() {
    checkException(
        SEAparams(treat_lfc=-1));
}

test_SEAparams_wrong_negDe <- function() {
    checkException(
        SEAparams(de_cutoff=-1));
}

test_SEAparams_wrong_over1De <- function() {
    checkException(
        SEAparams(de_cutoff=1.1));
}

test_SEAparams_wrong_fakeAdj <- function() {
    checkException(
        SEAparams(adjust_method="fakeAdj"));
}

test_SEAparams_wrong_emptybr <- function() {
    checkException(
        SEAparams(br=""));
}

test_SEAparams_wrong_twoEmptybr <- function() {
    checkException(
        SEAparams(br=c("", "")));
}

###### SEAparams tests

# It doesnt have any exported function to test

######## GenesetRes tests

###### GenesetRes tests

# It doesnt have any exported function to test

######## GSEAres tests

###### GSEAres tests

# It doesnt have any exported function to test

######## SEAres tests

###### SEAres tests

# It doesnt have any exported function to test

######## GenesetsRes tests

###### GenesetsRes tests

# It doesnt have any exported function to test

######## IGSAinput tests

###### IGSAinput-class tests

#### Correct ones

test_IGSAinput_ok_complete <- function() {
    name <- "myIgsaInput";
    nSamples <- 4;
    nGenes <- 10;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    checkTrue(
        validObject(IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts
        )));
}

test_IGSAinput_ok_oneGSet <- function() {
    require(GSEABase);
    
    name <- "myIgsaInput";
    nSamples <- 4;
    nGenes <- 10;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    myGs <- GeneSet(as.character(1:10), setName="fakeName1", 
        setIdentifier="fakeId1");
    
    gsets <- list(myGenesets=
        GeneSetCollection(list(myGs)));
    
    checkTrue(
        validObject(IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts,
            gene_sets_list=gsets
        )));
}

test_IGSAinput_ok_twoGSets <- function() {
    require(GSEABase);
    
    name <- "myIgsaInput";
    nSamples <- 4;
    nGenes <- 10;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    gsName1 <- "myGenesets1";
    gsName2 <- "myGenesets2";
    myGs <- GeneSet(as.character(1:10), setName="fakeName1", 
        setIdentifier="fakeId1");
    
    gsets <- list(
        myGenesets1=GeneSetCollection(list(myGs)),
        myGenesets2=GeneSetCollection(list(myGs)));
    
    checkTrue(
        validObject(IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts,
            gene_sets_list=gsets
        )));
}

#### Incorrect ones

test_IGSAinput_wrong_noName <- function() {
    nSamples <- 4;
    nGenes <- 10;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    checkException(
        IGSAinput(
            expr_data=exprData,
            fit_options=fitOpts
        ));
}

test_IGSAinput_wrong_emptyName <- function() {
    name <- "";
    nSamples <- 4;
    nGenes <- 10;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    checkException(
        IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts
        ));
}

test_IGSAinput_wrong_oneSubj <- function() {
    name <- "myIgsaInput";
    nSamples <- 1;
    nGenes <- 10;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c("C1", "C2"));
    
    checkException(
        IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts
        ));
}

test_IGSAinput_wrong_oneGene <- function() {
    name <- "myIgsaInput";
    nSamples <- 4;
    nGenes <- 1;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    checkException(
        IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts
        ));
}

test_IGSAinput_wrong_badFitExprData <- function() {
    name <- "myIgsaInput";
    nSamples <- 4;
    nGenes <- 10;
    
    # one more subject in exprData than in fit
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples+1, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    checkException(
        IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts
        ));
}

test_IGSAinput_wrong_repGSets <- function() {
    require(GSEABase);
    name <- "myIgsaInput";
    nSamples <- 4;
    nGenes <- 10;
    exprData <- new("MAList", list(M=matrix(0, ncol=nSamples, nrow=nGenes)));
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    
    myGs <- GeneSet(as.character(1:10), setIdentifier="fakeId1", 
        setName="fakeName1");
    
    gsets <- list(
        myGenesets=GeneSetCollection(list(myGs)),
        myGenesets=GeneSetCollection(list(myGs)));
    
    checkException(
        IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts,
            gene_sets_list=gsets
        ));
}

###### IGSAinput-getterSetters tests

test_IGSAinput_getterSetters_ok <- function() {
    exprData1 <- new("MAList", list(M=matrix(0,ncol=6, nrow=6)));
    fitOpts <- FitOptions(c(1,1,1,2,2,2));
    igsaInput <- IGSAinput(
        name="fakeName",
        expr_data=exprData1,
        fit_options=fitOpts
    );
    
    checkEquals(name(igsaInput), "fakeName");
    name(igsaInput) <- "fakeName2";
    checkEquals(name(igsaInput), "fakeName2");
    
    checkEquals(fitOptions(igsaInput), fitOpts);
    fitOpts2 <- FitOptions(c(1,1,1,1,1,2));
    fitOptions(igsaInput) <- fitOpts2;
    checkEquals(fitOptions(igsaInput), fitOpts2);
    
    checkEquals(exprData(igsaInput), exprData1);
    exprData2 <- new("MAList", list(M=matrix(1,ncol=6, nrow=6)));
    exprData(igsaInput) <- exprData2;
    checkEquals(exprData(igsaInput), exprData2);
    
    checkEquals(length(MIGSA::geneSetsList(igsaInput)), 0);
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(as.character(1:100),
        size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    geneSetsList(igsaInput) <- list(myGeneSets=myGSs);
    checkEquals(length(MIGSA::geneSetsList(igsaInput)), 1);
    
    checkEquals(gseaParams(igsaInput), GSEAparams());
    gseaParams(igsaInput) <- GSEAparams(perm_number=10);
    checkEquals(gseaParams(igsaInput), GSEAparams(perm_number=10));
    
    checkEquals(seaParams(igsaInput), SEAparams());
    seaParams(igsaInput) <- SEAparams(treat_lfc=1);
    checkEquals(seaParams(igsaInput), SEAparams(treat_lfc=1));
}

###### IGSAinput-getDEGenes tests

test_IGSAinput_getDEGenes_ok_deGenes <- function() {
    name <- "myIgsaInput";
    nSamples <- 40;
    nGenes <- 100;
    
    set.seed(8818);
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    # lets make some DE genes
    exprData[1:5, 1:(nSamples/2)] <- 
        matrix(6*abs(rnorm(5*nSamples/2)), ncol=nSamples/2);
    rownames(exprData) <- 1:nrow(exprData);
    exprData <- new("MAList", list(M=exprData));
    
    fitOpts <- FitOptions(c(rep("C1", nSamples/2), rep("C2", nSamples/2)));
    igsaInput <- IGSAinput(
            name=name,
            expr_data=exprData,
            fit_options=fitOpts
        );
    
    newIgsaInput <- getDEGenes(igsaInput);
    deGenes <- seaParams(newIgsaInput)@de_genes;
    
    checkEquals(deGenes, as.character(1:5));
}

###### IGSAinput tests

# It doesnt have any exported function to test

######## MGSZ tests

###### MGSZ tests

# It doesnt have any exported function to test

######## MIGSAmGSZ tests

###### MIGSAmGSZ tests

# this check is failing in bioconductor servers
# test_MIGSAmGSZ_ok_sameAsMgsz <- function() {
#     require(mGSZ);
#     perms <- 4;
#     nGenes <- 100;
#     nSamples <- 6;
#     geneNames <- paste("g", 1:nGenes, sep = "");
#     
#     ## Create random gene expression data matrix.
#     set.seed(8818);
#     exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
#     rownames(exprData) <- geneNames;
#     
#     ## There will be nGenes/25 differentialy expressed genes.
#     nDeGenes <- nGenes/25;
#     ## Lets generate the offsets to sum to the differentialy expressed genes.
#     deOffsets <- matrix(2*abs(rnorm(nDeGenes*nSamples/2)), ncol=nSamples/2);
#     
#     ## Randomly select which are the DE genes.
#     deIndexes <- sample(1:nGenes, nDeGenes, replace=FALSE);
#     exprData[deIndexes, 1:(nSamples/2)] <-
#         exprData[deIndexes, 1:(nSamples/2)] + deOffsets;
#     
#     ## half of each condition.
#     conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
#     
#     nGSets <- 4;
#     gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
#     names(gSets) <- paste("set", as.character(1:nGSets), sep="");
#     
#     set.seed(8818);
#     mGSZres <- mGSZ(exprData, gSets, conditions, p=perms);
#     mGSZres <- mGSZres$mGSZ;
#     
#     set.seed(8818);
#     MIGSAmGSZres <- MIGSAmGSZ(exprData, gSets, conditions, p=perms);
#     
#     mergedRes <- merge(mGSZres, MIGSAmGSZres, by="gene.sets",
#     suffixes=c("mGSZ", "MIGSAmGSZ"))
#     
#     # this check is failing in bioconductor servers
#     checkTrue(all.equal(round(mergedRes$gene.set.scores, 5),
#         round(abs(mergedRes$mGszScore), 5)));
#     checkTrue(all.equal(mergedRes$pvaluemGSZ, mergedRes$pvalueMIGSAmGSZ));
# }

test_MIGSAmGSZ_ok_validWithVoom <- function() {
    perms <- 10;
    nGenes <- 500;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = ""); # with names g1 ... g1000
    
    ## Create random gene expression data matrix.
    set.seed(8818);
    exprData <- matrix(rnbinom(nGenes*nSamples, mu=5, size=2),ncol=nSamples);
    rownames(exprData) <- geneNames;
    
    ## There will be 40 differentialy expressed genes.
    nDeGenes <- nGenes/25;
    ## Lets generate the offsets to sum to the differentialy expressed genes.
    deOffsets <- matrix(2*abs(rnbinom(nDeGenes*nSamples/2, mu=5, size=2)),
                            ncol=nSamples/2);
    
    ## Randomly select which are the DE genes.
    deIndexes <- sample(1:nGenes, nDeGenes, replace=FALSE);
    exprData[deIndexes, 1:(nSamples/2)] <-
        exprData[deIndexes, 1:(nSamples/2)] + deOffsets;
    
    ## 15 subjects with condition C1 and 15 with C2.
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 10;
    ## Lets create randomly 200 gene sets, of 10 genes each
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    
    set.seed(8818);
    MIGSAmGSZres <- MIGSAmGSZ(exprData, gSets, conditions, use.voom=!F,
                            p=perms);
    
    checkEquals(ncol(MIGSAmGSZres), 4);
    checkEquals(nrow(MIGSAmGSZres), nGSets);
}

######## DEnricher tests

###### DEnricher tests

# It doesnt have any exported function to test

######## IGSAres tests

###### IGSAres tests

# It doesnt have any exported function to test

######## IGSA tests

###### IGSA tests

# It doesnt have any exported function to test

###### IGSAinput-common tests

test_IGSAinput_common_ok_summary <- function() {
    require(GSEABase);
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get some DE genes
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInputName <- "igsaInput";
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=GeneSetCollection(myGSs)));
    
    igsaSumm <- summary(igsaInput);
    checkTrue(all(igsaSumm == c("igsaInput","6","C1VSC2","3","3","1",
        "200","0","0.3","fdr","1","briii","5","0.5")));
}

######## MIGSA tests

###### MIGSA tests

#### Correct ones

test_MIGSA_ok_oneExp <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get some DE genes
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInputName <- "MIGSA_ok_oneExp";
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    
    checkEquals(ncol(migsaRes), 4);
    checkEquals(nrow(migsaRes), nGSets);
    
    checkTrue(length(unique(migsaRes$GS_Name))==1);
    checkTrue(unique(migsaRes$GS_Name) == "myGeneSets");
    
    checkTrue(all(names(gSets) %in% migsaRes$id));
    
    checkTrue(colnames(migsaRes)[[4]] == igsaInputName);
}

test_MIGSA_ok_twoExp <- function() {
    set.seed(8818);
    nGenes <- 100;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData1 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData1) <- geneNames;
    exprData1 <- new("MAList",list(M=exprData1));
    
    exprData2 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData2) <- geneNames;
    exprData2 <- new("MAList",list(M=exprData2));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 4;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get some DE genes
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInput1Name <- "MIGSA_ok_twoExp1";
    igsaInput1 <- IGSAinput(name=igsaInput1Name, expr_data=exprData1,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    
    igsaInput2Name <- "MIGSA_ok_twoExp2";
    igsaInput2 <- IGSAinput(name=igsaInput2Name, expr_data=exprData2,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    
    experiments <- list(igsaInput1, igsaInput2);
    
    migsaRes <- MIGSA(experiments);
    
    checkEquals(ncol(migsaRes), 5);
    checkEquals(nrow(migsaRes), nGSets);
    
    checkTrue(length(unique(migsaRes$GS_Name))==1);
    checkTrue(unique(migsaRes$GS_Name) == "myGeneSets");
    
    checkTrue(all(names(gSets) %in% migsaRes$id));
    
    checkTrue(colnames(migsaRes)[[4]] == igsaInput1Name);
    checkTrue(colnames(migsaRes)[[5]] == igsaInput2Name);
}

test_MIGSA_ok_twoGSs <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));

    nGSets <- 10;
    
    gSets1 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets1) <- paste("set", as.character(1:nGSets), sep="");
    myGSs1 <- as.Genesets(gSets1);
    
    gSets2 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets2) <- paste("set", as.character(1:nGSets), sep="");
    myGSs2 <- as.Genesets(gSets2);
    
    fitOpts <- FitOptions(conditions);
    
    # to get some DE genes
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInputName <- "MIGSA_ok_twoGSs";
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets1=myGSs1, myGeneSets2=myGSs2));
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    
    checkEquals(ncol(migsaRes), 4);
    checkEquals(nrow(migsaRes), 2*nGSets);
    
    checkTrue(length(unique(migsaRes$GS_Name))==2);
    checkTrue(unique(migsaRes$GS_Name)[[1]] == "myGeneSets1");
    checkTrue(unique(migsaRes$GS_Name)[[2]] == "myGeneSets2");
    
    checkTrue(all(names(gSets1) %in% migsaRes$id));
    checkTrue(all(names(gSets2) %in% migsaRes$id));
    
    checkTrue(colnames(migsaRes)[[4]] == igsaInputName);
}

test_MIGSA_ok_twoExptwoGSs <- function() {
    set.seed(8818);
    nGenes <- 100;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData1 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData1) <- geneNames;
    exprData1 <- new("MAList",list(M=exprData1));
    
    exprData2 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData2) <- geneNames;
    exprData2 <- new("MAList",list(M=exprData2));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 4;
    
    gSets1 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets1) <- paste("set", as.character(1:nGSets), sep="");
    myGSs1 <- as.Genesets(gSets1);
    
    gSets2 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets2) <- paste("set", as.character(1:nGSets), sep="");
    myGSs2 <- as.Genesets(gSets2);
    
    fitOpts <- FitOptions(conditions);
    
    # to get some DE genes
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInput1Name <- "MIGSA_ok_twoExptwoGSs1";
    igsaInput1 <- IGSAinput(name=igsaInput1Name, expr_data=exprData1,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets1=myGSs1, myGeneSets2=myGSs2));
    
    igsaInput2Name <- "MIGSA_ok_twoExptwoGSs2";
    igsaInput2 <- IGSAinput(name=igsaInput2Name, expr_data=exprData2,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets1=myGSs1, myGeneSets2=myGSs2));
    
    experiments <- list(igsaInput1, igsaInput2);
    
    migsaRes <- MIGSA(experiments);
    
    checkEquals(ncol(migsaRes), 5);
    checkEquals(nrow(migsaRes), 2*nGSets);
    
    checkTrue(length(unique(migsaRes$GS_Name))==2);
    checkTrue(unique(migsaRes$GS_Name)[[1]] == "myGeneSets1");
    checkTrue(unique(migsaRes$GS_Name)[[2]] == "myGeneSets2");
    
    checkTrue(all(names(gSets1) %in% migsaRes$id));
    checkTrue(all(names(gSets2) %in% migsaRes$id));
    
    checkTrue(colnames(migsaRes)[[4]] == igsaInput1Name);
    checkTrue(colnames(migsaRes)[[5]] == igsaInput2Name);
}

test_MIGSA_ok_twoExpFourGSs <- function() {
    set.seed(8818);
    nGenes <- 100;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData1 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData1) <- geneNames;
    exprData1 <- new("MAList",list(M=exprData1));
    
    exprData2 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData2) <- geneNames;
    exprData2 <- new("MAList",list(M=exprData2));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 4;
    
    gSets1 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets1) <- paste("set", as.character(1:nGSets), sep="");
    myGSs1 <- as.Genesets(gSets1);
    
    gSets2 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets2) <- paste("set", as.character(1:nGSets), sep="");
    myGSs2 <- as.Genesets(gSets2);
    
    gSets3 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets3) <- paste("set", as.character((nGSets+1):(2*nGSets)), sep="");
    myGSs3 <- as.Genesets(gSets3);
    
    gSets4 <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets4) <- paste("set", as.character((nGSets+1):(2*nGSets)), sep="");
    myGSs4 <- as.Genesets(gSets4);
    
    fitOpts <- FitOptions(conditions);
    
    # to get some DE genes
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInput1Name <- "MIGSA_ok_twoExpFourGSs1";
    igsaInput1 <- IGSAinput(name=igsaInput1Name, expr_data=exprData1,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets1=myGSs1, myGeneSets3=myGSs3));
    
    igsaInput2Name <- "MIGSA_ok_twoExpFourGSs2";
    igsaInput2 <- IGSAinput(name=igsaInput2Name, expr_data=exprData2,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets2=myGSs2, myGeneSets4=myGSs4));
    
    experiments <- list(igsaInput1, igsaInput2);
    
    migsaRes <- MIGSA(experiments);
    
    checkEquals(ncol(migsaRes), 5);
    checkEquals(nrow(migsaRes), 4*nGSets);
    
    checkTrue(length(unique(migsaRes$GS_Name))==4);
    checkTrue(unique(migsaRes$GS_Name)[[1]] == "myGeneSets1");
    checkTrue(unique(migsaRes$GS_Name)[[2]] == "myGeneSets2");
    checkTrue(unique(migsaRes$GS_Name)[[3]] == "myGeneSets3");
    checkTrue(unique(migsaRes$GS_Name)[[4]] == "myGeneSets4");
    
    checkTrue(all(names(gSets1) %in% migsaRes$id));
    checkTrue(all(names(gSets2) %in% migsaRes$id));
    checkTrue(all(names(gSets3) %in% migsaRes$id));
    checkTrue(all(names(gSets4) %in% migsaRes$id));
    
    checkTrue(colnames(migsaRes)[[4]] == igsaInput1Name);
    checkTrue(colnames(migsaRes)[[5]] == igsaInput2Name);
}

test_MIGSA_ok_noErrorWNoDE <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get 0 DE genes
    seaParams <- SEAparams(de_cutoff=0);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInputName <- "MIGSA_ok_noErrorWNoDE";
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    checkTrue(validObject(migsaRes));
    checkEquals(ncol(migsaRes), 4);
    checkEquals(nrow(migsaRes), nGSets);
    
    checkTrue(all(migsaRes@migsa_res_all$SEA_pval == 1));
    checkTrue(all(is.na(migsaRes@migsa_res_all$SEA_score)));
    checkTrue(all(migsaRes@migsa_res_all$SEA_enriching_genes == ""));
}

test_MIGSA_ok_noErrorWtwoPerm <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get at least one DE gene
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=2);
    
    igsaInputName <- "MIGSA_ok_noErrorWtwoPerm";
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    checkTrue(validObject(migsaRes));
    checkEquals(ncol(migsaRes), 4);
    checkEquals(nrow(migsaRes), nGSets);
}

test_MIGSA_ok_manuallyDEGenes <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # de_cutoff=1 so if we use fit then all genes are DE. but we just set as DE 
    # the ones from gSet #1
    seaParams <- SEAparams(de_cutoff=1, de_genes=gSets[[1]]);
    gseaParams <- GSEAparams(perm_number=2);
    
    igsaInputName <- "MIGSA_ok_manuallyDEGenes";
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    checkEquals(sort(unique(unlist(strsplit(
        migsaRes@migsa_res_all$SEA_enriching_genes, ", ")))),
        sort(gSets[[1]]));
}

test_MIGSA_ok_usebri <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    # generate genes that are not in our experiment (so we get them in gSets).
    geneNames <- paste("g", 1:(10*nGenes), sep = "");
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get at least one DE gene
    seaParams1 <- SEAparams(de_cutoff=0.5); # with briii
    seaParams2 <- SEAparams(de_cutoff=0.5, br="bri"); # with bri
    gseaParams <- GSEAparams(perm_number=2, min_sz=0);
    
    igsaInputName1 <- "MIGSA_ok_usebri1";
    igsaInputName2 <- "MIGSA_ok_usebri2";
    
    igsaInput1 <- IGSAinput(name=igsaInputName1, expr_data=exprData,
        sea_params=seaParams1, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments1 <- list(igsaInput1);
    
    igsaInput2 <- IGSAinput(name=igsaInputName2, expr_data=exprData,
        sea_params=seaParams2, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments2 <- list(igsaInput2);
    
    set.seed(8818);
    migsaRes1 <- MIGSA(experiments1);
    
    set.seed(8818);
    migsaRes2 <- MIGSA(experiments2);
    
    # genes ranks are equal
    checkTrue(all(migsaRes1@genes_rank[[1]] == migsaRes2@genes_rank[[1]]));
    
    migsaRes1All <- as.data.frame(migsaRes1@migsa_res_all);
    migsaRes2All <- as.data.frame(migsaRes2@migsa_res_all);
    
    # results except for SEA are equal
    checkEquals(migsaRes1All[,-(c(1, 6:10))], migsaRes2All[,-(c(1, 6:10))]);
    
    # results for SEA are different
    checkTrue(any(migsaRes1All[,8] != migsaRes2All[,8]));
}

test_MIGSA_ok_useOwnbr <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    # generate genes that are not in our experiment (so we get them in gSets).
    geneNames <- paste("g", 1:(10*nGenes), sep = "");
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get at least one DE gene
    seaParams1 <- SEAparams(de_cutoff=0.5); # with briii
    seaParams2 <- SEAparams(de_cutoff=0.5,
        br=geneNames[1:(length(geneNames)/2)]);
    gseaParams <- GSEAparams(perm_number=2, min_sz=0);
    
    igsaInputName1 <- "MIGSA_ok_useOwnbr1";
    igsaInputName2 <- "MIGSA_ok_useOwnbr2";
    
    igsaInput1 <- IGSAinput(name=igsaInputName1, expr_data=exprData,
        sea_params=seaParams1, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments1 <- list(igsaInput1);
    
    igsaInput2 <- IGSAinput(name=igsaInputName2, expr_data=exprData,
        sea_params=seaParams2, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments2 <- list(igsaInput2);
    
    set.seed(8818);
    migsaRes1 <- MIGSA(experiments1);
    
    set.seed(8818);
    migsaRes2 <- MIGSA(experiments2);
    
    # genes ranks are equal
    checkTrue(all(migsaRes1@genes_rank[[1]] == migsaRes2@genes_rank[[1]]));
    
    migsaRes1All <- as.data.frame(migsaRes1@migsa_res_all);
    migsaRes2All <- as.data.frame(migsaRes2@migsa_res_all);
    
    # results except for SEA are equal
    checkEquals(migsaRes1All[,-(c(1, 6:10))], migsaRes2All[,-(c(1, 6:10))]);
    
    # results for SEA are different
    checkTrue(any(migsaRes1All[,8] != migsaRes2All[,8]));
}

test_MIGSA_ok_useDifferenttests <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    # generate genes that are not in our experiment (so we get them in gSets).
    geneNames <- paste("g", 1:(10*nGenes), sep = "");
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get at least one DE gene
    seaParams1 <- SEAparams(de_cutoff=0.5);
    seaParams2 <- SEAparams(de_cutoff=0.5, test="HypergeoTest");
    seaParams3 <- SEAparams(de_cutoff=0.5, test="BinomialTest");
    gseaParams <- GSEAparams(perm_number=2, min_sz=0);
    
    igsaInputName1 <- "MIGSA_ok_useDifferenttests1";
    igsaInputName2 <- "MIGSA_ok_useDifferenttests2";
    igsaInputName3 <- "MIGSA_ok_useDifferenttests3";
    
    igsaInput1 <- IGSAinput(name=igsaInputName1, expr_data=exprData,
        sea_params=seaParams1, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments1 <- list(igsaInput1);
    
    igsaInput2 <- IGSAinput(name=igsaInputName2, expr_data=exprData,
        sea_params=seaParams2, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments2 <- list(igsaInput2);
    
    igsaInput3 <- IGSAinput(name=igsaInputName3, expr_data=exprData,
        sea_params=seaParams3, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    # to test also putting gene sets in MIGSAinput
    experiments3 <- list(igsaInput3);
    
    set.seed(8818);
    migsaRes1 <- MIGSA(experiments1);
    
    set.seed(8818);
    migsaRes2 <- MIGSA(experiments2);
    
    
    set.seed(8818);
    migsaRes3 <- MIGSA(experiments3, geneSets=list(myGeneSets=myGSs));
    
    # genes ranks are equal
    checkTrue(all(migsaRes1@genes_rank[[1]] == migsaRes2@genes_rank[[1]]));
    checkTrue(all(migsaRes1@genes_rank[[1]] == migsaRes3@genes_rank[[1]]));
    
    migsaRes1All <- as.data.frame(migsaRes1@migsa_res_all);
    migsaRes2All <- as.data.frame(migsaRes2@migsa_res_all);
    migsaRes3All <- as.data.frame(migsaRes3@migsa_res_all);
    
    # results except for SEA are equal
    checkEquals(migsaRes1All[,-(c(1, 6:10))], migsaRes2All[,-(c(1, 6:10))]);
    checkEquals(migsaRes1All[,-(c(1, 6:10))], migsaRes3All[,-(c(1, 6:10))]);
    
    # results for SEA are different
    checkTrue(any(migsaRes1All[,8] != migsaRes2All[,8]));
    checkTrue(any(migsaRes1All[,8] != migsaRes3All[,8]));
    checkTrue(any(migsaRes2All[,8] != migsaRes3All[,8]));
}

test_MIGSA_ok_onlySeaOrGsea <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get some DE genes
    seaParams <- SEAparams(de_cutoff=0.3);
    gseaParams <- GSEAparams(perm_number=5);
    
    igsaInput1 <- IGSAinput(name="MIGSA_ok_onlySeaOrGsea1", expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    set.seed(8818);
    migsaRes1 <- MIGSA(list(igsaInput1));
    
    igsaInput2 <- IGSAinput(name="MIGSA_ok_onlySeaOrGsea2", expr_data=exprData,
        sea_params=seaParams, gsea_params=NULL, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    set.seed(8818);
    migsaRes2 <- MIGSA(list(igsaInput2));
    
    # todo: uncomment this and pass tests on Windows i386
#     igsaInput3 <- IGSAinput(name="MIGSA_ok_onlySeaOrGsea3", expr_data=exprData,
#         sea_params=NULL, gsea_params=gseaParams, fit_options=fitOpts, 
#         gene_sets_list=list(myGeneSets=myGSs));
#     set.seed(8818);
#     migsaRes3 <- MIGSA(list(igsaInput3));
    
    # check results from both SEA are equal, and the rest NA
    checkEquals(migsaRes1@migsa_res_all[,2:10], migsaRes2@migsa_res_all[,2:10]);
    checkTrue(all(is.na(migsaRes2@migsa_res_all[,11:15])));
    
    # check results from both GSEA are equal, and the rest NA
#     checkEquals(migsaRes1@migsa_res_all[,c(2:5, 11:15)], 
#         migsaRes3@migsa_res_all[,c(2:5, 11:15)]);
#     checkTrue(all(is.na(migsaRes3@migsa_res_all[,6:10])));
}

test_MIGSA_ok_useOwnbrNoExprData <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    # generate genes that are not in our experiment (so we get them in gSets).
    geneNames <- paste("g", 1:(10*nGenes), sep = "");
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get at least one DE gene
    seaParams1 <- SEAparams(de_cutoff=0.5); # with briii
    
    igsaInputName1 <- "MIGSA_ok_useOwnbrNoExprData1";
    igsaInputName2 <- "MIGSA_ok_useOwnbrNoExprData2";
    
    igsaInput1 <- IGSAinput(name=igsaInputName1, expr_data=exprData,
        sea_params=seaParams1, gsea_params=NULL, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments1 <- list(igsaInput1);
    
    myDeGenes <- seaParams(getDEGenes(igsaInput1));
    igsaInput2 <- IGSAinput(name=igsaInputName2,
        sea_params=SEAparams(de_genes=myDeGenes@de_genes, 
            br=rownames(exprData)),
        gsea_params=NULL,
        gene_sets_list=list(myGeneSets=myGSs));
    experiments2 <- list(igsaInput2);
    
    set.seed(8818);
    migsaRes1 <- MIGSA(experiments1);
    
    set.seed(8818);
    migsaRes2 <- MIGSA(experiments2);
    
    checkEquals(migsaRes1@migsa_res_all[,-1], migsaRes2@migsa_res_all[,-1]);
    checkTrue(all(migsaRes1@migsa_res_summary == migsaRes2@migsa_res_summary));
}

#### Incorrect ones

test_MIGSA_wrong_noExperiments <- function() {
    checkException(MIGSA(list()));
}

test_MIGSA_wrong_wrongExperiments <- function() {
    checkException(MIGSA(list(1:10)));
}

test_MIGSA_wrong_repExpNames <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    fitOpts <- FitOptions(conditions);
    
    igsaInput <- IGSAinput(name="exp", expr_data=exprData,
        fit_options=fitOpts);
    
    checkException(MIGSA(list(
        igsaInput,
        igsaInput
    )));
}

test_MIGSA_wrong_noGSets <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    fitOpts <- FitOptions(conditions);
    
    igsaInputName <- "MIGSA_wrong_noGSets";
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        fit_options=fitOpts, gene_sets_list=list());
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    
    checkTrue(is.na(migsaRes));
}

# test_MIGSA_ok_wrongbr <- function() {
#     set.seed(8818);
#     nGenes <- 200;
#     nSamples <- 6;
#     geneNames <- paste("g", 1:nGenes, sep = "");
#     
#     exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
#     rownames(exprData) <- geneNames;
#     exprData <- new("MAList",list(M=exprData));
#     
#     conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
#     
#     # generate genes that are not in our experiment (so we get them in gSets).
#     geneNames <- paste("g", 1:(10*nGenes), sep = "");
#     nGSets <- 10;
#     gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
#     names(gSets) <- paste("set", as.character(1:nGSets), sep="");
#     myGSs <- as.Genesets(gSets);
#     
#     fitOpts <- FitOptions(conditions);
#     
#     # to get at least one DE gene
#     seaParams <- SEAparams(de_cutoff=0.5,
#         br=paste("gene", 1:100));
#     gseaParams <- GSEAparams(perm_number=2, min_sz=0);
#     
#     igsaInputName <- "igsaInput";
#     
#     igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
#         sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
#         gene_sets_list=list(myGeneSets=myGSs));
#     experiments <- list(igsaInput);
#     
#     migsaRes <- MIGSA(experiments);
#     checkTrue(is.na(migsaRes));
# }

test_MIGSA_ok_wrongbrOption <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    # generate genes that are not in our experiment (so we get them in gSets).
    geneNames <- paste("g", 1:(10*nGenes), sep = "");
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    # to get at least one DE gene
    seaParams <- SEAparams(de_cutoff=0.5, br="wrongbr");
    gseaParams <- GSEAparams(perm_number=2, min_sz=0);
    
    igsaInputName <- "MIGSA_ok_wrongbrOption";
    
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=seaParams, gsea_params=gseaParams, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    checkTrue(is.na(migsaRes));
}

test_MIGSA_ok_paramsCantbebothNull <- function() {
    set.seed(8818);
    nGenes <- 200;
    nSamples <- 6;
    geneNames <- paste("g", 1:nGenes, sep = "");
    
    exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
    rownames(exprData) <- geneNames;
    exprData <- new("MAList",list(M=exprData));
    
    conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
    
    # generate genes that are not in our experiment (so we get them in gSets).
    geneNames <- paste("g", 1:(10*nGenes), sep = "");
    nGSets <- 10;
    gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
    names(gSets) <- paste("set", as.character(1:nGSets), sep="");
    myGSs <- as.Genesets(gSets);
    
    fitOpts <- FitOptions(conditions);
    
    igsaInputName <- "MIGSA_ok_paramsCantbebothNull";
    
    igsaInput <- IGSAinput(name=igsaInputName, expr_data=exprData,
        sea_params=NULL, gsea_params=NULL, fit_options=fitOpts, 
        gene_sets_list=list(myGeneSets=myGSs));
    experiments <- list(igsaInput);
    
    migsaRes <- MIGSA(experiments);
    checkTrue(is.na(migsaRes));
}

######## MIGSAres tests

###### MIGSAres-class tests

# It doesnt have any exported function to test

###### MIGSAres-common tests

test_MIGSAres_common_ok_summaryMinimPval <- function() {
    data(migsaRes);
    
    migsaResShow <- show(migsaRes);
    migsaResPvals <- as.data.frame(migsaRes@migsa_res_all);
    
    checkEquals(rep(migsaResShow$id, 2), migsaResPvals$id);
    checkEquals(rep(migsaResShow$Name, 2), migsaResPvals$name);
    checkEquals(rep(migsaResShow$GS_Name, 2), migsaResPvals$gene_set_name);
    
    checkEquals(migsaResShow[,4], as.numeric(
        apply(migsaResPvals[1:nrow(migsaResShow), c(8,13)], 1, min)));
    checkEquals(migsaResShow[,5], as.numeric(
        apply(migsaResPvals[(nrow(migsaResShow)+1):nrow(migsaResPvals),
        c(8,13)], 1, min)));
}

test_MIGSAres_common_ok_tailHead <- function() {
    data(migsaRes);
    
    checkTrue(class(head(migsaRes)) == "MIGSAres");
    checkTrue(class(tail(migsaRes)) == "MIGSAres");
    
    checkEquals(nrow(head(migsaRes,10)), 10);
    checkEquals(nrow(tail(migsaRes,10)), 10);
    
    checkEquals(head(migsaRes,10), migsaRes[1:10,]);
    checkEquals(tail(migsaRes,10),
        migsaRes[(nrow(migsaRes)-9):nrow(migsaRes),]);
}

test_MIGSAres_common_ok_summary <- function() {
    data(migsaRes);
    
    migsaResSummary <- summary(migsaRes);
    checkEquals(c(migsaResSummary), c(6,15,31,8,23,37));
    migsaRes <- setEnrCutoff(migsaRes, 0.01);
    migsaResSummary <- summary(migsaRes);
    
    checkEquals(mean(migsaResSummary$consensusGeneSets), 100);
    checkEquals(unlist(c(migsaResSummary$enrichmentIntersections)),
        c(6,0,0,8));
}

test_MIGSAres_common_ok_asDframe <- function() {
    data(migsaRes);
    migsaResDframe <- as.data.frame(migsaRes);
    checkEquals(dim(migsaResDframe), c(200,5));
}

test_MIGSAres_common_ok_merge <- function() {
    data(migsaRes);
    data(bcMigsaRes);
    migsaResMerge <- merge(migsaRes, bcMigsaRes);
    
    checkEquals(nrow(migsaResMerge), nrow(migsaRes)+nrow(bcMigsaRes));
    checkEquals(ncol(migsaResMerge), ncol(migsaRes)+ncol(bcMigsaRes)-3);
    
    migsaResMergeDframe <- as.data.frame(migsaResMerge);
    migsaResDframe <- as.data.frame(migsaRes);
    bcMigsaResDframe <- as.data.frame(bcMigsaRes);
    
    checkTrue(all(migsaResDframe$id %in% migsaResMergeDframe$id));
    checkTrue(all(bcMigsaResDframe$id %in% migsaResMergeDframe$id));
    
    checkTrue(all(migsaResDframe$Name %in% migsaResMergeDframe$Name));
    checkTrue(all(bcMigsaResDframe$Name %in% migsaResMergeDframe$Name));
    
    checkTrue(all(migsaResDframe$GS_Name %in% migsaResMergeDframe$GS_Name));
    checkTrue(all(bcMigsaResDframe$GS_Name %in% migsaResMergeDframe$GS_Name));
    
    checkEquals(sort(migsaResMergeDframe[,4]),
        sort(migsaResDframe[,4]));
    checkEquals(sort(migsaResMergeDframe[,5]),
        sort(migsaResDframe[,5]));
    
    checkEquals(sort(migsaResMergeDframe[,6]),
        sort(bcMigsaResDframe[,4]));
    checkEquals(sort(migsaResMergeDframe[,7]),
        sort(bcMigsaResDframe[,5]));
    checkEquals(sort(migsaResMergeDframe[,8]),
        sort(bcMigsaResDframe[,6]));
    checkEquals(sort(migsaResMergeDframe[,9]),
        sort(bcMigsaResDframe[,7]));
    checkEquals(sort(migsaResMergeDframe[,10]),
        sort(bcMigsaResDframe[,8]));
    checkEquals(sort(migsaResMergeDframe[,11]),
        sort(bcMigsaResDframe[,9]));
    checkEquals(sort(migsaResMergeDframe[,12]),
        sort(bcMigsaResDframe[,10]));
    checkEquals(sort(migsaResMergeDframe[,13]),
        sort(bcMigsaResDframe[,11]));
}

test_MIGSAres_common_wrong_merge <- function() {
    data(migsaRes);
    
    # experiment names can not be the same
    checkException(merge(migsaRes, migsaRes));
}

###### MIGSAres-setEnrCutoff tests

test_MIGSAres_setEnrCutoff_ok <- function() {
    data(migsaRes);
    
    pvals <- migsaRes[,4:5];
    
    checkTrue(all(unlist(lapply(seq(0,1,by=0.05), function(actCoff) {
        actMigsaRes <- setEnrCutoff(migsaRes, actCoff);
        actEnr <- actMigsaRes[,4:5];
        all(actEnr == (pvals < actCoff));
    }))));
}

test_MIGSAres_setEnrCutoff_ok_cOffNoEnr <- function() {
    data(migsaRes);
    
    noEnrMigsaRes <- setEnrCutoff(migsaRes, 0);
    checkTrue(all(!noEnrMigsaRes[,4:5]));
}

test_MIGSAres_setEnrCutoff_ok_cOffAllEnr <- function() {
    data(migsaRes);
    
    allEnrMigsaRes <- setEnrCutoff(migsaRes, 1);
    checkTrue(all(!!allEnrMigsaRes[,4:5]));
}

test_MIGSAres_setEnrCutoff_ok_NACoff <- function() {
    data(migsaRes);
    
    checkEquals(class(migsaRes[,4, drop=!FALSE]), "numeric");
    checkEquals(class(migsaRes[,5, drop=!FALSE]), "numeric");
    
    migsaResCoff <- setEnrCutoff(migsaRes, 0.01);
    
    checkEquals(class(migsaResCoff[,4, drop=!FALSE]), "logical");
    checkEquals(class(migsaResCoff[,5, drop=!FALSE]), "logical");
}

###### MIGSAres-genesInSets tests

test_MIGSAres_genesInSets_ok <- function() {
    data(migsaRes);
    migsaRes <- migsaRes[1:10,];
    
    additionalInfo <- getAdditionalInfo(migsaRes);
    genesInfo <- additionalInfo[,
        c("experiment_name", "gene_set_name", "id", "SEA_pval", "GSEA_pval",
        "SEA_enriching_genes", "GSEA_enriching_genes")];
    
    checkTrue(all(unlist(lapply(seq(0,1,by=0.25), function(actCoff) {
        actMigsaRes <- setEnrCutoff(migsaRes, actCoff);
        actGInSets <- genesInSets(actMigsaRes);
        
        totalRes <- apply(unique(genesInfo[,c("gene_set_name", "id")]), 1, 
            function(actGSet) {
            actGSetInfo <- genesInfo[
                            genesInfo$gene_set_name == actGSet[[1]] &
                            genesInfo$id == actGSet[[2]]
                            ,]
            exprGenes <- table(unlist(apply(actGSetInfo, 1, 
                function(actActGSetInfo) {
                enrGenes <- NA;
                if (actActGSetInfo[["SEA_pval"]] < actCoff) {
                    enrGenes <- union(enrGenes, unlist(
                        strsplit(actActGSetInfo[["SEA_enriching_genes"]], ", ")
                        ));
                }
                if (actActGSetInfo[["GSEA_pval"]] < actCoff) {
                    enrGenes <- union(enrGenes, unlist(
                        strsplit(actActGSetInfo[["GSEA_enriching_genes"]], ", ")
                        ));
                }
                
                enrGenes <- enrGenes[!is.na(enrGenes)];
                return(enrGenes)
            })));
            
            actGISEprGenes <- actGInSets[
                    paste(actGSet[[1]], actGSet[[2]], sep="_"), ];
            res <- all(actGISEprGenes[names(exprGenes)] == exprGenes);
            res <- res && sum(actGISEprGenes) == sum(exprGenes);
            return(res);
        })
        return(totalRes);
    }))));
}

test_MIGSAres_genesInSets_ok_noEnriched <- function() {
    data(migsaRes);
    
    # So I get 0 enriched gene sets
    noEnrMigsaRes <- setEnrCutoff(migsaRes, 0);
    gInSets <- genesInSets(noEnrMigsaRes);
    
    checkEquals(nrow(gInSets), 200);
    checkEquals(ncol(gInSets), 0);
}

###### MIGSAres-filterByGenes tests

#### Correct ones

test_MIGSAres_filterByGenes_ok_someGenes <- function() {
    data(migsaRes);
    migsaRes <- migsaRes[1:40,];
    
    # to enrich every gene set
    migsaRes <- setEnrCutoff(migsaRes, 1);
    gInSets <- genesInSets(migsaRes);
    impGenes <- colnames(gInSets)[1:10];
    
    filtMigsaRes <- filterByGenes(migsaRes, impGenes);
    checkEquals(ncol(filtMigsaRes), 5);
}

test_MIGSAres_filterByGenes_ok_oneGene <- function() {
    data(migsaRes);
    migsaRes <- migsaRes[1:40,];
    
    # to enrich every gene set
    migsaRes <- setEnrCutoff(migsaRes, 1);
    gInSets <- genesInSets(migsaRes);
    impGenes <- colnames(gInSets)[[1]];
    
    filtMigsaRes <- filterByGenes(migsaRes, impGenes);
    checkEquals(ncol(filtMigsaRes), 5);
}

test_MIGSAres_filterByGenes_ok_emptyGenes <- function() {
    data(migsaRes);
    migsaRes <- migsaRes[1:40,];
    
    # to enrich every gene set
    migsaRes <- setEnrCutoff(migsaRes, 1);
    
    filtMigsaRes <- filterByGenes(migsaRes, "");
    checkEquals(ncol(filtMigsaRes), 0);
    checkEquals(nrow(filtMigsaRes), 0);
}

test_MIGSAres_filterByGenes_ok_fakeGenes <- function() {
    data(migsaRes);
    
    # to enrich every gene set
    migsaRes <- setEnrCutoff(migsaRes, 1);
    
    filtMigsaRes <- filterByGenes(migsaRes, c("fakeGene1", "fakeGene2"));
    checkEquals(ncol(filtMigsaRes), 0);
    checkEquals(nrow(filtMigsaRes), 0);
}

test_MIGSAres_filterByGenes_ok_noEnriched <- function() {
    data(migsaRes);
    
    # to enrich every gene set
    migsaRes <- setEnrCutoff(migsaRes, 1);
    gInSets <- genesInSets(migsaRes);
    impGenes <- colnames(gInSets);
    
    # to enrich 0 gene sets
    migsaRes <- setEnrCutoff(migsaRes, 0);
    
    filtMigsaRes <- filterByGenes(migsaRes, impGenes);
    checkEquals(ncol(filtMigsaRes), 0);
    checkEquals(nrow(filtMigsaRes), 0);
}

###### MIGSAres-genesBarplot tests

test_MIGSAres_genesBarplot_ok_simplePlot <- function() {
    data(migsaRes);
    
    plotRes <- genesBarplot(migsaRes);
    checkTrue(is(plotRes, "gg"));
#     checkEqualsNumeric(summary(plotRes$data$number), c(1,1,1,1.213,1,4));
    checkEqualsNumeric(length(unique(plotRes$data$id)), 47);
}

###### MIGSAres-genesHeatmap tests

test_MIGSAres_genesHeatmap_ok_simplePlot <- function() {
    data(migsaRes);
    
    plotRes <- genesHeatmap(migsaRes);
    checkTrue(is(plotRes, "list"));
    checkEqualsNumeric(plotRes$rowInd, 
        c(13, 8, 5, 6, 2, 12, 11, 4, 14, 1, 10, 9, 3, 7));
    checkEqualsNumeric(summary(plotRes$colInd), c(1,12.5,24,24,35.5,47));
#     checkEqualsNumeric(summary(c(plotRes$data)), c(0,0,0,0.08662614,0,1));
}

###### MIGSAres-geneSetBarplot tests

test_MIGSAres_geneSetBarplot_ok_simplePlot <- function() {
    data(migsaRes);
    
    plotRes <- geneSetBarplot(migsaRes);
    checkTrue(is(plotRes, "gg"));
    checkEqualsNumeric(summary(plotRes$data$number), c(1,1,1,1,1,1));
    checkEqualsNumeric(length(unique(plotRes$data$id)), 14);
    checkEqualsNumeric(length(unique(plotRes$data$GS_Name)), 1);
}

###### MIGSAres-getAdditionalInfo tests

# quite a irrelevant test, as the function does the same.
test_MIGSAres_getAdditionalInfo_ok <- function() {
    data(migsaRes);
    
    additionalInfo <- getAdditionalInfo(migsaRes);
    allInfo <- migsaRes@migsa_res_all;
    
    checkEquals(ncol(allInfo), ncol(additionalInfo));
    
    # to avoid NAs
    allInfo[is.na(allInfo)] <- 0;
    additionalInfo[is.na(additionalInfo)] <- 0;
    
    checkTrue(all(unlist(lapply(1:ncol(allInfo), function(j) {
        all(allInfo[[j]] == additionalInfo[,j]);
    }))));
}

###### MIGSAres-migsaHeatmap tests

test_MIGSAres_migsaHeatmap_ok_simplePlot <- function() {
    data(migsaRes);
    migsaRes <- migsaRes[1:50, ];
    migsaRes <- setEnrCutoff(migsaRes, 0.01);
    plotRes <- migsaHeatmap(migsaRes);
    
    checkTrue(is(plotRes, "list"));
    checkEqualsNumeric(plotRes$rowInd, c(1,3,2,4));
    checkEqualsNumeric(plotRes$colInd, 1:2);
#     checkEqualsNumeric(summary(c(plotRes$data)), c(0,0,0.5,0.5,1,1));
}

###### MIGSAres tests

# It doesnt have any exported function to test

######## GoAnalysis tests

###### GoAnalysis-getHeights tests

if (testAll) {
    test_GoAnalysis_getHeights <- function() {
        heights <- getHeights(
            c("GO:0008150", "GO:0007610", "GO:0050789", "fakeId"));
        checkEquals(heights[1:3], c(0,1,1));
        checkTrue(is.na(heights[[4]]));
    }
}

test_GoAnalysis_getHeights_maxHeights <- function() {
    heights <- getHeights(
        c("GO:0008150", "GO:0007610", "GO:0050789", "fakeId"),
        minHeight=FALSE);
    checkEquals(heights[1:3], c(0,1,2));
    checkTrue(is.na(heights[[4]]));
}

###### GoAnalysis-migsaGoTree tests

#### Correct ones

test_GoAnalysis_migsaGoTree_ok_simplePlot <- function() {
    data(bcMigsaRes);
    bcMigsaRes <- bcMigsaRes[5:6,];
    
    plotRes <- migsaGoTree(bcMigsaRes, ont="MF");
    checkTrue(all(sort(unlist(lapply(plotRes$gotree, nrow))) == 
        sort(table(bcMigsaRes$GS_Name)[-3])));
}

#### Incorrect ones

test_GoAnalysis_migsaGoTree_wrong_noGO <- function() {
    data(bcMigsaRes);
    keggMigsaRes <- bcMigsaRes[ bcMigsaRes$GS_Name == "KEGG_2015", ]
    
    checkException(migsaGoTree(keggMigsaRes));
}
