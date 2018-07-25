### R code from vignette source 'gettingPbcmcData.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: gettingData
###################################################
library(limma);
library(pbcmc);

# datasets included in BioConductor repository
libNames <- c("mainz", "nki", "transbig", "unt", "upp", "vdx");

# let's load them!
pbcmcData <- lapply(libNames, function(actLibName) {
    print(actLibName);
    
    # the pbcmc package provides an easy way to download and classify them
    actLib <- loadBCDataset(Class=PAM50, libname=actLibName, verbose=FALSE);
    actLibFilt <- filtrate(actLib, verbose=FALSE);
    actLibFilt <- classify(actLibFilt, std="none", verbose=FALSE);
    actSubtypes <- classification(actLibFilt)$subtype;
    
    # get the expression matrix and the annotation
    actExprs <- exprs(actLib);
    actAnnot <- annotation(actLib);
    
    # we recommend working allways with Entrez IDs, let's match them with 
    # expression matrix rownames (and modify them)
    if (all(actAnnot$probe == rownames(actExprs))) {
        actExprs <- actExprs[!is.na(actAnnot$EntrezGene.ID),];
        actAnnot <- actAnnot[!is.na(actAnnot$EntrezGene.ID),];
        rownames(actExprs) <- as.character(actAnnot$EntrezGene.ID);
    } else {
        matchedEntrez <- match(rownames(actExprs), actAnnot$probe);
        # all(rownames(actExprs) %in% actAnnot$probe == !is.na(matchedEntrez));
        
        stopifnot(all(
            actAnnot$probe[!is.na(matchedEntrez)] ==
            rownames(actExprs)[!is.na(matchedEntrez)]));
        
        actExprs <- actExprs[!is.na(matchedEntrez),];
        actAnnot <- actAnnot[!is.na(matchedEntrez),];
        stopifnot(all(actAnnot$probe == rownames(actExprs)));
        actExprs <- actExprs[!is.na(actAnnot$EntrezGene.ID),];
        actAnnot <- actAnnot[!is.na(actAnnot$EntrezGene.ID),];
        rownames(actExprs) <- as.character(actAnnot$EntrezGene.ID);
    }
    
    # average repeated genes expression
    actExprs <- avereps(actExprs);
    
    stopifnot(all(colnames(actExprs) == names(actSubtypes)));
    # filtrate only these two conditions
    actExprs <- actExprs[, actSubtypes %in% c("Basal", "LumA")];
    actSubtypes <- as.character(
        actSubtypes[actSubtypes %in% c("Basal", "LumA")]);
    
    return(list(geneExpr=actExprs, subtypes=actSubtypes));
})
names(pbcmcData) <- libNames;


###################################################
### code chunk number 2: validateData
###################################################
# save the just created pbcmcData to newPbcmcData
newPbcmcData <- pbcmcData;

library(MIGSAdata);

# and load the MIGSAdata one.
data(pbcmcData);
all.equal(newPbcmcData, pbcmcData);


###################################################
### code chunk number 3: Session Info
###################################################
sessionInfo()


