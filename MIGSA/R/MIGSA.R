#'MIGSA execution
#'
#'\code{MIGSA} runs a MIGSA execution. Functional analysis is done for each 
#'experiment by means of dEnricher and mGSZ.
#'
#'@param igsaInputs list of IGSAinput objects to execute.
#'@param geneSets (optional) named list of GeneSetCollection objects to be 
#'tested for enrichment (names must be unique). If provided then it will be 
#'tested in every IGSAinput, if not, each IGSAinput object must have its own 
#'list of GeneSetCollection.
#'@param ... not in use.
#'
#'@return A MIGSAres object.
#'
#'@docType methods
#'@name MIGSA
#'@rdname MIGSA
#'@seealso \code{\link{IGSAinput-class}}
#'@seealso \code{\link{MIGSAres-class}}
#'
#'@exportMethod MIGSA
#'
setGeneric(name="MIGSA", def=function(igsaInputs, ...) {
    standardGeneric("MIGSA")
})

#'@inheritParams MIGSA
#'@rdname MIGSA
#'@aliases MIGSA,list
#'
#'@importFrom futile.logger flog.info flog.warn
#'@include IGSA.R
#'@include IGSAinput.R
#'@include IGSAinput-class.R
#'@include IGSAinput-getterSetters.R
#'@include IGSAres.R
#'@include MIGSAres-class.R
#'@examples
#'## Lets simulate two expression matrices of 1000 genes and 30 subjects.
#'nGenes <- 1000; # 1000 genes
#'nSamples <- 30; # 30 subjects
#'geneNames <- paste("g", 1:nGenes, sep = ""); # with names g1 ... g1000
#'## Create random gene expression data matrix.
#'set.seed(8818);
#'exprData1 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
#'rownames(exprData1) <- geneNames;
#'exprData2 <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
#'rownames(exprData2) <- geneNames;
#'
#'## There will be 40 differentialy expressed genes.
#'nDeGenes <- nGenes/25;
#'## Lets generate the offsets to sum to the differentialy expressed genes.
#'deOffsets <- matrix(2*abs(rnorm(nDeGenes*nSamples/2)), ncol=nSamples/2);
#'
#'## Randomly select which are the DE genes.
#'deIndexes1 <- sample(1:nGenes, nDeGenes, replace=FALSE);
#'exprData1[deIndexes1, 1:(nSamples/2)] <-
#'exprData1[deIndexes1, 1:(nSamples/2)] + deOffsets;
#'
#'deIndexes2 <- sample(1:nGenes, nDeGenes, replace=FALSE);
#'exprData2[deIndexes2, 1:(nSamples/2)] <-
#'exprData2[deIndexes2, 1:(nSamples/2)] + deOffsets;
#'
#'exprData1 <- new("MAList",list(M=exprData1));
#'exprData2 <- new("MAList",list(M=exprData2));
#'
#'## 15 subjects with condition C1 and 15 with C2.
#'conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
#'
#'nGSets <- 200; # 200 gene sets
#'## Lets create randomly 200 gene sets, of 10 genes each
#'gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
#'names(gSets) <- paste("set", as.character(1:nGSets), sep="");
#'myGSs <- as.Genesets(gSets);
#'
#'fitOpts <- FitOptions(conditions);
#'
#'igsaInput1 <- IGSAinput(name="igsaInput1", expr_data=exprData1, 
#'fit_options=fitOpts, gene_sets_list=list(myGSs=myGSs));
#'igsaInput2 <- IGSAinput(name="igsaInput2", expr_data=exprData2, 
#'fit_options=fitOpts, gene_sets_list=list(myGSs=myGSs));
#'
#'experiments <- list(igsaInput1, igsaInput2);
#'
#'## Finally run MIGSA!
#'\dontrun{
#'migsaRes <- MIGSA(experiments);
#'}
#'
setMethod(
    f="MIGSA",
    signature=c("list"),
    definition=function(igsaInputs, geneSets=list()) {
        flog.info("*************************************");
        flog.info("Starting MIGSA analysis.");
        
        if (length(igsaInputs) <= 0) {
            stop("At least one IGSAinput must be provided");
        }
        
        if (!all(unlist(lapply(igsaInputs, function(x) is(x, "IGSAinput"))))) {
            stop("igsaInputs must be a list of IGSAinput objects");
        }
        
        # IGSAinput names must be unique
        exprs_names <- unlist(lapply(igsaInputs, name));
        if (length(exprs_names) != length(unique(exprs_names))) {
            stop("IGSAinput names must be unique");
        }
        
        # todo: check that geneSets is a list of GeneSetCollection, however 
        # this might be already evaluated by IGSAinput's setter
        
        resDir <- tempdir();
        
        # for each IGSAinput
        actRes <- lapply(seq_along(igsaInputs), function(i) {
            igsaInput <- igsaInputs[[i]];
            resFile <- paste0(resDir, '/', igsaInput@name, '.RData');
            
            # if gene sets were provided then these must be used
            if (length(geneSets) > 0) {
                geneSetsList(igsaInput) <- geneSets;
            }
            
            # if it was already ran ok, then just return
            if (file.exists(resFile)) return(TRUE);
            
            igsaRes <- try({ IGSA(igsaInput); });
            
            ranOk <- !inherits(igsaRes, 'try-error');
            if (ranOk) {
                # save intermediate results
                save(igsaRes, file=resFile);
            }
            return(ranOk);
        });
        
        if (!all(unlist(actRes))) {
            flog.warn(paste(sum(!unlist(actRes)), 'experiments had errors'));
        }
        
        actRes <- lapply(seq_along(igsaInputs), function(i) {
            igsaInput <- igsaInputs[[i]];
            resFile <- paste0(resDir, '/', igsaInput@name, '.RData');
            igsaRes <- NA;
            if (file.exists(resFile)) {
                igsaRes <- get(load(resFile));
            } else {
                flog.warn(paste(igsaInput@name, 'had errors'));
            }
            
            return(igsaRes);
        })
        
        # delete results which gave errors
        actRes <- actRes[!is.na(unlist(actRes))];
        
        migsaRes <- NA;
        if (length(actRes) > 0) {
            # if we have any result then create the MIGSAres
            migsaRes <- MIGSAres(actRes);
        }
        
        return(migsaRes);
    }
)
