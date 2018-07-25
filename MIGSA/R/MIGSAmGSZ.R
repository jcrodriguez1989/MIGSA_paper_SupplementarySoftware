#'MIGSAmGSZ
#'
#'\code{MIGSAmGSZ} is an optimized mGSZ version. It runs much faster than the 
#'original mGSZ version, moreover it can run in multicore technology.
#'It allows to analyze RNAseq data by using \code{\link[limma]{voom}} function.
#'mGSZ: Gene set analysis based on Gene Set Z scoring function and asymptotic 
#'p-value.
#'
#'@param x gene expression data matrix (rows as genes and columns as samples).
#'@param y gene set data (list).
#'@param l vector of response values (example:c("Cond1","Cond1","Cond2", 
#'"Cond2","Cond2")).
#'@param use.voom logical indicating wether use voom or not (if RNAseq data we 
#'recommend using use.voom=TRUE).
#'@param rankFunction internal use.
#'@param min.sz minimum size of gene sets (number of genes in a gene set) to 
#'be included in the analysis.
#'@param pv estimate of the variance associated with each observation.
#'@param w1 weight 1, parameter used to calculate the prior variance obtained 
#'with class size var.constant. This penalizes especially small classes and 
#'small subsets. Values around 0.1 - 0.5 are expected to be reasonable.
#'@param w2 weight 2, parameter used to calculate the prior variance obtained 
#'with the same class size as that of the analyzed class. This penalizes small 
#'subsets from the gene list. Values around 0.3 and 0.5 are expected to be 
#'reasonable.
#'@param vc size of the reference class used with wgt1.
#'@param p number of permutations for p-value calculation.
#'@param ... not in use.
#'
#'@return A data.frame with gene sets p-values and additional information.
#'
#'@docType methods
#'@name MIGSAmGSZ
#'@rdname MIGSAmGSZ
#'
#'@exportMethod MIGSAmGSZ
#'
setGeneric(name="MIGSAmGSZ", def=function(x, y, l, ...) {
    standardGeneric("MIGSAmGSZ")
})

#'@inheritParams MIGSAmGSZ
#'@rdname MIGSAmGSZ
#'@aliases MIGSAmGSZ,matrix,list,vector-method
#'
#'@importFrom BiocParallel bpparam
#'@importClassesFrom edgeR DGEList
#'@importFrom edgeR DGEList
#'@importClassesFrom limma MAList
#'@include FitOptions-class.R
#'@include GSEAparams.R
#'@include MGSZ.R
#'
#'@examples
#'nGenes <- 1000; # 1000 genes
#'nSamples <- 30; # 30 subjects
#'geneNames <- paste("g", 1:nGenes, sep=""); # with names g1 ... g1000
#'## Create random gene expression data matrix.
#'set.seed(8818);
#'exprData <- matrix(rnorm(nGenes*nSamples),ncol=nSamples);
#'rownames(exprData) <- geneNames;
#'
#'## There will be 40 differentialy expressed genes.
#'nDeGenes <- nGenes/25;
#'## Lets generate the offsets to sum to the differentialy expressed genes.
#'deOffsets <- matrix(2*abs(rnorm(nDeGenes*nSamples/2)), ncol=nSamples/2);
#'
#'## Randomly select which are the DE genes.
#'deIndexes <- sample(1:nGenes, nDeGenes, replace=FALSE);
#'exprData[deIndexes, 1:(nSamples/2)] <-
#'exprData[deIndexes, 1:(nSamples/2)] + deOffsets;
#'
#'## 15 subjects with condition C1 and 15 with C2.
#'conditions <- rep(c("C1", "C2"),c(nSamples/2,nSamples/2));
#'
#'nGSets <- 200; # 200 gene sets
#'## Lets create randomly 200 gene sets, of 10 genes each
#'gSets <- lapply(1:nGSets, function(i) sample(geneNames, size=10));
#'names(gSets) <- paste("set", as.character(1:nGSets), sep="");
#'
#'\dontrun{
#'mGSZres <- MIGSAmGSZ(exprData, gSets, conditions);
#'}
#'
setMethod(
    f="MIGSAmGSZ",
    signature=c("matrix", "list", "vector"),
    definition=function (x, y, l, use.voom=FALSE, rankFunction=NA,
    min.sz=5, pv=0, w1=0.2, w2=0.5, vc=10, p=200) {
        # it formats almost the same inputs as mGSZ and uses MIGSAs mGSZ.
        # setting all MIGSA parameters
        
        if (use.voom) {
            exprData <- DGEList(counts=x);
            if (!is(rankFunction, "function")) {
                rankFunction <- voomLimaRank;
            }
        } else {
            exprData <- new("MAList", list(M=x));
            if (!is(rankFunction, "function")) {
                rankFunction <- mGszEbayes;
            }
        }
        fitOptions <- FitOptions.default(l);
        
        params <- GSEAparams(
            perm_number=p,
            min_sz=min.sz,
            pv=pv,
            w1=w1,
            w2=w2,
            vc=vc
        );
        gSets <- y;
        
        mgszRes <- MIGSA_mGSZ(exprData, fitOptions, gSets,
            rankFunction, params);
        return(mgszRes);
    }
)
