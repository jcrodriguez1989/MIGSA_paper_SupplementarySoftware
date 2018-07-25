#'IGSAinput S4 class implementation in R
#' 
#'This S4 class contains all the necessary inputs to execute a functional 
#'analysis (SEA and GSEA) on one experiment.
#'Important: Make sure that gene IDs are concordant between the expression 
#'matrix and the provided gene sets.
#'
#'@slot name character indicating the name of this experiment.
#'@slot expr_data ExprData object with the expression data (MicroArray or 
#'RNAseq). Note: expr_data can be a 0x0 matrix only if gsea_params=NULL and 
#'de_genes and br slots from sea_params are correctly set (vectors of gene 
#'names), in this case only SEA will be run.
#'@slot fit_options FitOptions object with the parameters to be used when 
#'fitting the model.
#'@slot gene_sets_list named list of GeneSetCollection objects to be 
#'tested for enrichment (names must be unique).
#'@slot sea_params SEAparams object with the parameters to be used 
#'by SEA, if NULL then SEA wont be run.
#'@slot gsea_params GSEAparams object with the parameters to be used 
#'by GSEA, if NULL then GSEA wont be run.
#'
#'@docType methods
#'@name IGSAinput-class
#'@rdname IGSAinput-class
#'@seealso \code{\link{ExprData-class}}
#'@seealso \code{\link{SEAparams-class}}
#'@seealso \code{\link{GSEAparams-class}}
#'@seealso \code{\link{IGSAinput-getterSetters}}
#'@seealso \code{\link{getDEGenes}}
#'@seealso \code{\link{MIGSA}}
#'@seealso \code{\link{summary}}
#'
#'@importFrom futile.logger flog.error
#'@importFrom GSEABase GeneSet GeneSetCollection
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include SEAparams.R
#'@include GSEAparams.R
#'@include Genesets.R
#'@export IGSAinput
#'@examples
#'## Lets create a basic IGSAinput object.
#'## First create a expression matrix.
#'maData <- matrix(rnorm(10000),ncol=4);
#'rownames(maData) <- 1:nrow(maData); # It must have rownames (gene names).
#'maExprData <- new("MAList",list(M=maData));
#'
#'## Now lets create the FitOptions object.
#'myFOpts <- FitOptions(c("Cond1", "Cond1", "Cond2", "Cond2"));
#'
#'## Finally lets create the Genesets to test for enrichment.
#'library(GSEABase);
#'myGs1 <- GeneSet(as.character(1:10), setIdentifier="fakeId1", 
#'     setName="fakeName1");
#'myGs2 <- GeneSet(as.character(7:15), setIdentifier="fakeId2", 
#'     setName="fakeName2");
#'myGSs <- GeneSetCollection(list(myGs1, myGs2));
#'
#'## And now we can create our IGSAinput ready for MIGSA.
#'igsaInput <- IGSAinput(name="igsaInput", expr_data=maExprData, 
#'fit_options=myFOpts, gene_sets_list=list(myGSs=myGSs));
#'
#'## Valid IGSAinput object with no expr_data (only run SEA).
#'igsaInput <- IGSAinput(name="igsaInput", gene_sets_list=list(myGSs=myGSs),
#'gsea_params=NULL, 
#'sea_params=SEAparams(de_genes=rownames(maExprData)[1:10],
#'br=rownames(maExprData)));
#'validObject(igsaInput);
#'
IGSAinput <- setClass(
    Class="IGSAinput",
    slots=c(
        name="character",
        expr_data="ExprData",
        fit_options="FitOptions",
        gene_sets_list="list",
        sea_params="SEAparamsOrNULL",
        gsea_params="GSEAparamsOrNULL"
    ),
    prototype=list(
        sea_params=SEAparams(),
        gsea_params=GSEAparams()
    ),
    validity=function(object) {
        # must have name
        name_ok <- length(object@name) == 1 && object@name != "";
        
        # must have genes and samples
        expr_data_ok <- (!is.null(object@expr_data)) &&
                            ncol(object@expr_data) > 1 &&
                            nrow(object@expr_data) > 1;
        
        # check that the FitOptions and the ExprData are concordant
        fit_opts_ok <- nrow(MIGSA:::designMatrix(object@fit_options)) ==
                                                    ncol(object@expr_data);
        
        # gene_sets_list is a list of GeneSetCollection
        gene_sets_list_ok <- all(unlist(lapply(object@gene_sets_list,
                                    function(x) is(x, "GeneSetCollection"))));
        
        only_sea <- FALSE;
        if (!is.null(object@sea_params)) {
            only_sea <- length(MIGSA:::br(object@sea_params)) > 1 &&
            length(MIGSA:::de_genes(object@sea_params)) > 1;
        }
        
        # GeneSetCollection names must be unique
        if (gene_sets_list_ok) {
            gss_names <- names(object@gene_sets_list);
            gene_sets_list_ok <- length(gss_names) ==
                                        length(unique(gss_names));
            
            # every gene set collection must have a name
            gene_sets_list_ok <- gene_sets_list_ok && 
                length(gss_names) == length(object@gene_sets_list);
            
            if (!gene_sets_list_ok) {
                flog.error("GeneSetCollection names must be unique");
            }
        }
        
        return(name_ok && (only_sea || (expr_data_ok && gene_sets_list_ok)) &&
            fit_opts_ok);
    }
)
