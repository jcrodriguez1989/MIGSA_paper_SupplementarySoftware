#'SEAparams S4 class implementation in R
#' 
#'This S4 class contains the parameters to provide for SEA.
#'
#'@slot treat_lfc numeric, lfc parameter passed to 
#'\code{treat} function (default: 0).
#'@slot de_cutoff numeric, cutoff value to define a gene as differentialy
#'expressed (default: 0.01).
#'@slot adjust_method character, method parameter passed to 
#'\code{\link[stats]{p.adjust}} function (default: "fdr").
#'@slot de_genes (optional) character vector of differentialy expressed 
#'genes,
#'if not provided it will be filled using the other parameters.
#'@slot br background reference to use, there are two possible options 
#'for this slot (default: "briii"):
#'\itemize{
#'\item character indicating to use "bri" or "briii".
#'\item character vector indicating the genes to use as background 
#'reference.
#'}
#'@slot test character indicating which test to use, it must be one of 
#'"FisherTest", "HypergeoTest" or "BinomialTest" (default: "FisherTest").
#'
#'@docType methods
#'@name SEAparams-class
#'@rdname SEAparams-class
#'@seealso \code{GSEAparams-class}
#'@seealso \code{summary}
#'
#'@export SEAparams
#'@examples
#'## Lets create the default SEAparams object.
#'mySeaParams <- SEAparams();
#'
#'## Lets create another SEAparams object with my default br genes. Make sure 
#'## these genes are present in the analyzed expression matrix rownames.
#'myBr <- as.character(1:10);
#'mySeaParamsOwnBr <- SEAparams(br=myBr);
#'
#'## Lets create another SEAparams object with my default differentialy 
#'## expressed genes. Make sure these genes are present in the analyzed 
#'## expression matrix rownames.
#'myDEGenes <- as.character(3:15);
#'mySeaParamsOwnDEG <- SEAparams(de_genes=myDEGenes);
#'
#'## Lets create another SEAparams object changing the differentialy 
#'## expressed genes parameters.
#'mySeaParamsOwnParams <- SEAparams(treat_lfc=0.25, de_cutoff=0.05, 
#'adjust_method="none");
#'
SEAparams <- setClass(
    Class="SEAparams",
    slots=c(
        treat_lfc="numeric",
        de_cutoff="numeric",
        adjust_method="character",
        de_genes="character",
        br="character",
        test="character"
    ),
    prototype=list(
        treat_lfc=0,
        de_cutoff=0.01,
        adjust_method="fdr",
        de_genes=as.character(NA),
        br="briii",
        test="FisherTest"
    ),
    validity=function(object) {
        treat_lfc_ok  <- object@treat_lfc >=0;
        de_cutoff_ok  <- object@de_cutoff >= 0 && object@de_cutoff <= 1;
        adjust_method_ok <- object@adjust_method %in% p.adjust.methods;
        br_ok <- object@br[[1]] %in% c("bri", "briii");
        
        # br must be bri, briii or a list of genes
        if (!br_ok) {
            br_ok <- length(unique(object@br)) > 0 && 
                        any(!unique(object@br) == "");
        }
        
        # checking if test is valid
        test_ok <- length(object@test) == 1 && 
            object@test %in% c("FisherTest", "HypergeoTest", "BinomialTest");
        
        return(treat_lfc_ok && de_cutoff_ok && adjust_method_ok && br_ok && 
            test_ok);
    }
)

setClassUnion("SEAparamsOrNULL", c("SEAparams", "NULL"));
