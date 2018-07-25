#'GSEAparams S4 class implementation in R
#' 
#'This S4 class contains the parameters to provide for GSEA.
#'
#'@slot perm_number number of permutations for p-value calculation 
#'(default: 200).
#'@slot min_sz minimum size of gene sets (number of genes in a gene set) 
#'to be included in the analysis (default: 5).
#'@slot pv estimate of the variance associated with each observation 
#'(default: 0).
#'@slot w1 weight 1, parameter used to calculate the prior variance 
#'obtained with class size var.constant. This penalizes especially small 
#'classes and small subsets. Values around 0.1 - 0.5 are expected to be 
#'reasonable. (default: 0.2).
#'@slot w2 weight 2, parameter used to calculate the prior variance 
#'obtained with the same class size as that of the analyzed class. This 
#'penalizes small subsets from the gene list. Values around 0.3 and 0.5 are 
#'expected to be reasonable (default: 0.5).
#'@slot vc size of the reference class used with wgt1. (default: 10).
#'
#'@docType methods
#'@name GSEAparams-class
#'@rdname GSEAparams-class
#'@seealso \code{\link{SEAparams-class}}
#'@seealso \code{\link{summary}}
#'@export GSEAparams
#'@examples
#'## Lets create the default GSEAparams object.
#'myGseaParams <- GSEAparams();
#'
#'## Lets create another GSEAparams object with 500 permutations.
#'myGseaParams500Perms <- GSEAparams(perm_number=500);
#'
GSEAparams <- setClass(
    Class="GSEAparams",
    slots=c(
        perm_number="numeric",
        min_sz="numeric",
        pv="numeric",
        w1="numeric",
        w2="numeric",
        vc="numeric"
    ),
    prototype=list(
        perm_number=200,
        min_sz=5,
        pv=0,
        w1=0.2,
        w2=0.5,
        vc=10
    ),
    validity=function(object) {
        perm_number_ok <- object@perm_number > 1;
        min_sz_ok <- object@min_sz >= 0;
        pv_ok <- object@pv >= 0;
        w1_ok <- object@w1 >= 0;
        w2_ok <- object@w2 >= 0;
        vc_ok <- object@vc >= 0;
        
        return(perm_number_ok && min_sz_ok && pv_ok && w1_ok && w2_ok &&
            vc_ok);
    }
)

setClassUnion("GSEAparamsOrNULL", c("GSEAparams", "NULL"));
