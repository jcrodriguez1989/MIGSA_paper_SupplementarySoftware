% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GSEAparams-class.R
\docType{methods}
\name{GSEAparams-class}
\alias{GSEAparams-class}
\alias{GSEAparams}
\title{GSEAparams S4 class implementation in R}
\description{
This S4 class contains the parameters to provide for GSEA.
}
\section{Slots}{

\describe{
\item{\code{perm_number}}{number of permutations for p-value calculation 
(default: 200).}

\item{\code{min_sz}}{minimum size of gene sets (number of genes in a gene set) 
to be included in the analysis (default: 5).}

\item{\code{pv}}{estimate of the variance associated with each observation 
(default: 0).}

\item{\code{w1}}{weight 1, parameter used to calculate the prior variance 
obtained with class size var.constant. This penalizes especially small 
classes and small subsets. Values around 0.1 - 0.5 are expected to be 
reasonable. (default: 0.2).}

\item{\code{w2}}{weight 2, parameter used to calculate the prior variance 
obtained with the same class size as that of the analyzed class. This 
penalizes small subsets from the gene list. Values around 0.3 and 0.5 are 
expected to be reasonable (default: 0.5).}

\item{\code{vc}}{size of the reference class used with wgt1. (default: 10).}
}}

\examples{
## Lets create the default GSEAparams object.
myGseaParams <- GSEAparams();

## Lets create another GSEAparams object with 500 permutations.
myGseaParams500Perms <- GSEAparams(perm_number=500);

}
\seealso{
\code{\link{SEAparams-class}}

\code{\link{summary}}
}
