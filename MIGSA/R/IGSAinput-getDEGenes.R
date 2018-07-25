#'Calculate differentialy expressed genes of an IGSAinput object
#'
#'\code{getDEGenes} calculates the differentialy expressed genes of an 
#'IGSAinput object using the expression matrix, the FitOtions and the 
#'SEAparams..
#'
#'@param igsaInput a valid IGSAinput object.
#'
#'@return A IGSAinput object with the updated slots and DE genes calculated.
#'
#'@docType methods
#'@name getDEGenes
#'@rdname IGSAinput-getDEGenes
#'@seealso \code{\link{IGSAinput-class}}
#'
#'@include IGSAinput-class.R
#'@include IGSAinput.R
#'@include SEAparams.R
#'@exportMethod getDEGenes
#'
setGeneric(name="getDEGenes", def=function(igsaInput) {
    standardGeneric("getDEGenes")
})

#'@inheritParams getDEGenes
#'@rdname IGSAinput-getDEGenes
#'@aliases getDEGenes,IGSAinput
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
#'## Lets create our IGSAinput object.
#'myIgsaInput <- IGSAinput(name="myIgsaInput", expr_data=maExprData, 
#'fit_options=myFOpts);
#'
#'## And check how many differentialy expressed genes does it get with actual
#'## parameters.
#'aux <- getDEGenes(myIgsaInput);
#'
setMethod(f="getDEGenes",
    signature=c("IGSAinput"),
    definition=function(igsaInput) {
        validObject(igsaInput);
        sea_par <- igsaInput@sea_params;
        
        if (is.null(sea_par)) return(sea_par);
        
        # get differentialy expressed genes
        dif <- igsaGetDEGenes(sea_par, igsaInput@expr_data, 
                    igsaInput@fit_options);
        
        # update sea_par with this dif genes, so we dont need to recalculate
        de_genes(sea_par) <- dif;
        igsaInput@sea_params <- sea_par;
        
        return(igsaInput);
    }
)
