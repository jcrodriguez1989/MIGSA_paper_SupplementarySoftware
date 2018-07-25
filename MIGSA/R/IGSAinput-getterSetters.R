#'Accessors for IGSAinput class
#'
#'Getters and setters functions to access IGSAinput object slots.
#'
#'@param object IGSAinput object.
#'@param value value to replace in the slot.
#'
#'@return Modified IGSAinput object or desired slot.
#'
#'@docType methods
#'@name IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@seealso \code{\link{IGSAinput-class}}
#'@include IGSAinput-class.R
#'
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
#'## And now we can create our IGSAinput ready for MIGSA.
#'igsaInput <- IGSAinput(name="igsaInput", expr_data=maExprData, 
#'fit_options=myFOpts);
#'
#'## Lets get igsaInput values, and modify its name.
#'name(igsaInput);
#'name(igsaInput) <- "newName";
#'fitOptions(igsaInput); 
#'exprData(igsaInput);
#'geneSetsList(igsaInput);
#'gseaParams(igsaInput);
#'seaParams(igsaInput);
#'
setGeneric(name="IGSAinput-getterSetters", def=function(object) {
    standardGeneric("IGSAinput-getterSetters")
})

#'@name name
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases name,IGSAinput-method
#'@exportMethod name
#'
setGeneric(name="name", def=function(object) {
    standardGeneric("name")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases name,IGSAinput-method
#'
setMethod(f="name", signature="IGSAinput",
    definition=function(object) {
        return(object@name)
    }
)

#'@name name<-
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases name<-,IGSAinput-method
#'@exportMethod name<-
#'
setGeneric(name="name<-", def=function(object, value) {
    standardGeneric("name<-")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases name<-,IGSAinput-method
#'
setReplaceMethod(f="name", signature="IGSAinput",
    definition=function(object, value) {
        object@name <- value;
        validObject(object);
        return(object);
    }
)

#'@name fitOptions
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases fitOptions,IGSAinput-method
#'@exportMethod fitOptions
#'
setGeneric(name="fitOptions", def=function(object) {
    standardGeneric("fitOptions")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases fitOptions,IGSAinput-method
#'
setMethod(f="fitOptions", signature="IGSAinput",
    definition=function(object) {
        return(object@fit_options)
    }
)

#'@name fitOptions<-
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases fitOptions<-,IGSAinput-method
#'@exportMethod fitOptions<-
#'
setGeneric(name="fitOptions<-", def=function(object, value) {
    standardGeneric("fitOptions<-")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases fitOptions<-,IGSAinput-method
#'
setReplaceMethod(f="fitOptions", signature="IGSAinput",
    definition=function(object, value) {
        object@fit_options <- value;
        validObject(object);
        return(object);
    }
)

#'@name exprData
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases exprData,IGSAinput-method
#'@exportMethod exprData
#'
setGeneric(name="exprData", def=function(object) {
    standardGeneric("exprData")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases exprData,IGSAinput-method
#'
setMethod(f="exprData", signature="IGSAinput",
    definition=function(object) {
        return(object@expr_data)
    }
)

#'@name exprData<-
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases exprData<-,IGSAinput-method
#'@exportMethod exprData<-
#'
setGeneric(name="exprData<-", def=function(object, value) {
    standardGeneric("exprData<-")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases exprData<-,IGSAinput-method
#'
setReplaceMethod(f="exprData", signature="IGSAinput",
    definition=function(object, value) {
        object@expr_data <- value;
        validObject(object);
        return(object);
    }
)

#'@name geneSetsList
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases geneSetsList,IGSAinput-method
#'@exportMethod geneSetsList
#'
setGeneric(name="geneSetsList", def=function(object) {
    standardGeneric("geneSetsList")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases geneSetsList,IGSAinput-method
#'
setMethod(f="geneSetsList", signature="IGSAinput",
    definition=function(object) {
        return(object@gene_sets_list)
    }
)

#'@name geneSetsList<-
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases geneSetsList<-,IGSAinput-method
#'@exportMethod geneSetsList<-
#'
setGeneric(name="geneSetsList<-", def=function(object, value) {
    standardGeneric("geneSetsList<-")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases geneSetsList<-,IGSAinput-method
#'
setReplaceMethod(f="geneSetsList", signature="IGSAinput",
    definition=function(object, value) {
        object@gene_sets_list <- value;
        validObject(object);
        return(object);
    }
)

#'@name gseaParams
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases gseaParams,IGSAinput-method
#'@exportMethod gseaParams
#'
setGeneric(name="gseaParams", def=function(object) {
    standardGeneric("gseaParams")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases gseaParams,IGSAinput-method
#'
setMethod(f="gseaParams", signature="IGSAinput",
    definition=function(object) {
        return(object@gsea_params)
    }
)

#'@name gseaParams<-
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases gseaParams<-,IGSAinput-method
#'@exportMethod gseaParams<-
#'
setGeneric(name="gseaParams<-", def=function(object, value) {
    standardGeneric("gseaParams<-")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases gseaParams<-,IGSAinput-method
#'
setReplaceMethod(f="gseaParams", signature="IGSAinput",
    definition=function(object, value) {
        object@gsea_params <- value;
        validObject(object);
        return(object);
    }
)

#'@name seaParams
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases seaParams,IGSAinput-method
#'@exportMethod seaParams
#'
setGeneric(name="seaParams", def=function(object) {
    standardGeneric("seaParams")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases seaParams,IGSAinput-method
#'
setMethod(f="seaParams", signature="IGSAinput",
    definition=function(object) {
        return(object@sea_params)
    }
)

#'@name seaParams<-
#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases seaParams<-,IGSAinput-method
#'@exportMethod seaParams<-
#'
setGeneric(name="seaParams<-", def=function(object, value) {
    standardGeneric("seaParams<-")
})

#'@inheritParams IGSAinput-getterSetters
#'@rdname IGSAinput-getterSetters
#'@aliases seaParams<-,IGSAinput-method
#'
setReplaceMethod(f="seaParams", signature="IGSAinput",
    definition=function(object, value) {
        object@sea_params <- value;
        validObject(object);
        return(object);
    }
)
