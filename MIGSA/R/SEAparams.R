#'@include SEAparams-class.R
setGeneric(name="treat_lfc", def=function(object) {
    standardGeneric("treat_lfc")
})

setMethod(f="treat_lfc", signature="SEAparams",
    definition=function(object) {
        return(object@treat_lfc)
    }
)

setGeneric(name="br", def=function(object) {
    standardGeneric("br")
})

setMethod(f="br", signature="SEAparams",
    definition=function(object) {
        return(object@br)
    }
)

setGeneric(name="adjust_method", def=function(object) {
    standardGeneric("adjust_method")
})

setMethod(f="adjust_method", signature="SEAparams",
    definition=function(object) {
        return(object@adjust_method)
    }
)

setGeneric(name="de_genes", def=function(object) {
    standardGeneric("de_genes")
})

setMethod(f="de_genes", signature="SEAparams",
    definition=function(object) {
        return(object@de_genes)
    }
)

setGeneric(name="de_genes<-", def=function(object, value) {
    standardGeneric("de_genes<-")
})

setReplaceMethod(f="de_genes", signature="SEAparams",
    definition=function(object, value) {
        object@de_genes <- value
        validObject(object)
        return(object)
    }
)

setGeneric(name="de_cutoff", def=function(object) {
    standardGeneric("de_cutoff")
})

setMethod(f="de_cutoff", signature="SEAparams",
    definition=function(object) {
        return(object@de_cutoff)
    }
)

setGeneric(name="test", def=function(object) {
    standardGeneric("test")
})

setMethod(f="test", signature="SEAparams",
    definition=function(object) {
        return(object@test)
    }
)
