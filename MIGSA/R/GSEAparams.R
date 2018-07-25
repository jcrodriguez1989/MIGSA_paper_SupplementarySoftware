#'@include GSEAparams-class.R
setGeneric(name="perm_number", def=function(object) {
    standardGeneric("perm_number")
})

setMethod(f="perm_number", signature="GSEAparams",
    definition=function(object) {
        return(object@perm_number)
    }
)

setGeneric(name="minSz", def=function(object) {
    standardGeneric("minSz")
})

setMethod(f="minSz", signature="GSEAparams",
    definition=function(object) {
        return(object@min_sz)
    }
)

setGeneric(name="pv", def=function(object) {
    standardGeneric("pv")
})

setMethod(f="pv", signature="GSEAparams",
    definition=function(object) {
        return(object@pv)
    }
)

setGeneric(name="w1", def=function(object) {
    standardGeneric("w1")
})

setMethod(f="w1", signature="GSEAparams",
    definition=function(object) {
        return(object@w1)
    }
)

setGeneric(name="w2", def=function(object) {
    standardGeneric("w2")
})

setMethod(f="w2", signature="GSEAparams",
    definition=function(object) {
        return(object@w2)
    }
)

setGeneric(name="vc", def=function(object) {
    standardGeneric("vc")
})

setMethod(f="vc", signature="GSEAparams",
    definition=function(object) {
        return(object@vc)
    }
)
