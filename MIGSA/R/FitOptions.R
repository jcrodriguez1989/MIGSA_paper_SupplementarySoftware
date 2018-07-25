#'@include FitOptions-class.R
setGeneric(name="designMatrix", def=function(object) {
    standardGeneric("designMatrix")
})

setMethod(f="designMatrix", signature="FitOptions",
    definition=function(object) {
        return(object@design_matrix)
    }
)

setGeneric(name="designMatrix<-", def=function(object, value) {
    standardGeneric("designMatrix<-")
})

setReplaceMethod(f="designMatrix", signature="FitOptions",
    definition=function(object, value) {
        object@design_matrix <- value
        validObject(object)
        return(object)
    }
)

setGeneric(name="contrast", def=function(object) {
    standardGeneric("contrast")
})

setMethod(f="contrast", signature="FitOptions",
    definition=function(object) {
        return(object@contrast)
    }
)

setGeneric(name="col_data", def=function(object) {
    standardGeneric("col_data")
})

setMethod(f="col_data", signature="FitOptions",
    definition=function(object) {
        return(object@col_data)
    }
)
