GenesetRes <- setClass(
    Class="GenesetRes",
    slots=c(
        id="character",
        name="character",
        enriched="logical",
        score="numeric",
        pvalue="numeric",
        genes="character",
        enriching_genes="character"
    ),
    prototype=list(
        score=as.numeric(NA),
        enriched=FALSE
    ),
    validity=function(object) {
        id_ok <- object@id != "";
        pvalue_ok <- (is.na(object@pvalue) ||
                            (object@pvalue >=0 && object@pvalue <= 1));
        genes_ok <- length(object@genes) > 0;
        
        return(id_ok && pvalue_ok && genes_ok);
    }
)

setGeneric(name="asCharacter", def=function(x, ...) {
    standardGeneric("asCharacter")
})

setMethod(
    f="asCharacter",
    signature=c("GenesetRes"),
    definition=function(x, ...) {
        # translate a gene set result into a character vector
        actName <- ifelse(length(x@name) == 0, "", x@name);
        to <- c(x@id, actName, x@enriched, x@score, x@pvalue);
        to <- c(to, paste(x@enriching_genes, collapse=", "),
                    paste(x@genes, collapse=", "));
        
        # to avoid rare characters
        to <- iconv(to, to="ASCII", sub="");
        
        return(to)
    }
)

setGeneric(name="id", def=function(object) {
    standardGeneric("id")
})

setMethod(f="id", signature="GenesetRes", definition=function(object) {
    return(object@id)
})

setGeneric(name="id<-", def=function(object, value) {
    standardGeneric("id<-")
})

setReplaceMethod(f="id", signature="GenesetRes",
    definition=function(object, value) {
        object@id <- value
        validObject(object)
        return(object)
    }
)

setGeneric(name="getName", def=function(object) {
    standardGeneric("getName")
})

setMethod(f="getName", signature="GenesetRes", definition=function(object) {
    return(object@name)
})

setGeneric(name="getName<-", def=function(object, value) {
    standardGeneric("getName<-")
})

setReplaceMethod(f="getName", signature="GenesetRes",
    definition=function(object, value) {
        object@name <- value
        validObject(object)
        return(object)
    }
)
