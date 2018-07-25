#'Set enrichment cutoff for a MIGSAres object
#'
#'\code{setEnrCutoff} sets the enrichment cutoff value for a MIGSAres object in
#'order to get detailed enrichment results.
#'
#'@param object a MIGSAres object.
#'@param newEnrCutoff numeric value in range [0,1] or NA to set the new
#'enrichment cutoff. If NA then enrichment cutoff is unset for object.
#'
#'@return A MIGSAres object with the updated slots.
#'
#'@docType methods
#'@name setEnrCutoff
#'@rdname MIGSAres-setEnrCutoff
#'
#'@include MIGSAres-class.R
#'@exportMethod setEnrCutoff
#'
setGeneric(name="setEnrCutoff", def=function(object, newEnrCutoff) {
    standardGeneric("setEnrCutoff");
})

#'@inheritParams setEnrCutoff
#'@rdname MIGSAres-setEnrCutoff
#'@aliases setEnrCutoff,MIGSAres,numeric-method
#'@examples
#'data(migsaRes);
#'## After executing MIGSA, the MIGSAres object does not have an enrichment 
#'## cutoff set, so it will be shown as follows:
#'head(migsaRes);
#'
#'## Now lets set an enrichment cutoff, and then show again the object.
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.01);
#'head(migsaResWCoff);
#'
#'## So now we can do stuff like check how many gene sets are enriched in both 
#'## experiments, in each one, etc. Moreover now we can do different MIGSAres 
#'## plots available in MIGSA package.
#'table(migsaResWCoff[, c("igsaInput1", "igsaInput2")]);
#'
#'## And with different cutoffs.
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'table(migsaResWCoff[, c("igsaInput1", "igsaInput2")]);
#'
setMethod("setEnrCutoff",
    signature=c("MIGSAres", "numeric"),
    definition=function(object, newEnrCutoff) {
        stopifnot(validObject(object));
        
        res <- object;
        enrCutoff(res) <- newEnrCutoff;
        
        return(res);
    }
)
