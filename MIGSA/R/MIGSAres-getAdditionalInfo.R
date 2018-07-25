#'Gets additional information about enrichment results
#'
#'\code{getAdditionalInfo} gets additional enrichment information of the 
#'analysis done.
#'
#'@param migsaRes MIGSAres object.
#'
#'@return data.frame with additional information for each analyzed gene set.
#'
#'@docType methods
#'@name getAdditionalInfo
#'@rdname MIGSAres-getAdditionalInfo
#'
#'@exportMethod getAdditionalInfo
setGeneric(name="getAdditionalInfo", def=function(migsaRes) {
    standardGeneric("getAdditionalInfo")
})

#'@inheritParams getAdditionalInfo
#'@rdname MIGSAres-getAdditionalInfo
#'@aliases getAdditionalInfo,MIGSAres
#'@include MIGSAres-class.R
#'@include MIGSAres-setEnrCutoff.R
#'@seealso \code{\link{setEnrCutoff}}
#'@examples
#'data(migsaRes);
#'
#'## Lets get additional enrichment information of the MIGSAres object.
#'adtnlInfo <- getAdditionalInfo(migsaRes);
#'dim(adtnlInfo); # it is a huge data.frame
#'
#'## This huge amount of information commonly is not that interesting. Lets 
#'## keep the gene sets that were enriched in every experiment and check that 
#'## information. Lets set a cutoff of 0.1.
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'migsaResFiltered <- migsaResWCoff[rowSums(migsaResWCoff[,4:5]) == 2,];
#'adtnlInfo <- getAdditionalInfo(migsaResFiltered);
#'dim(adtnlInfo);
#'
setMethod(
    f="getAdditionalInfo",
    signature=c("MIGSAres"),
    definition=function(migsaRes) {
        res <- as.data.frame(migsaRes@migsa_res_all);
        return(res);
    }
)
