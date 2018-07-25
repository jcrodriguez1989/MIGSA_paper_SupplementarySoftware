#'@include GenesetRes.R
GSEAres <- setClass(
    Class="GSEAres",
    slots=c(
        gene_sets_res="list"
    ),
    prototype=list(
    ),
    validity=function(object) {
        gene_sets_res_ok <- all(unlist(lapply(object@gene_sets_res, function(x)
            is(x, "GenesetRes"))));

        return(gene_sets_res_ok);
    }
)

#'@importFrom data.table as.data.table data.table setkey
as.data.table.GSEAres <- function(x, ...) {
    # convert each gene set result to character vector
    to <- do.call(rbind, lapply(x@gene_sets_res, function(res) {
        actInfo <- asCharacter(res);
        return(actInfo);
    }))
    
    if (is.null(to)) return(
        data.table(id=character(),
        name=character(),
        GSEA_enriched=logical(),
        GSEA_score=logical(),
        GSEA_pval=logical(),
        GSEA_enriching_genes=character(),
        GSEA_GS_genes=character()
        ));
    
    # make it a data.table and put the colnames
    to <- data.table(to);
    colnames(to) <- c("id", "name", "GSEA_enriched", "GSEA_score",
                        "GSEA_pval", "GSEA_enriching_genes",
                        "GSEA_GS_genes");
    
    # set data.table keys, to merge faster
    setkey(to, id, name);
    
    return(to)
}

setGeneric(name="gene_sets_res", def=function(object) {
    standardGeneric("gene_sets_res")
})

setMethod(f="gene_sets_res", signature="GSEAres",
    definition=function(object) {
        return(object@gene_sets_res)
    }
)
