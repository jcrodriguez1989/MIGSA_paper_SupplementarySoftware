#'@include GenesetsRes.R
IGSAres <- setClass(
    Class="IGSAres",
    slots=c(
        name="character",
        gene_sets_res="list",
        genes_rank="data.frame"
    ),
    prototype=list(
        genes_rank=data.frame()
    ),
    validity=function(object) {
        # gene_sets_res is a list of GenesetsRes
        gene_sets_res_ok <- all(unlist(lapply(object@gene_sets_res,
            function(x) is(x, "GenesetsRes"))));
        return(gene_sets_res_ok);
    }
)

# to avoid R CMD check errors we set them as NULL
experiment_name = gene_set_name = NULL;

#'@importFrom data.table as.data.table data.table setkey
#'@include GenesetsRes.R
as.data.table.IGSAres <- function(x, ...) {
    # convert to data.table each GenesetsRes
    to <- do.call(rbind, lapply(x@gene_sets_res, function(gsetsRes) {
        actRes <- as.data.table(gsetsRes);
        return(actRes);
    }));
    
    # add some IGSA information
    to <- data.table(experiment_name=x@name, to);
    # set this keys for faster merging
    setkey(to, experiment_name, gene_set_name, id, name);
    
    return(to);
}

setGeneric(name="genes_rank", def=function(object) {
    standardGeneric("genes_rank")
})

setMethod(f="genes_rank", signature="IGSAres",
    definition=function(object) {
        return(object@genes_rank)
    }
)
