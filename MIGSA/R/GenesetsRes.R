#'@include SEAres.R
#'@include GSEAres.R
GenesetsRes <- setClass(
    Class="GenesetsRes",
    slots=c(
        gene_sets_name="character",
        is_GO="logical",
        sea_res="SEAres",
        gsea_res="GSEAres"
    ),
    prototype=list(
        is_GO=FALSE
    )
)

#'@importFrom data.table as.data.table data.table
as.data.table.GenesetsRes <- function(x, ...) {
    # convert to data.table SEA and GSEA results
    seaRes  <- as.data.table(x@sea_res);
    gseaRes <- as.data.table(x@gsea_res);
    
    # merge them and fill the data.table with Genesets info
    to <- merge(seaRes, gseaRes, by=c("id", "name"), all=!FALSE);
    to <- data.table(gene_set_name=x@gene_sets_name, is_GO=x@is_GO, to);
    
    return(to);
}
