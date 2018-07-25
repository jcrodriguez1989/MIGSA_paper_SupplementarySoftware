#'MIGSAres S4 class implementation in R
#'
#'This S4 class represents the results of one MIGSA execution in R.
#'
#'@importFrom data.table data.table
#'@name MIGSAres-class
#'@rdname MIGSAres-class
#'@seealso \code{\link{MIGSA}}
#'@seealso \code{\link{getAdditionalInfo}}
#'@seealso \code{\link{MIGSAres-common}}
#'@seealso \code{\link{MIGSAres-genes}}
#'@seealso \code{\link{MIGSAres-GOanalysis}}
#'@seealso \code{\link{MIGSAres-plots}}
#'@seealso \code{\link{setEnrCutoff}}
#'@seealso \code{\link{summary}}
#'@exportClass MIGSAres
#'
setClass(
    Class="MIGSAres",
    slots=c(
        migsa_res_all="data.table",
        migsa_res_summary="data.frame",
        enr_cutoff="numeric",
        genes_rank="list"
    ),
    prototype=list(
        enr_cutoff=as.numeric(NA)
    ),
    validity=function(object) {
        # check that data frames are correct
        migsa_res_all_ok <- all(colnames(object@migsa_res_all) %in%
        c("experiment_name", "gene_set_name", "id", "name", "SEA_GS_genes",
            "SEA_enriched", "SEA_score", "SEA_pval", "SEA_enriching_genes",
            "GSEA_GS_genes", "GSEA_enriched", "GSEA_score", "GSEA_pval",
            "GSEA_enriching_genes", "is_GO"));
        migsa_res_all_ok <- migsa_res_all_ok &&
                        (ncol(object@migsa_res_all) == 15);
        
        migsa_res_summary_ok <- all(c("id") %in%
                                        colnames(object@migsa_res_summary));
        
        enr_cutoff_ok <- is.na(object@enr_cutoff) || 
                            (0 <= object@enr_cutoff && object@enr_cutoff <= 1);
        genes_rank_ok <- all(unlist(lapply(object@genes_rank, is,
                                                            "data.frame")));
        
        return(migsa_res_all_ok && enr_cutoff_ok &&
                    migsa_res_summary_ok && genes_rank_ok);
    }
)
