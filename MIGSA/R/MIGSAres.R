#'@include MIGSAres-class.R

# to avoid R CMD check errors we set them as NULL
experiment_name = gene_set_name = NULL;

############################ Creation functions

MIGSAres <- function(x, ... ) {
    UseMethod("MIGSAres", x);
}

#'@importFrom data.table data.table
MIGSAres.data.table <- function(x, genes_rank=list(), ... ) {
    migsa_res_summ <- summary_from_all(x);
    
    .Object <- new("MIGSAres", migsa_res_all=x,
                    migsa_res_summary=migsa_res_summ,
                    genes_rank=genes_rank);
    return(.Object);
}

#'@include IGSAres.R
# creates a MIGSAres from a list of IGSAres
MIGSAres.list <- function(x, ... ) {
    # check that is a list of IGSAres
    stopifnot(all(unlist(lapply(x, function(y) is(y, "IGSAres")))));
    
    # convert each IGSAres to a data.table and rbind them
    migsa_res_dtable <- do.call(rbind, lapply(x, function(igsaRes) {
        actRes <- as.data.table(igsaRes);
        return(actRes);
    }));
    
    # create a list of every genes_rank
    genes_rank <- lapply(x, function(igsaRes) {
        actRank <- genes_rank(igsaRes);
        return(actRank);
    })
    
    .Object <- MIGSAres(migsa_res_dtable, genes_rank=genes_rank);
    return(.Object);
}

setGeneric(name="summary_from_all", def=function(x, ... ) {
    standardGeneric("summary_from_all")
})

# to avoid R CMD check errors we set them as NULL
SEA_pval = GSEA_pval = Name = GS_Name = NULL;

#'@importFrom data.table data.table setkey
setMethod(
    f="summary_from_all",
    signature=c("data.table"),
    definition=function(x, by_method=FALSE) {
        stopifnot(all(colnames(x) %in%
            c("experiment_name", "gene_set_name", "id", "name", "SEA_GS_genes",
                "SEA_enriched", "SEA_score", "SEA_pval", "SEA_enriching_genes",
                "GSEA_GS_genes", "GSEA_enriched", "GSEA_score", "GSEA_pval",
                "GSEA_enriching_genes", "is_GO")));
        stopifnot(ncol(x) == 15);
        
        setkey(x, experiment_name);
        experiments <- as.character(unique(x$experiment_name));
        
        # merge each experiments minimum pvalue by id, Name and GS_Name
        by_gset <- Reduce(
            function( ... ) merge( ... , by=c("id", "Name", "GS_Name"),
            all=!FALSE), lapply(experiments, function(exp_name) {
                # for each experiment
                actRes <- x[ exp_name, ];
                
                enrichments <- actRes[, list(SEA_pval, GSEA_pval)];
                enrichments <- sapply(enrichments, as.numeric);
                colnames(enrichments) <- paste(exp_name, c("SEA", "GSEA"),
                                                sep="_");
                
                if (!by_method) {
                    # we will get the min pvalue as the result pvalue
                    enrichments <- suppressWarnings(
                        apply(enrichments, 1, min, na.rm=!FALSE));
                    enrichments[ is.infinite(enrichments) ] <- NA;
                    enrichments <- data.frame(enrichments);
                    colnames(enrichments) <- exp_name;
                }
                
                # create the new data table with all the info, plus the minimum 
                # pvalue of SEA/GSEA
                actRes <- data.table(actRes[, list(id, name, gene_set_name)],
                                        enrichments);
                colnames(actRes)[1:3] <- c("id", "Name", "GS_Name");
                setkey(actRes, id, Name, GS_Name);
                
                return(actRes);
            })
        )
        by_gset <- as.data.frame(by_gset);
        
        return(by_gset);
    }
)

setMethod(
    f="summary_from_all",
    signature=c("MIGSAres"),
    definition=function(x, by_method) {
        if (missing(by_method)) {
            by_method=FALSE;
        }
        
        by_gset <- summary_from_all(x@migsa_res_all, by_method);
        return(by_gset);
    }
)

############################ Exploratory functions

setGeneric(name="get_summary", def=function(migsa_obj) {
    standardGeneric("get_summary")
})

setMethod(
    f="get_summary",
    signature=c("MIGSAres"),
    definition=function(migsa_obj) {
        stopifnot(validObject(migsa_obj));
        
        enr_cutoff <- migsa_obj@enr_cutoff;
        res <- migsa_obj@migsa_res_summary;
        
        # when cutoff is 1 then it must enrich everything.
        enr_cutoff <- ifelse(enr_cutoff == 1, 1.1, enr_cutoff);
        
        if (!is.na(enr_cutoff)) {
            # if we have cutoff then lets return a logical (enriched or not)
            res[, -(1:3)] <- res[, -(1:3)] < enr_cutoff;
        }
        return(res);
    }
)

# following was made for giving different colors discriminating enrichment by
# SEA, GSEA, both
# setMethod(
#   f="migsaHeatmap",
#   signature=c("MIGSAres"),
#   definition=function(x, flevel=0, penrich=0, col.dist="jaccard",
#                       row.dist=col.dist, ...) {
#     require(gplots);
#     require(vegan);
#     allRes_byMeth <- by.gene.sets(x, by_method=!F);
#     allRes <- do.call(cbind, lapply(seq(3, ncol(allRes_byMeth), by=2), 
#         function(i) {
#       enrichments <- cbind(allRes_byMeth[,i:(i+1)]);
#       anyEnr <- apply(enrichments, 1, function(y) {
#         enr <- NA;
#         if (sum(is.na(y)) == 2) {
#           enr <- NA; # both results are NA
#         } else {
#           y[is.na(y)] <- F;
#           enr <- "none";
#           if (any(y)) {
#             enr <- "GSEA";
#             if (all(y)) {
#               enr <- "both";
#             } else if (y[[1]]) {
#               enr <- "SEA";
#             }
#           }
#         }
#         return(enr);
#       })
#       return(anyEnr);
#     }));
#     
#     # terms filtering by flevel and penrich
#     keepRows <- rowSums(allRes!="none", na.rm=!F)/ncol(allRes) >= flevel;
#     keepCols <- colSums(allRes!="none", na.rm=!F)/nrow(allRes) >= penrich;
#     allRes <- allRes[ keepRows, keepCols ];
#     allRes_byMeth <- allRes_byMeth[ keepRows, keepCols ];
#     
#     # however rows with no enrichment by any method are deleted
#     keepRows <- rowSums(allRes!="none", na.rm=!F) > 0;
#     allRes <- allRes[keepRows,];
#     allRes_byMeth <- allRes_byMeth[keepRows,];
#     rm(keepRows); rm(keepCols);
#     
#     # get color for each different database
#     rowColors <- as.character(allRes_byMeth$GS_Name);
#     pal <- rainbow(length(unique(rowColors)));
#     names(pal) <- unique(rowColors);
#     rowColors <- pal[rowColors]; rm(pal);
#     
#     # jaccard clustering per column (database)
#     numRes <- apply(allRes!="none", 2, as.numeric);
#     
#     dd.s <- suppressWarnings(vegdist(t(numRes), col.dist, na.rm=!F));
#     dd.s[is.na(dd.s)] <- 0;
#     h.s <- hclust(dd.s, method="average");
#     
#     # jaccard clustering per row (term)
#     ddr.s <- suppressWarnings(vegdist(numRes, row.dist, na.rm=!F));
#     ddr.s[is.na(ddr.s)] <- 0;
#     hr.s <- hclust(ddr.s, method="average");
#     
#     # give the colors to each slot of the migsaHeatmap
#     colors <- character();
#     if (any(is.na(allRes))) {
#       numRes[is.na(allRes)] <- 0; colors <- c(colors, "black");
#     }
#     if (any(allRes == "none")) {
#       numRes[allRes == "none"] <-  1; colors <- c(colors, "white");
#     }
#     if (any(allRes == "SEA")) {
#       numRes[allRes == "SEA"]  <-  2; colors <- c(colors, "red");
#     }
#     if (any(allRes == "GSEA")) {
#       numRes[allRes == "GSEA"] <- 3; colors <- c(colors, "green");
#     }
#     if (any(allRes == "both")) {
#       numRes[allRes == "both"] <- 4; colors <- c(colors, "violetred4");
#     }
#     
#     heatmap.2(as.matrix(numRes), Rowv=as.dendrogram(hr.s), trace="none",
#               labRow=rep("", nrow(numRes)), Colv=as.dendrogram(h.s),
#               colsep=1:(ncol(numRes)-1), sepwidth=c(0.025, 0.025), key=F,
#               RowSideColors=rowColors, breaks=length(colors)+1, col=colors,
#                 ...);
#     
#     # return selected terms and the clustering
#     return(list(row.dendrogram=as.dendrogram(hr.s),
#                 col.dendrogram=as.dendrogram(h.s)));
#   }
# )
# 

setGeneric(name="setDefaultEnrCutoff", def=function(object) {
    standardGeneric("setDefaultEnrCutoff")
})

# if there is not enr_cutoff set then it puts the default one (0.01)
setMethod(f="setDefaultEnrCutoff", signature="MIGSAres",
    definition=function(object) {
        if (is.na(enrCutoff(object))) {
            warning("No enrichment cutoff set. Using 0.01");
            enrCutoff(object) <- 0.01;
        }
        return(object)
    }
)

setGeneric(name="enrCutoff", def=function(object) {
    standardGeneric("enrCutoff")
})

setMethod(f="enrCutoff", signature="MIGSAres",
    definition=function(object) {
        return(object@enr_cutoff)
    }
)

setGeneric(name="enrCutoff<-", def=function(object, value) {
    standardGeneric("enrCutoff<-")
})

setReplaceMethod(f="enrCutoff", signature="MIGSAres",
    definition=function(object, value) {
        object@enr_cutoff <- value
        validObject(object)
        return(object)
    }
)

setGeneric(name="migsaResAll", def=function(object) {
    standardGeneric("migsaResAll")
})

setMethod(f="migsaResAll", signature="MIGSAres",
    definition=function(object) {
        return(object@migsa_res_all)
    }
)
