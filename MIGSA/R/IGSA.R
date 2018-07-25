setGeneric(name="IGSA", def=function(igsaInput) {
    standardGeneric("IGSA")
})

#'@importFrom futile.logger flog.info
#'@include DEnricher.R
#'@include IGSAinput.R
#'@include IGSAinput-class.R
#'@include IGSAinput-getterSetters.R
#'@include IGSAres.R
#'@include MGSZ.R
setMethod(
    f="IGSA",
    signature=c("IGSAinput"),
    definition=function(igsaInput) {
        flog.info("*************************************");
        flog.info(paste(name(igsaInput), ": Starting IGSA analysis."));
        
        fit_options <- fitOptions(igsaInput);
        expr_data   <- exprData(igsaInput);
        
        actGeneSets <- geneSetsList(igsaInput);
        if (length(actGeneSets) == 0) {
            stop("No gene sets provided.");
        }
        
        # merging all gene sets in order to run just one mGSZ and one dEnricher,
        # its much faster, but after we will have to split these results
        merged_gene_sets <- merge_gene_sets(actGeneSets);
        flog.info(paste(length(merged_gene_sets), "Gene Sets."));
        
        if (is.null(seaParams(igsaInput)) && is.null(gseaParams(igsaInput))) {
            stop("SEAparams and GSEAparams are both NULL");
        }
        
        deRes <- SEAres();
        if (!is.null(seaParams(igsaInput))) {
            flog.info(paste(name(igsaInput), ": dEnricher starting."));
            # run DEnricher
            deRes <- DEnricher(seaParams(igsaInput), expr_data, 
                            fit_options, merged_gene_sets);
            flog.info(paste(name(igsaInput), ": dEnricher finnished."));
        }
        
        mgszRes <- GSEAres();
        if (!is.null(gseaParams(igsaInput))) {
            flog.info(paste(name(igsaInput), ": mGSZ starting."));
            # run MGSZ
            mgszRes <- MGSZ(gseaParams(igsaInput), expr_data, fit_options,
                            merged_gene_sets);
            flog.info(paste(name(igsaInput), ": mGSZ finnished."));
        }
        
        # splitting all results. it is a list of GenesetsRes objects
        splitted_res <- split_results(deRes, mgszRes, actGeneSets);
        
        genes_rank <- data.frame();
        
        if (ncol(expr_data) > 0) {
            # get the genes rank to give more information in the results object
            genes_rank <- get_fit(expr_data, fit_options, SEAparams());
            genes_rank <- data.frame(geneID=rownames(genes_rank),
                                        rank=genes_rank$t);
            colnames(genes_rank)[2] <- name(igsaInput);
        }
        
        # create the IGSAres object
        igsaRes <- IGSAres(name=name(igsaInput), gene_sets_res=splitted_res,
                            genes_rank=genes_rank);
        
        flog.info(paste(name(igsaInput), ": IGSA analysis ended."));
        
        return(igsaRes);
    }
)

setGeneric(name="merge_gene_sets", def=function(geneSetsList) {
    standardGeneric("merge_gene_sets")
})

#'@importClassesFrom GSEABase GeneSetCollection
#'@importFrom GSEABase GeneSetCollection setName setName<- setIdentifier 
#'setIdentifier<-
# input: list of Genesets
# output: a Genesets object
setMethod(
    f="merge_gene_sets",
    signature=c("list"),
    definition=function(geneSetsList) {
        stopifnot(all(unlist(lapply(geneSetsList, function(x)
                    is(x, "GeneSetCollection")))));
        
        # each individual gene set id will be its gene set name "_MIGSA_"
        # its id. e.g. BP gene set GO:10000 will be BP_MIGSA_GO:10000
        # Not the best idea haha
        merged <- Reduce(function(...) append(...),
            lapply(names(geneSetsList), function(act_name) {
                gene_sets <- geneSetsList[[act_name]];
                act_gene_sets <- lapply(gene_sets, function(gene_set) {
                    new_gene_set <- gene_set;
                    
                    # I have to save and re copy it because if no it overwrites
                    #  ID, I think it is GSEABase bug
                    old_set_id <- setIdentifier(new_gene_set);
                    
                    new_set_name <- paste(act_name, "MIGSA",
                        setName(new_gene_set), sep="_");
                    setName(new_gene_set) <- new_set_name;
                    
                    setIdentifier(new_gene_set) <- old_set_id;
                    
                    return(new_gene_set);
                })
            })
        )
        return(GeneSetCollection(merged));
    }
)

setGeneric(name="split_results",
    def=function(sea_res, gsea_res, geneSetsList) {
    standardGeneric("split_results")
})

#'@importClassesFrom GSEABase GOCollection
#'@importFrom GSEABase collectionType
#'@include GenesetRes.R
#'@include GenesetsRes.R
#'@include GSEAres.R
#'@include SEAres.R
# input geneSetsList is a list of GeneSetCollection
setMethod(
    f="split_results",
    signature=c("SEAres", "GSEAres", "list"),
    definition=function(sea_res, gsea_res, geneSetsList) {
        
        # So now lets split the merged results!!
        # and return to each its additional info they had (as is_GO)
        
        splitted_res <- lapply(names(geneSetsList), function(act_name) {
            act_GS <- geneSetsList[[act_name]];
            act_is_GO <- any(unlist(lapply(act_GS, function(x) 
                is(collectionType(x), "GOCollection"))));
            # grep the ones which start with gene set name followed by "_MIGSA_"
            act_pattern <- paste("^", act_name, "_MIGSA_", sep="");
            act_sea_res <- lapply(gene_sets_res(sea_res), function(act_res) {
                if (grepl(act_pattern, id(act_res))) {
                    id(act_res) <- gsub(act_pattern, "", id(act_res));
                    getName(act_res) <- gsub(act_pattern, "", getName(act_res));
                    return(act_res);
                } else {
                    return(NA);
                }
            })
            act_sea_res <- SEAres(
                            gene_sets_res=act_sea_res[!is.na(act_sea_res)]);
            
            act_gsea_res <- lapply(gene_sets_res(gsea_res), function(act_res) {
                if (grepl(act_pattern, id(act_res))) {
                    id(act_res) <- gsub(act_pattern, "", id(act_res));
                    getName(act_res) <- gsub(act_pattern, "", getName(act_res));
                    return(act_res);
                } else {
                    return(NA);
                }
            })
            act_gsea_res <- GSEAres(
                            gene_sets_res=act_gsea_res[!is.na(act_gsea_res)]);
            
            return(GenesetsRes(gene_sets_name=act_name, is_GO=act_is_GO,
                                sea_res=act_sea_res, gsea_res=act_gsea_res));
        })
        return(splitted_res);
    }
)
