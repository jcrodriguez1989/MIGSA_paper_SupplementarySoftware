#'Get or edit MIGSAres with/by genes that contributed to enrichment
#'
#'\code{genesInSets} returns a data.frame with gene sets as rows, genes as 
#'columns, and as value the number of experiments in which each gene 
#'contributed to enrich each gene set. 
#'If it was enriched only by SEA then it returns the genes that contributed in 
#'SEA. If it was enriched only by GSEA then it returns the genes that 
#'contributed in GSEA. If it was enriched by both then it returns the genes 
#'that contributed in both.
#'\code{filterByGenes} returns a MIGSAres object with only the gene sets which 
#'resulted enriched by at least one gene from a provided list.
#'
#'@param migsaRes MIGSAres object.
#'@param genes character vector of the interest genes for MIGSAres filtering.
#'
#'@return If genesInSets: a data.frame with the number of experiments in which 
#'each gene contributed to enrich each gene set. If filterByGenes: A MIGSAres 
#'object containing only the gene sets in which provided genes contributed to 
#'enrichment.
#'
#'@docType methods
#'@name MIGSAres-genes
#'@rdname MIGSAres-genesManipulation
#'@seealso \code{\link{setEnrCutoff}}
#'
#'@examples
#'data(migsaRes);
#'
#'###### genesInSets
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'gInSets <- genesInSets(migsaResWCoff);
#'class(gInSets); # matrix
#'
#'## Now we can do stuff as check which genes enriched a gene set in all (two) 
#'## experiments.
#'gInSets[rowSums(gInSets==2) > 0, colSums(gInSets==2) > 0];
#'
#'###### filterByGenes
#'## Suppose we are interested in studying these genes:
#'intGenes <- c("g10", "g91", "g388", "g742", "g874");
#'migsaResIntGenes <- filterByGenes(migsaResWCoff, intGenes);
#'
#'## Now in migsaResIntGenes we have the MIGSA results of the gene sets in 
#'## which at least one gene of our list contributed to enrich.
#'migsaResIntGenes
#'
setGeneric(name="MIGSAres-genes", def=function(migsaRes) {
    standardGeneric("MIGSAres-genes")
})

#'@name genesInSets
#'@inheritParams MIGSAres-genes
#'@rdname MIGSAres-genesManipulation
#'@aliases genesInSets,MIGSAres-method
#'@exportMethod genesInSets
#'
setGeneric(name="genesInSets", def=function(migsaRes) {
    standardGeneric("genesInSets")
})

#'@inheritParams MIGSAres-genes
#'@rdname MIGSAres-genesManipulation
#'@aliases genesInSets,MIGSAres-method
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="genesInSets",
    signature=c("MIGSAres"),
    definition=function(migsaRes) {
        stopifnot(validObject(migsaRes));
        
        data_frame_all <- migsaRes@migsa_res_all;
        setkey(data_frame_all, gene_set_name, id);
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        enr_cutoff <- enrCutoff(migsaRes);
        
        # when cutoff is 1 then it must enrich everything.
        enr_cutoff <- ifelse(enr_cutoff == 1, 1.1, enr_cutoff);
        
        # get the gene sets name and id
        gsets <- apply(unique(data_frame_all[,list(gene_set_name, id)]),2,
                        as.character);
        res <- lapply(split(gsets, seq(nrow(gsets))), function(gset) {
            # for each gene set get the data
            actData <- data_frame_all[ gene_set_name == gset[[1]] &
                                        id == gset[[2]], ]
            enr_genes <- apply(actData, 1, function(y) {
                # for each gene set data
                gsea_enr_genes <- NULL;
                sea_enr_genes <- NULL;
                
                gsea_pval <- as.numeric(y[["GSEA_pval"]]);
                sea_pval  <- as.numeric(y[["SEA_pval"]]);
                
                if (!is.na(gsea_pval) && gsea_pval < enr_cutoff) {
                    # if was enriched by gsea get its enriching genes
                    gsea_enr_genes <- gsub(" ", "", unlist(
                        strsplit(as.character(y[["GSEA_enriching_genes"]]),
                                ",")));
                }
                
                if (!is.na(sea_pval) && sea_pval < enr_cutoff) {
                    # if was enriched by sea get its enriching genes
                    sea_enr_genes <- gsub(" ", "", unlist(
                        strsplit(as.character(y[["SEA_enriching_genes"]]),
                                ",")));
                }
                
                enr_genes <- union(gsea_enr_genes, sea_enr_genes);
                return(enr_genes);
            })
            
            if (is.list(enr_genes)) {
                enr_genes <- do.call(c, enr_genes);
            } else {
                enr_genes <- c(enr_genes);
            }
            
            res <- NA;
            if (length(enr_genes) > 0) {
                # if there was any gene then sort them
                res <- sort(table(enr_genes), decreasing=!FALSE);
            }
            
            return(res);
        })
        names(res) <- paste(gsets[,1], gsets[,2], sep="_");
        
        gsets <- unique(unlist(lapply(res, names)));
        # lets turn these results (list) into a data frame
        finalRes <- do.call(rbind, lapply(names(res), function(y) {
            actRes <- res[[y]];
            newRes <- rep(0, length(gsets));
            if (!is.na(actRes[[1]])) {
                names(newRes) <- gsets;
                newRes[ names(actRes) ] <- actRes;
            }
            return(newRes);
        }))
        rownames(finalRes) <- names(res);
        
        return(finalRes);
    }
)

#'@name filterByGenes
#'@inheritParams MIGSAres-genes
#'@rdname MIGSAres-genesManipulation
#'@aliases filterByGenes,MIGSAres,character-method
#'@exportMethod filterByGenes
#'
setGeneric(name="filterByGenes", def=function(migsaRes, genes) {
    standardGeneric("filterByGenes")
})

#'@inheritParams MIGSAres-genes
#'@rdname MIGSAres-genesManipulation
#'@aliases filterByGenes,MIGSAres,character-method
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="filterByGenes",
    signature=c("MIGSAres", "character"),
    definition=function(migsaRes, genes) {
        stopifnot(validObject(migsaRes));
        stopifnot(length(genes) > 0);
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # lets get just gene sets enriched in at least one dataset
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=!FALSE) > 0, ];
        
        if (!is(actRes, "MIGSAres")) {
            warning("No enriched gene set with used cutOff.");
            return(data.frame());
        }
        migsaRes_genes <- genesInSets(actRes);
        
        if (any(genes %in% colnames(migsaRes_genes))) {
            # if we have any of our interest genes in the enriching ones,
            # then remove the other genes
            migsaRes_genes <- migsaRes_genes[, genes, drop=FALSE];
            
            # keet just the gene sets enriched by at least one of these genes
            migsaRes_genes <- migsaRes_genes[
                        rowSums(migsaRes_genes, na.rm=!FALSE)>0,, drop=FALSE];
            
            # complete gene sets names
            gsNames <- paste(migsaRes[,3,drop=TRUE],
                                migsaRes[,1,drop=TRUE], sep="_");
            
            # so return the subsetted MIGSAres object
            res <- migsaRes[ gsNames %in% rownames(migsaRes_genes),,
                                drop=FALSE];
            
            if (nrow(res) != nrow(migsaRes_genes)) {
                warning("Different number of rows after filtering by genes.");
            }
        } else {
            warning("Provided genes did not contribute to enrich any gene set."
                );
            res <- data.frame();
        }
        
        return(res);
    }
)
