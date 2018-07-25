#'Summary functions for some MIGSA classes
#'
#'R base summary overwritten functions to manipulate MIGSA objects.
#'
#'@param object SEAparams, GSEAparams, IGSAinput or MIGSAres object.
#'@param ... not in use.
#'
#'@return A summary of the object.
#'
#'@docType methods
#'@name summary
#'@rdname summaryFunctions
#'
#'@examples
#'## Lets get the summary of the default SEAparams object
#'seaParams <- SEAparams();
#'summary(seaParams);
#'
#'## Lets get the summary of the default GSEAparams object
#'gseaParams <- GSEAparams();
#'summary(gseaParams);
#'
#'## Lets create a basic valid IGSAinput object to get its summary.
#'## First create a expression matrix.
#'maData <- matrix(rnorm(10000),ncol=4);
#'rownames(maData) <- 1:nrow(maData); # It must have rownames (gene names).
#'maExprData <- new("MAList",list(M=maData));
#'
#'## Now lets create the FitOptions object.
#'myFOpts <- FitOptions(c("Cond1", "Cond1", "Cond2", "Cond2"));
#'
#'## And now we can create our IGSAinput ready for MIGSA.
#'igsaInput <- IGSAinput(name="myIgsaInput", expr_data=maExprData, 
#'fit_options=myFOpts);
#'summary(igsaInput);
#'
#'## Now lets get the summary of out migsaRes data object.
#'data(migsaRes);
#'
#'### As enrichment cutoff is not set then we will get for each experiment the 
#'### number of enriched gene sets at different cutoff values.
#'summary(migsaRes);
#'
#'### Lets set the enrichment cutoff at 0.01
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.01);
#'
#'### Now as summary we will get the number of enriched gene sets per 
#'### experiment and their intersections.
#'summary(migsaResWCoff);
#'
#'@include SEAparams-class.R
#'@method summary SEAparams
#'@aliases summary,SEAparams-method
#'@export summary.SEAparams
summary.SEAparams <- function(object, ...) {
    stopifnot(validObject(object));
    
    br <- object@br;
    if (length(br) > 1) {
        br <- "UserDefined";
    }
    
    res <- c(object@treat_lfc,
            object@de_cutoff,
            object@adjust_method,
            length(object@de_genes),
            br);
    names(res) <- c("treat_lfc", "de_cutoff", "adjust_method", "#de_genes",
                                                                    "br");
    
    return(res);
}

#'@inheritParams summary
#'@rdname summaryFunctions
#'@include GSEAparams-class.R
#'@method summary GSEAparams
#'@aliases summary,GSEAparams-method
#'@export summary.GSEAparams
#'
summary.GSEAparams <- function(object, ...) {
    stopifnot(validObject(object));
    
    # not showing the other params, as they are from mGSZ, I dont know if 
    # anyone uses them
    res <- object@perm_number;
    names(res) <- "perm_number";
    
    return(res);
}

#'@inheritParams summary
#'@rdname summaryFunctions
#'@include IGSAinput-class.R
#'@method summary IGSAinput
#'@aliases summary,IGSAinput-method
#'@export summary.IGSAinput
#'
summary.IGSAinput <- function(object, ...) {
    validObject(object);
    
    deGenes <- getDEGenes(object);
    # number of samples of each contrast
    ctrst <- table(col_data(deGenes@fit_options));
    sea_params <- summary(deGenes@sea_params);
    gsea_params <- summary(deGenes@gsea_params);
    
    res <- c(deGenes@name,
            ncol(deGenes@expr_data),
            paste(names(ctrst), collapse="VS"),
            ctrst[[1]],
            ctrst[[2]],
            length(deGenes@gene_sets_list),
            nrow(deGenes@expr_data),
            sea_params,
            gsea_params,
            round(100*(as.numeric(sea_params[[4]]) / 
                nrow(deGenes@expr_data)), 2)
    );
    names(res) <- c("exp_name",
                    "#samples",
                    "contrast",
                    "#C1",
                    "#C2",
                    "#gene_sets",
                    "#genes",
                    names(sea_params),
                    names(gsea_params),
                    "%de_genes");
    return(res);
}

#'@inheritParams summary
#'@rdname summaryFunctions
#'@include MIGSAres-class.R
#'@importFrom futile.logger flog.info
#'@method summary MIGSAres
#'@aliases summary,MIGSAres-method
#'@export summary.MIGSAres
#'
summary.MIGSAres <- function(object, ...) {
    stopifnot(validObject(object));
    
    pvals <- object@migsa_res_summary[,-(1:3), drop=FALSE];
    
    if (is.na(object@enr_cutoff)) {
        # if there is no cutoff then give results with these three cutoffs
        res <- rbind(
            enr_at_0_01=colSums(pvals < 0.01, na.rm=!FALSE),
            enr_at_0_05=colSums(pvals < 0.05, na.rm=!FALSE),
            enr_at_0_1=colSums(pvals < 0.1, na.rm=!FALSE)
        )
    } else {
        # if we have a cutoff set then give some statistics
        consGsets <- table(rowSums(pvals < object@enr_cutoff,
            na.rm=!FALSE));
        invisible(lapply(names(consGsets), function(actName) {
            flog.info(paste("Gene sets enriched in", actName,
                "experiments:", consGsets[actName]));
        }))
        
        numExps <- ncol(pvals);
        # lets get the gene sets enriched between each pair of experiments
        enrInters <- do.call(rbind, 
        lapply(1:ncol(pvals), function(actExp1) {
            lapply(1:ncol(pvals), function(actExp2) {
                sum(
                    pvals[,actExp1] < object@enr_cutoff &
                    pvals[,actExp2] < object@enr_cutoff,
                    na.rm=!FALSE
                );
            })
        }))
        colnames(enrInters) <- colnames(pvals);
        rownames(enrInters) <- colnames(pvals);
        res <- list(consensusGeneSets=consGsets, 
                    enrichmentIntersections=enrInters);
    }
    
    return(res);
}
