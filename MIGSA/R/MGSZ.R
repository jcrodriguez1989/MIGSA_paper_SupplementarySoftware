setGeneric(name="MGSZ", def=function(params, M, fit_options, genesets) {
    standardGeneric("MGSZ")
})

# params <- gseaParams(igsaInput); M <- exprData(igsaInput);
# genesets <- merged_gene_sets;
#'@importFrom GSEABase GeneSetCollection setIdentifier setName
#'@include FitOptions.R
#'@include GenesetRes.R
#'@include Genesets.R
#'@include GSEAparams.R
#'@include GSEAres.R
#'@include IGSAinput.R
setMethod(
    f="MGSZ",
    signature=c("GSEAparams", "ExprData", "FitOptions", "GeneSetCollection"),
    definition=function(params, M, fit_options, genesets) {
        # mGSZ uses the gene sets as a list
        mgsz.gene.sets <- asList(genesets);
        
        # check the ranking function depending if we use voom or not
        if(is(M, "DGEList")) {
            rankFunc <- voomLimaRank;
        } else {
            rankFunc <- mGszEbayes;
        }
        
        # run faster version of mGSZ
        mgszRes <- MIGSA_mGSZ(M, fit_options, mgsz.gene.sets, rankFunc, params);
        
        # convert mGSZ results (data.frame) to GenesetRes
        mgszRes <- lapply(genesets, function(actGset) {
            mgszGSres <- mgszRes[ mgszRes$gene.sets == setName(actGset), ];
            if (nrow(mgszGSres) == 0) {
                pval <- as.numeric(NA);
                actScore <- as.numeric(NA);
                actImpGenes <- "";
            } else {
                pval <- mgszGSres$pvalue;
                actScore <- mgszGSres$mGszScore;
                actImpGenes <- as.character(mgszGSres$impGenes);
                # todo: put better genes separator
                actImpGenes <- strsplit(actImpGenes, ", ")[[1]];
            }
            
            # We put name as setIdentifier and id as setName because we 
            # disagree with GSEABase's usage as they use Names as unique, and 
            # identifiers not.
            actRes <- GenesetRes(id=setName(actGset), 
                                name=setIdentifier(actGset), score=actScore, 
                                pvalue=pval, genes=geneIds(actGset),
                                enriching_genes=actImpGenes);
        })
        
        # this is a GSEAres
        gseaRes <- GSEAres(gene_sets_res=mgszRes);
        
        return(gseaRes);
    }
)

setGeneric(name="MIGSA_mGSZ",
    def=function(exprData, fitOptions, gSets, rankFunction, params) {
    standardGeneric("MIGSA_mGSZ")
})

# gene.sets <- mgsz.gene.sets; rankFunction <- rankFunc
# M <- igsaInput@expr_data; fit_options <- igsaInput@fit_options;
#'@importClassesFrom edgeR DGEList
#'@importFrom BiocParallel bplapply bpparam
#'@importFrom data.table data.table
#'@importFrom edgeR [.DGEList
#'@importFrom futile.logger flog.debug flog.info
#'@importFrom ismev gum.fit
#'@importFrom matrixStats colSds rowMedians colVars
#'@include FitOptions.R
#'@include GSEAparams.R
#'@include IGSAinput.R
setMethod(
    f="MIGSA_mGSZ",
    signature=c("ExprData", "FitOptions", "list", "function", "GSEAparams"),
    definition=function(exprData, fitOptions, gSets, rankFunction, params) {
        if (is(exprData, "DGEList")) {
            stopifnot(all(rownames(exprData$samples) == 
                                colnames(exprData$counts)));
        }
        
        # filtering inputs, as required
        filteredInputs <- filterInputs(exprData, gSets, minSz(params));
        exprData <- filteredInputs$exprData;
        gSets <- filteredInputs$gSets;
        rm(filteredInputs);
        
        # get gene rankings (also for permuted data)
        nPerm <- perm_number(params);
        preVar <- pv(params);
        wgt2 <- w2(params);
        varConstant <- vc(params);
        
        rankings <- getRankings(exprData, fitOptions, nPerm, rankFunction);
        rownames(rankings) <- rownames(exprData);
        rm(nPerm);
        
        # calculate normalizing values
        setSizes <- unlist(lapply(gSets, length));
        uniqSizes <- unique(setSizes);
        normValues <- countHygeVarMean(nrow(rankings), uniqSizes);
        normMean <- normValues$mean; colnames(normMean) <- uniqSizes;
        normVar <- normValues$var; colnames(normVar) <- uniqSizes;
        rm(normValues);
        
        probSum <- do.call(cbind, lapply(seq_along(uniqSizes), function(j) {
            countProbSum(nrow(rankings), normMean[, j], normVar[, j]);
        }));
        colnames(probSum) <- uniqSizes;
        normFactors <- list(normMean=normMean, normVar=normVar, 
                                                probSum=probSum);
        
        # get extra info for real data
        realRank <- rankings[,1];
        realRank <- sort(realRank, decreasing=!FALSE);
        sVMC_dec <- sumVarMeanCalc(realRank, preVar, normFactors);
        sVMC_inc <- sumVarMeanCalc(rev(realRank), preVar, normFactors);
        zVars_dec <- zVarCalc(sVMC_dec$Z_vars, wgt2, varConstant);
        zVars_inc <- zVarCalc(sVMC_inc$Z_vars, wgt2, varConstant);
        
        enrichScores <- do.call(c, bplapply(uniqSizes, function(actSize) {
    #             actSize <- uniqSizes[[1]];
            actGsets <- gSets[setSizes == actSize];
            actSize_c <- as.character(actSize);
            
            zM_d <- sVMC_dec$Z_means[, actSize_c];
            zM_i <- sVMC_inc$Z_means[, actSize_c];
            zV_d <- zVars_dec[, actSize_c];
            zV_i <- zVars_inc[, actSize_c];
            
            actESs <- lapply(actGsets, function(actGset) {
                    diffScores <- getEnrichScore_c(realRank, actGset, 
                                                zM_d, zM_i, zV_d, zV_i);
                    # real data
                    maxI <- which.max(abs(diffScores));
                    rawES <- diffScores[[maxI]];
                    nGenes <- length(realRank);
                    
                    if (maxI > nGenes) {
                        impGenes <- names(diffScores)[(nGenes+1):maxI];
                    } else {
                        impGenes <- names(diffScores)[1:maxI];
                    }
                    impGenes <- intersect(impGenes, actGset);
                    actES <- abs(rawES);
                    
                    return(list(actES=actES, rawES=rawES, impGenes=impGenes));
                });
            names(actESs) <- names(actGsets);
            return(actESs);
        }))
        rm(sVMC_dec); rm(sVMC_inc); rm(zVars_dec); rm(zVars_inc);
        
        realESs <- do.call(c, lapply(enrichScores, function(x) x$actES));
        rawESs <- do.call(c, lapply(enrichScores, function(x) x$rawES));
        impGenes <- do.call(c, lapply(enrichScores, function(x) 
                                            paste(x$impGenes, collapse=', ')));
        
        permESs <- do.call(cbind, bplapply(seq_len(ncol(rankings)-1)+1,
            function(i) {
    #         i <- 2;
            ranking <- rankings[,i];
            ranking <- sort(ranking, decreasing=!FALSE);
            sVMC_dec <- sumVarMeanCalc(ranking, preVar, normFactors);
            sVMC_inc <- sumVarMeanCalc(rev(ranking), preVar, normFactors);
            
            zVars_dec <- zVarCalc(sVMC_dec$Z_vars, wgt2, varConstant);
            zVars_inc <- zVarCalc(sVMC_inc$Z_vars, wgt2, varConstant);
            
            # as norm factors are separated by gene set sizes, then in order to
            # consume less ram, lets separate by sizes.
            # In some cases it will not use the best of cores, but its for ram!
            enrichScores <- do.call(c, lapply(uniqSizes, 
                function(actSize, sVMC_dec, sVMC_inc, zVars_dec, zVars_inc) {
    #             actSize <- uniqSizes[[1]];
                actGsets <- gSets[setSizes == actSize];
                actSize_c <- as.character(actSize);
                
                act_z_means_dec <- sVMC_dec$Z_means[, actSize_c];
                act_z_means_inc <- sVMC_inc$Z_means[, actSize_c];
                act_zVars_dec <- zVars_dec[, actSize_c];
                act_zVars_inc <- zVars_inc[, actSize_c];
#                 rm(sVMC_dec); rm(sVMC_inc); rm(zVars_dec); rm(zVars_inc);
                
                actESs <- do.call(c, lapply(actGsets,
                    function(actGset, zM_d, zM_i, zV_d, zV_i) {
                        actES <- max(abs(getEnrichScore_c(ranking, actGset, 
                                                    zM_d, zM_i, zV_d, zV_i)));
                        return(actES);
                    },
                    zM_d=act_z_means_dec, zM_i=act_z_means_inc,
                    zV_d=act_zVars_dec, zV_i=act_zVars_inc
                ));
                names(actESs) <- names(actGsets);
                return(actESs);
            }, sVMC_dec=sVMC_dec, sVMC_inc=sVMC_inc, 
                zVars_dec=zVars_dec, zVars_inc=zVars_inc))
            return(enrichScores);
        }));
        
        stopifnot(all(names(realESs) == rownames(permESs)));
        pvals <- getPvalues(realESs, t(permESs));
        names(pvals) <- names(realESs);
        
        result <- data.frame(gene.sets=names(pvals), pvalue=pvals, 
                                mGszScore=rawESs, impGenes=impGenes);
        result <- result[order(pvals),];
        
        return(result);
    }
)

#'@importFrom limma treat contrasts.fit lmFit
#'@include FitOptions.R
mGszEbayes <- function(exprMatrix, fit_options) {
#     flog.info("Using ebayes");
    design <- designMatrix(fit_options);
    contrast <- contrast(fit_options);
    fit1 <- lmFit(exprMatrix, design);
    fit2 <- treat(contrasts.fit(fit1, contrast));
    
    res <- fit2$t;
    return(res);
}

## voom + limma
#'@importFrom limma treat contrasts.fit lmFit voom
#'@include FitOptions.R
voomLimaRank <- function(exprMatrix, fit_options) {
#     flog.info("Using voom+limma");
    design <- designMatrix(fit_options);
    contrast <- contrast(fit_options);
    newExpr <- voom(exprMatrix, design);
    
    # Adjust the model
    fit1 <- lmFit(newExpr, design);
    fit2 <- treat(contrasts.fit(fit1, contrast));
    
    res <- fit2$t;
    return(res);
}

filterInputs <- function(exprData, gSets, minGsetSize) {
    # keep only genes in exprData and gene sets
    exprDataGenes <- rownames(exprData);
    geneSetsGenes <- do.call(c, gSets);
    commonGenes <- intersect(exprDataGenes, geneSetsGenes);
    
    exprData <- exprData[commonGenes,];
    gSets <- lapply(gSets, intersect, commonGenes);
    
    # remove gene sets with less than minSz genes (detected in the experiment)
    gSets <- gSets[lapply(gSets, length) > minGsetSize];
    return(list(exprData=exprData, gSets=gSets));
}

getRankings <- function(exprData, fitOptions, nPerm, rankFunction) {
    # generate permutations
    perms <- replicate(nPerm, sample(1:ncol(exprData), replace=FALSE));
    conds <- col_data(fitOptions);
    permsConds <- do.call(cbind, lapply(seq_len(ncol(perms)), function(i) {
        as.character(conds[perms[,i],]);
    }))
    perms <- t(perms[,!duplicated(t(permsConds))]);
    rm(permsConds);
    
    nPerm <- nrow(perms);
    flog.info(paste("Number of unique permutations:", nPerm));
    
    rankings <- matrix(0, nrow=nrow(exprData), ncol=nPerm+1);
    # gene scores for real data
    rankings[,1] <- rankFunction(exprData, fitOptions);
    
    flog.info(paste("Getting ranking at cores:", bpparam()$workers));
    
    # gene scores for permuted data
    # I am passing MIGSA's not exported functions to bplapply to avoid
    # SnowParam environment errors
    rankings[,-1] <- do.call(cbind,
    bplapply(seq_len(nPerm), function(i, designMatrix) {
        # modify the design matrix using the order given by the permutation
        permDesign <- designMatrix(fitOptions)[perms[i,],];
        newFitOptions <- fitOptions;
        newFitOptions@design_matrix <- permDesign;
        
        actRank <- rankFunction(exprData, newFitOptions);
        return(actRank);
    }, designMatrix=designMatrix));
    
    # todo: maybe the best alternative would be delete the gene from
    # everywhere
    rankings[is.na(rankings)] <- mean(rankings, na.rm=TRUE);
    # genes which are NA are given the mean value of all genes, 
    # so they dont contribute to enrichment score.
    
    return(rankings);
}

countHygeVarMean <- function (M, K) {
    N <- matrix(rep(seq_len(M), length(K)), ncol=length(K));
    K <- matrix(rep(K, M), byrow=TRUE, ncol=length(K));
    mean <- N * K/M;
    var <- N * (K * (M - K)/M^2) * ((M - N)/(M - 1));
    return(list(mean=mean, var=var));
}

countProbSum <- function (M, hyge.mean, hyge.var) {
    tulos <- matrix(0, nrow=length(hyge.mean));
    N_tab <- matrix(c(2:M), byrow=FALSE);
    tulos[1] <- hyge.mean[1];
    tulos[2:M] <- N_tab/(N_tab-1) * hyge.mean[2:M] - (hyge.mean[2:M]^2 + 
        hyge.var[2:M])/(N_tab-1);
    return(tulos);
}

sumVarMeanCalc <- function(ranking, preVar, normFactors) {
#     preVar <- pv(params);
    divider <- seq_along(ranking);
    mean_table <- cumsum(ranking)/divider;
    mean_table_sq <- mean_table^2;
    
    # not sure if it is +2*preVar or just +preVar , in mGSZ it does +preVar
    # and again +preVar
    var_table <- cumsum(ranking^2)/divider - (mean_table)^2 + 2*preVar;
    
    z_means <- do.call(cbind, lapply(seq_len(ncol(normFactors$normMean)),
    function(j) {
        mean_table * (2 * normFactors$normMean[,j] - divider);
    })); colnames(z_means) <- colnames(normFactors$normMean);
    
    z_vars <- do.call(cbind, lapply(seq_len(ncol(normFactors$normVar)), 
    function(j) {
        4 * (var_table * normFactors$probSum[,j] + mean_table_sq * 
                                                    normFactors$normVar[,j]);
    })); colnames(z_vars) <- colnames(normFactors$normVar);
    
    return(list(Z_means=z_means, Z_vars=z_vars));
}

zVarCalc <- function(z_vars, wgt2, varConstant) {
    medianPart <- matrix(matrixStats::rowMedians(t(z_vars)) * 
                    wgt2, nrow=nrow(z_vars), ncol=ncol(z_vars), byrow=!FALSE);
    z_var <- (z_vars + medianPart + varConstant)^0.5;
    return(z_var);
}

getEnrichScore <- function(ranking, actGset, zM_d, zM_i, zV_d, zV_i) {
    startVal <- 5; # hard coded
    nMG <- !names(ranking) %in% actGset;
    ranking[nMG] <- -ranking[nMG];
    
    es_dec <- cumsum(ranking) - zM_d;
    es_inc <- cumsum(rev(ranking)) - zM_i;
    
    es_dec[1:startVal] <- 0; es_inc[1:startVal] <- 0;
    
    es_dec <- es_dec/zV_d;
    es_inc <- es_inc/zV_i;
    
    return(c(es_dec, es_inc));
#     return(max(abs(c(es_dec, es_inc))));
}

#'@importFrom compiler cmpfun
getEnrichScore_c <- cmpfun(getEnrichScore);

logEVcdf <- function (x, fitParams) {
    mu <- fitParams[[1]];
    sigma <- fitParams[[2]];
    x <- (mu - x)/sigma;
    out <- rep(0, length(x));
    if (!(is.vector(x))) {
        out <- matrix(out, nrow(x), ncol(x));
    }
    po1 <- which(x < 5);
    out[po1] <- -log(1 - exp(-exp(x[po1])));
    x <- x[-po1];
    out[-po1] <- -x + exp(x)/2 - exp(2 * x)/24 + exp(4 * x)/2880;
    out <- out/log(10);
    out <- 10^(-out);
    return(out);
}

# realES <- allESs[,1]; permESs <- t(allESs[,-1])
getPvalues <- function(realES, permESs) {
#     pos.data <- realES; perm.data <- permESs;
    # some normalization
    mean.prof <- colMeans(permESs);
    std.prof <- matrixStats::colSds(permESs);
    std.prof[std.prof == 0] <- 0.1;
    realES <- (realES - mean.prof)/std.prof;
    permESs <- do.call(cbind, lapply(seq_len(ncol(permESs)), function(k) {
        (permESs[, k] - mean.prof[k])/std.prof[k];
    }))
    rm(mean.prof); rm(std.prof);

    col.ind <- colVars(permESs) > 0;
    permESs <- permESs[, col.ind];
    realES <- realES[col.ind];
    ev.p.val.class <- rep(0, length(realES))
    
    pvals <- do.call(c, lapply(seq_along(realES), function(k) {
        ev.param.class <- ismev::gum.fit(as.vector(permESs[, k]), 
            show=FALSE)$mle
        logEVcdf(realES[[k]], ev.param.class);
    }))
    return(pvals);
}
