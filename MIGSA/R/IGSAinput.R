#'@include IGSAinput-class.R
setGeneric(name="get_fit", def=function(M, fit_options, params) {
    standardGeneric("get_fit")
})

#'@importFrom limma contrasts.fit lmFit treat voom
#'@importFrom stats p.adjust
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include SEAparams.R
setMethod(f="get_fit",
    signature=c("ExprData", "FitOptions", "SEAparams"),
    definition=function(M, fit_options, params) {
        act_treat_lfc <- treat_lfc(params);
        act_design <- designMatrix(fit_options);
        act_contrast <- contrast(fit_options);
        act_adj_meth <- adjust_method(params);
        
        # apply voom
        if (is(M, "DGEList")) {
            M <- voom(M, design=act_design);
        }
        
        # Adjust the model
        fit <- lmFit(M, act_design);
        # treat correction
        fit2 <- treat(contrasts.fit(fit, act_contrast), lfc=act_treat_lfc);
        # Adjusted pvalues
        fit2$p.adjust <- apply(fit2$p.value, 2, p.adjust, method=act_adj_meth);
        
        return(fit2);
    }
)

setGeneric(name="igsaGetDEGenes",
    def=function(seaParams, exprData, fitOptions) {
    standardGeneric("igsaGetDEGenes")
})

#'@importFrom futile.logger flog.info
#'@include ExprData-class.R
#'@include FitOptions-class.R
#'@include IGSAinput-class.R
setMethod(
    f="igsaGetDEGenes",
    signature=c("SEAparams", "ExprData", "FitOptions"),
    definition=function(seaParams, exprData, fitOptions) {
        # get the fit
        act_fit <- get_fit(exprData, fitOptions, seaParams);
        
        de_coff <- de_cutoff(seaParams);
        # get the DE genes depending on the cutoff value
        dif <- act_fit$p.adjust[,,drop=FALSE] <= de_coff;
        dif <- unique(rownames(dif)[dif]);
        
        flog.info(paste("DE genes", length(dif), "of a total of", 
                    nrow(exprData), "(", 
                    round(length(dif)/nrow(exprData)*100,2), "%)"));
        
        return(dif);
    }
)
