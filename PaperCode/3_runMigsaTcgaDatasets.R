# alias R=R3.5.0;
# R
# stopifnot(R.Version()$version.string ==
#     'R Under development (unstable) (2017-12-13 r73907)');
setwd('~/NatureMIGSA/'); # or any directory where Supplementary Data was downloaded to

library('edgeR');
library('MIGSA');
library('pbcmc');

# Download TCGAcohorts.xz from
# https://figshare.com/articles/TCGAcohorts_xz/6741686
TCGA <- readRDS('TCGAcohorts.xz');

rnaSeq <- TCGA$RNAseq;
geneExpr <- TCGA$MA;
protExpr <- TCGA$prot;
permSubtypes <- TCGA$pbcmc;
rm(TCGA);

stopifnot(all(colnames(geneExpr) == colnames(protExpr)));
stopifnot(all(colnames(geneExpr) == colnames(rnaSeq)));
stopifnot(all(colnames(geneExpr) == rownames(permSubtypes)));

table(permSubtypes$Permuted);
#     Assigned Not Assigned    Ambiguous 
#           68           15           14

permSubtypes$Permuted[ permSubtypes$Permuted == 'Ambiguous' ] <- 'Assigned';
# we will use Ambiguous subjects as if they were correctly assigned with their 
# original genefu classification

subtypes <- as.character(permSubtypes$PAM50);
names(subtypes) <- rownames(permSubtypes);
sum(permSubtypes$Permuted == 'Not Assigned');
# [1] 15
subtypes <- subtypes[permSubtypes$Permuted != 'Not Assigned'];

table(subtypes);
# Basal  Her2  LumA  LumB 
#    25    25    12    20

contrasts <- combn(c('Basal', 'LumA', 'Her2', 'LumB'), 2);
# actCtrst <- contrasts[,1];

# filtrate tcga matrices
allExperiments <- apply(contrasts, 2, function(actCtrst) {
    actC1subj <- names(subtypes)[ subtypes == actCtrst[[1]] ]; # cond1 subjects
    actC2subj <- names(subtypes)[ subtypes == actCtrst[[2]] ]; # cond2 subjects
    
    print(c(actCtrst, length(actC1subj), length(actC2subj)));
    
    actGeneExpr <- geneExpr[, c(actC1subj, actC2subj) ];
    actRnaSeq <- rnaSeq[, c(actC1subj, actC2subj) ];
    actProtExpr <- protExpr[, c(actC1subj, actC2subj) ];
    actSubtypes <- as.character(subtypes[ c(actC1subj, actC2subj) ]);
    
    stopifnot(all(colnames(actGeneExpr) == c(actC1subj, actC2subj)));
    stopifnot(all(colnames(actRnaSeq) == c(actC1subj, actC2subj)));
    stopifnot(all(colnames(actProtExpr) == c(actC1subj, actC2subj)));
    stopifnot(all(colnames(actProtExpr) == colnames(actRnaSeq)));
    
    # Filtrate genes with less than 30% per condition
    actRnaSeq <- actRnaSeq[
        rowSums(is.na(actRnaSeq[, actC1subj])) < .3*length(actC1subj) &
        rowSums(is.na(actRnaSeq[, actC2subj])) < .3*length(actC2subj)
    , ]
    actGeneExpr <- actGeneExpr[
        rowSums(is.na(actGeneExpr[, actC1subj])) < .3*length(actC1subj) &
        rowSums(is.na(actGeneExpr[, actC2subj])) < .3*length(actC2subj)
    , ]
    actProtExpr <- actProtExpr[
        rowSums(is.na(actProtExpr[, actC1subj])) < .3*length(actC1subj) &
        rowSums(is.na(actProtExpr[, actC2subj])) < .3*length(actC2subj)
    , ]
    
    # Filtrate genes with mean counts < 10
    actRnaSeq <- actRnaSeq[rowMeans(actRnaSeq[, actC1subj], na.rm=!F) > 10 &
        rowMeans(actRnaSeq[, actC2subj], na.rm=!F) > 10,];
    
    actRnaSeq <- DGEList(actRnaSeq, group=actSubtypes);
    actRnaSeq <- calcNormFactors(actRnaSeq, method=c('TMM'));
    actGeneExpr <- new('MAList', list(M=actGeneExpr));
    actProtExpr <- new('MAList', list(M=actProtExpr));
    fitOpts <- FitOptions(actSubtypes);
    
    actComp <- paste(actCtrst, collapse='_');
    
    actRnaSeq <- IGSAinput(name=paste('rnaSeq', actComp, sep='_'),
        expr_data=actRnaSeq, fit_options=fitOpts);
    actGeneExpr <- IGSAinput(name=paste('geneExpr', actComp, sep='_'),
        expr_data=actGeneExpr, fit_options=fitOpts);
    actProtExpr <- IGSAinput(name=paste('protExpr', actComp, sep='_'),
        expr_data=actProtExpr, fit_options=fitOpts);
    
    experiments <- list(actRnaSeq, actGeneExpr, actProtExpr);
    
    return(experiments);
})

# [1] "Basal" "LumA"  "25"    "12"   
# [1] "Basal" "Her2"  "25"    "25"   
# [1] "Basal" "LumB"  "25"    "20"   
# [1] "LumA" "Her2" "12"   "25"  
# [1] "LumA" "LumB" "12"   "20"  
# [1] "Her2" "LumB" "25"   "20"

# check how many contrasts for tcga
# xx <- unlist(lapply(allExperiments, function(a) {
#     unlist(lapply(a, function(x) {
#         aa <- x@fit_options;
#         paste(colnames(aa@design_matrix), collapse='_')
#     }))
# })); table(xx);
# xx <- unlist(lapply(allExperiments, function(a) {
#     unlist(lapply(a, function(x) {
#         aa <- x@fit_options;
#         paste(aa@contrast, collapse='_')
#     }))
# })); table(xx);

# igsaInput <- allExperiments[[1]][[1]];
gettreatVal <- function(igsaInput, dePerc=0.05, treatStep=0.05) {
    igsaInput <- getDEGenes(igsaInput);
    actDePerc <- length(MIGSA:::de_genes(seaParams(igsaInput))) / nrow(exprData(igsaInput));
    if (actDePerc < 0.04) {
        actSeaP <- SEAparams(adjust_method = 'none');
        seaParams(igsaInput) <- actSeaP;
        igsaInput <- getDEGenes(igsaInput);
        actDePerc <- length(MIGSA:::de_genes(seaParams(igsaInput))) / nrow(exprData(igsaInput));
        if (actDePerc < 0.04) {
            actSeaP <- SEAparams(adjust_method = 'none', de_cutoff=0.05);
            seaParams(igsaInput) <- actSeaP;
            igsaInput <- getDEGenes(igsaInput);
            actDePerc <- length(MIGSA:::de_genes(seaParams(igsaInput))) / nrow(exprData(igsaInput));
            if (actDePerc < 0.04) {
                actSeaP <- SEAparams(adjust_method = 'none', de_cutoff=0.1);
                seaParams(igsaInput) <- actSeaP;
                igsaInput <- getDEGenes(igsaInput);
                actDePerc <- length(MIGSA:::de_genes(seaParams(igsaInput))) / nrow(exprData(igsaInput));
            }
        }
    }
    oldDePerc <- 100;
    for (i in seq(0, 5, by=treatStep)) {
        actSeaP <- seaParams(igsaInput);
        seaParams(igsaInput) <- SEAparams(treat_lfc=i, adjust_method=actSeaP@adjust_method, de_cutoff=actSeaP@de_cutoff);
        igsaInput <- getDEGenes(igsaInput);
        
        actDePerc <- length(MIGSA:::de_genes(seaParams(igsaInput))) / nrow(exprData(igsaInput))-dePerc;
        print(c(i, actDePerc));
        if (actDePerc < 0) {
            treatRes <- i;
            if (abs(oldDePerc) < abs(actDePerc)) {
                treatRes <- i - treatStep;
                actDePerc <- oldDePerc;
                seaParams(igsaInput) <- SEAparams(treat_lfc=treatRes, adjust_method=actSeaP@adjust_method, de_cutoff=actSeaP@de_cutoff);
            }
            
            break;
        }
        oldDePerc <- actDePerc;
    }
    
    return(igsaInput);
}

allExperiments <- lapply(allExperiments, function(actExperiment) {
    lapply(actExperiment, gettreatVal);
})
allExperiments <- do.call(c, allExperiments);

aa <- do.call(rbind, lapply(allExperiments, summary));
# there are some datasets with DE genes percentages far from 5. For these, let's 
# try changing treat step for calculation

allExperimentsOK <- allExperiments[ as.numeric(4 < aa[, '%de_genes']) & 
    as.numeric(aa[, '%de_genes']) < 6 ];
allExperimentsBad <- allExperiments[ !(as.numeric(4 < aa[, '%de_genes']) & 
    as.numeric(aa[, '%de_genes']) < 6) ];

stopifnot(length(allExperiments) == length(allExperimentsOK) + length(allExperimentsBad));

allExperimentsBad <- lapply(allExperimentsBad, function(actExperiment) {
    gettreatVal(actExperiment, treatStep=0.01);
})

allExperimentsOK <- c(allExperimentsOK, allExperimentsBad);
stopifnot(length(allExperiments) == length(allExperimentsOK))

# Now every experiment has its parameters to get between 4 and 6 percent of their genes DE

## Supplementary: table SubjectsParams (only TCGA subjects)
allExperiments <- allExperimentsOK;
aa <- do.call(rbind, lapply(allExperiments, summary)); aa[order(aa[,1]),];

rm(list=setdiff(ls(), c('allExperiments')));

# load the gene sets to analyze
gSets <- list(
    BP=loadGo('BP'),
    CC=loadGo('CC'),
    MF=loadGo('MF'));

unlist(lapply(gSets, length))
#    BP    CC    MF 
# 15796  1902  4604
# org.Hs.eg.db_3.5.0

sum(unlist(lapply(gSets, length)))
# [1] 22302

migsaTCGAinput <- list(allExperiments=allExperiments, gSets=gSets);
# save(migsaTCGAinput, file='migsaTCGAinput.RData');

#################### And run MIGSA!

# set.seed(8818);
# load('migsaTCGAinput.RData');

library('MIGSA');
allExperiments <- migsaTCGAinput$allExperiments;
gSets <- migsaTCGAinput$gSets;
rm(migsaTCGAinput);

library('BiocParallel');
register(MulticoreParam(workers=2, threshold='DEBUG', progressbar=TRUE));

# In this step be careful not to give too many cores, as it consumes RAM
migsaTCGAres <- MIGSA(allExperiments, geneSets=gSets);
# saveRDS(migsaTCGAres, file='TCGAres.xz', compress='xz'); # File given in Suppl Data
