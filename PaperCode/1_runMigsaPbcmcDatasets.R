# alias R=R3.5.0;
# R
# stopifnot(R.Version()$version.string ==
#     'R Under development (unstable) (2017-12-13 r73907)');
setwd('~/NatureMIGSA/'); # or any directory where Supplementary Data was downloaded to

library('MIGSA');
# set.seed(8818);

# Download HKdatasets.xz from
# https://figshare.com/articles/HKdatasets_xz/6741698
HaibeKains <- readRDS('HKdatasets.xz');

# get the subtypes according to pbcmc supplementary material
## Supplementary Table 1
HaibeKainsWPAM50 <- lapply(names(HaibeKains), function(actExp) {
#     actExp <- names(HaibeKains)[[1]];
    actDs <- HaibeKains[[actExp]];
    
    permSubtypes <- actDs$demo;
    # we are going to use Ambiguous, so we put them as Assigned
    permSubtypes$Permutation.Results[ permSubtypes$Permutation.Results == 'Ambiguous' ] <- 'Assigned';
    subtypes <- as.character(permSubtypes$Genefu);
    subtypes[ permSubtypes$Permutation.Results == 'Not Assigned' ] <- 'NotAssigned';
    names(subtypes) <- permSubtypes$Samplename;
    
    print(actExp);
    print(c(table(subtypes), length(subtypes)));
    
    stopifnot(all(actDs$demo$samplename == colnames(actDs$exprMatrix)));

    return(list(exprMatrix=actDs$exprMatrix, classes=subtypes));
})
names(HaibeKainsWPAM50) <- names(HaibeKains);

# get all possible contrasts, each subtype must have more than 8 subjects to be used
minSamples <- 8;
contrasts <- combn(c('Basal', 'LumA', 'Her2', 'LumB'), 2);
HaibeKainsWPAM50Igsas <- lapply(names(HaibeKainsWPAM50), function(actName) {
#     actName <- names(HaibeKainsWPAM50)[[1]];
    actExp <- HaibeKainsWPAM50[[ actName ]];
    actIgsas <- apply(contrasts, 2, function (actCtrst) {
#         actCtrst <- contrasts[,1];
        c1 <- actCtrst[[1]];
        c2 <- actCtrst[[2]];
        if (sum(actExp$classes == c1, na.rm=!F) >= minSamples && 
            sum(actExp$classes == c2, na.rm=!F) >= minSamples) {
            actExpr <- actExp$exprMatrix[,actExp$classes %in% actCtrst];
            actClasses <- actExp$classes[actExp$classes %in% actCtrst];
            
            actC1subj <- names(actClasses)[ actClasses == c1 ]; # condition 1 subjects
            actC2subj <- names(actClasses)[ actClasses == c2 ]; # condition 2 subjects
            
            # leave genes with at least 30% samples reads
            actExpr <- actExpr[
                        rowSums(is.na(actExpr[, actC1subj])) < .3*length(actC1subj) &
                        rowSums(is.na(actExpr[, actC2subj])) < .3*length(actC2subj)
                    , ]
            
            actExpr <- new('MAList',list(M=actExpr));
            fit_options <- FitOptions(as.character(actClasses));
            stopifnot(all(fit_options@contrast == c(-1,1)));
            actIgsa <- IGSAinput(
                name=paste(actName, paste(actCtrst, collapse='Vs'), sep='_'),
                expr_data=actExpr,
                fit_options=fit_options);
            return(actIgsa);
        } else {
            return(NA);
        }
    })
    actIgsas <- actIgsas[!is.na(actIgsas)];
    return(actIgsas);
})
HaibeKainsWPAM50Igsas <- do.call(c, unlist(HaibeKainsWPAM50Igsas, recursive=FALSE));
HaibeKainsWPAM50Igsas <- HaibeKainsWPAM50Igsas[!is.na(HaibeKainsWPAM50Igsas)];

# check how many contrasts for HK
# xx <- unlist(lapply(HaibeKainsWPAM50Igsas, function(x) { 
#     aa <- x@fit_options;
#     paste(colnames(aa@design_matrix), collapse='_')
# })); table(xx);
# xx <- unlist(lapply(HaibeKainsWPAM50Igsas, function(x) { 
#     aa <- x@fit_options;
#     paste(aa@contrast, collapse='_')
# })); table(xx);

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

# calculate the parameters that get 5% of the genes as differentially expressed
allExperiments <- lapply(HaibeKainsWPAM50Igsas, function(actExperiment) {
    gettreatVal(actExperiment);
})

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

aa <- do.call(rbind, lapply(allExperimentsBad, summary));
# there are a few more datasets far from %5. Repeat with lower step

allExperimentsOK <- c(allExperimentsOK, allExperimentsBad[ 
    as.numeric(4 < aa[, '%de_genes']) & as.numeric(aa[, '%de_genes']) < 6 ]);
allExperimentsBad <- allExperimentsBad[ !(as.numeric(4 < aa[, '%de_genes']) & 
    as.numeric(aa[, '%de_genes']) < 6) ];

stopifnot(length(allExperiments) == length(allExperimentsOK) + length(allExperimentsBad));

allExperimentsBad <- lapply(allExperimentsBad, function(actExperiment) {
    gettreatVal(actExperiment, treatStep=0.001);
})

do.call(rbind, lapply(allExperimentsBad, summary));

allExperimentsOK <- c(allExperimentsOK, allExperimentsBad);
stopifnot(length(allExperiments) == length(allExperimentsOK))

# Now every experiment has its parameters to get between 4 and 6 percent of their genes DE

## Supplementary Table 2 (only HaibeKains subjects)
allExperiments <- allExperimentsOK;
aa <- do.call(rbind, lapply(allExperiments, summary)); aa[order(aa[,1]),];

rm(list=setdiff(ls(), 'allExperiments'));

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

migsaHKinput <- list(allExperiments=allExperiments, gSets=gSets);
# save(migsaHKinput, file='migsaHKinput.RData');

#################### And run MIGSA!

# set.seed(8818);
# load('migsaHKinput.RData');

library('MIGSA');
allExperiments <- migsaHKinput$allExperiments;
gSets <- migsaHKinput$gSets;
rm(migsaHKinput);

library('BiocParallel');
register(MulticoreParam(workers=2, threshold='DEBUG', progressbar=TRUE));

# In this step be careful not to give too many cores, as it consumes RAM
migsaHKres <- MIGSA(allExperiments, geneSets=gSets);
# saveRDS(migsaHKres, file='HKres.xz', compress='xz'); # File given in Suppl Data
