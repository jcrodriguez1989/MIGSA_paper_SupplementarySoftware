# alias R=R3.5.0;
# R
# stopifnot(R.Version()$version.string ==
#     'R Under development (unstable) (2017-12-13 r73907)');
setwd('~/NatureMIGSA/'); # or any directory where Supplementary Data was downloaded to

library('MIGSA');
library('AnnotationDbi');

# dataMatrix must be a matrix with GO ids as rownames (can have not GO).
toCsv <- function(dataMatrix, outFile) {
    res <- data.frame(id=rownames(dataMatrix));
    
    # with Ontology(rownames(dataMatrix)) or Term sometimes it crashed
    res$ontology <- NA;
    res$ontology[grepl('GO:', rownames(dataMatrix))] <- Ontology(rownames(dataMatrix)[grepl('GO:', rownames(dataMatrix))]);
    
    res$name <- NA;
    res$name[grepl('GO:', rownames(dataMatrix))] <- Term(rownames(dataMatrix)[grepl('GO:', rownames(dataMatrix))]);
    
#     res$ontology <- lapply(rownames(dataMatrix), Ontology);
#     res$name <- lapply(rownames(dataMatrix), Term);
    res$depth <- getHeights(rownames(dataMatrix));
    res <- cbind(res, data.frame(dataMatrix));
    
    write.table(res, file=outFile, sep='\t', row.names=F, col.names=!F, quote=F);
    invisible(return(res));
}


# Download TCGAres.xz from
# https://figshare.com/articles/TCGAres_xz/6741695
migsaTCGAres <- readRDS('TCGAres.xz');

# Download CCETs.xz from
# https://figshare.com/articles/CCETs_xz/6801263
CCETs <- readRDS('CCETs.xz');

## consensus heatmap against tcga results
stypeGSetsFull <- CCETs$EGs;
subtypes <- names(stypeGSetsFull);
allGsets <- unique(Reduce(union, lapply(stypeGSetsFull, function(aa) Reduce(union, aa))));
CCET_EGS <- matrix(F, nrow=length(allGsets), ncol=ncol(combn(subtypes,2)));
rownames(CCET_EGS) <- allGsets;
colnames(CCET_EGS) <- apply(combn(subtypes,2),2,paste0, collapse='');

invisible(lapply(names(stypeGSetsFull), function(actStype) {
    actCtrst <- stypeGSetsFull[[actStype]];
    lapply(names(actCtrst), function(actGset) {
        actGsets <- actCtrst[[actGset]];
        if (paste0(actStype, actGset)%in% colnames(CCET_EGS)) {
            CCET_EGS[actGsets, paste0(actStype, actGset)] <<- !F;
        } else {
            CCET_EGS[actGsets, paste0(actGset, actStype)] <<- !F;
        }
    });
}))

apply(CCET_EGS,2,sum);
# BaB BaH BaA  BH  BA  HA 
# 488 437 768 184 983 707

# set enrichment cutoff of 0.01
migsaTCGAres <- setEnrCutoff(migsaTCGAres, 0.01);

dim(migsaTCGAres);
# [1] 22302    21
migsaTCGAres <- migsaTCGAres[ rowSums(migsaTCGAres[,-(1:3)], na.rm=!F) > 0, ];
dim(migsaTCGAres); # enriched in at least one experiment
# [1] 4694   21

TCGA_EGS <- as.data.frame(migsaTCGAres);
TCGA_EGS <- TCGA_EGS[,c(1, 4:ncol(TCGA_EGS))]; # discard Name and GS_Name
dim(TCGA_EGS);
# [1] 4694   19

colnames(TCGA_EGS);
#  [1] "id"                  "geneExpr_Basal_Her2" "geneExpr_Basal_LumA"
#  [4] "geneExpr_Basal_LumB" "geneExpr_Her2_LumB"  "geneExpr_LumA_Her2" 
#  [7] "geneExpr_LumA_LumB"  "protExpr_Basal_Her2" "protExpr_Basal_LumA"
# [10] "protExpr_Basal_LumB" "protExpr_Her2_LumB"  "protExpr_LumA_Her2" 
# [13] "protExpr_LumA_LumB"  "rnaSeq_Basal_Her2"   "rnaSeq_Basal_LumA"  
# [16] "rnaSeq_Basal_LumB"   "rnaSeq_Her2_LumB"    "rnaSeq_LumA_Her2"   
# [19] "rnaSeq_LumA_LumB"
unique(sub('.*?_', '', colnames(TCGA_EGS)));
# [1] "id"         "Basal_Her2" "Basal_LumA" "Basal_LumB" "Her2_LumB" 
# [6] "LumA_Her2"  "LumA_LumB"

CCET_EGS <- data.frame(CCET_EGS);
colnames(CCET_EGS);
# [1] "BaB" "BaH" "BaA" "BH"  "BA"  "HA" 
colnames(CCET_EGS) <- sub('BA', 'AB', sub('HA', 'AH', sub('BH', 'HB', colnames(CCET_EGS))));
colnames(CCET_EGS) <- paste0('C_', colnames(CCET_EGS));
CCET_EGS$id <- rownames(CCET_EGS);

# shorten technology
colnames(TCGA_EGS) <- sub('geneExpr', 'M', colnames(TCGA_EGS));
colnames(TCGA_EGS) <- sub('protExpr', 'I', colnames(TCGA_EGS));
colnames(TCGA_EGS) <- sub('rnaSeq', 'R', colnames(TCGA_EGS));

# shorten contrast (and put it as in consensus data.frame)
colnames(TCGA_EGS) <- sub('Basal_Her2', 'BaH', colnames(TCGA_EGS));
colnames(TCGA_EGS) <- sub('Basal_LumA', 'BaA', colnames(TCGA_EGS));
colnames(TCGA_EGS) <- sub('Basal_LumB', 'BaB', colnames(TCGA_EGS));
colnames(TCGA_EGS) <- sub('Her2_LumB', 'HB', colnames(TCGA_EGS));
colnames(TCGA_EGS) <- sub('LumA_Her2', 'AH', colnames(TCGA_EGS));
colnames(TCGA_EGS) <- sub('LumA_LumB', 'AB', colnames(TCGA_EGS));

dim(CCET_EGS); dim(TCGA_EGS);
# [1] 1728    7
# [1] 4694   19
CCET_Wtcga_EGS <- merge(CCET_EGS, TCGA_EGS, all=!F, by='id');
dim(CCET_Wtcga_EGS);
# [1] 4885   25
gs_ids <- CCET_Wtcga_EGS[,1];
CCET_Wtcga_EGS <- CCET_Wtcga_EGS[,-1];
CCET_Wtcga_EGS <- apply(CCET_Wtcga_EGS,2,as.logical);
rownames(CCET_Wtcga_EGS) <- gs_ids;

# categories
istcga <- rep('C', ncol(CCET_Wtcga_EGS));
istcga[!grepl('C_', colnames(CCET_Wtcga_EGS))] <- 'T';

tech <- rep('M', ncol(CCET_Wtcga_EGS));
tech[grep('I_', colnames(CCET_Wtcga_EGS))] <- 'I';
tech[grep('R_', colnames(CCET_Wtcga_EGS))] <- 'R';

ctrst <- sub('.*?_', '', colnames(CCET_Wtcga_EGS));

cond1 <- unlist(lapply(ctrst, function(char) substr(char, nchar(char), nchar(char))));
cond2 <- unlist(lapply(ctrst, function(char) substr(char, 1, 1)));
cond2[grep('Ba', ctrst)] <- 'Ba';

stopifnot(all(paste0(cond2, cond1) == ctrst));

# lines to get nice colors in heatmap
aux <- unique(c(cond1, cond2));
c1 <- setdiff(aux, unique(cond1));
c2 <- setdiff(aux, unique(cond2));
cond1[which(ctrst==paste0(c1,c2))[[1]]] <- c1;
cond2[which(ctrst==paste0(c1,c2))[[1]]] <- c2;

wtcgaCateg <- list(istcga, tech, ctrst, cond2, cond1);

## Fig. 2
pdf('CCET_Wtcga.pdf');
aa <- migsaHeatmap(CCET_Wtcga_EGS, categories=wtcgaCateg, categLabels=F);
dev.off();

# save csv
## Supplementary File 6 EGS
aux <- aa$data[rownames(aa$data) %in% TCGA_EGS$id,
                colnames(aa$data) %in% colnames(TCGA_EGS)];
invisible(toCsv(aux, 'MigsaTCGA_EGS.csv'));
dim(aux); dim(TCGA_EGS); # CCET was removed for Suppl File
# [1] 4694   18
# [1] 4694   19

## get important genes for tcga enrichments
library('org.Hs.eg.db');
migsaTCGAresIEG <- lapply(4:ncol(migsaTCGAres), function(i) {
    print(colnames(migsaTCGAres)[[i]]);
    actTCGAres <- migsaTCGAres[,c(1:3,i)];
    genesinsets <- genesInSets(actTCGAres);
    genes <- lapply(1:nrow(genesinsets), function(j) {
        actGset <- genesinsets[j,];
        gsetGenes <- names(actGset)[actGset == 1];
        
        res <- list();
        if (length(gsetGenes) > 0) {
            res <- suppressMessages(select(org.Hs.eg.db, gsetGenes, "SYMBOL", "ENTREZID"));
            stopifnot(nrow(res) == length(gsetGenes));
            res <- res$SYMBOL;
        }
        return(res);
    })
    names(genes) <- rownames(genesinsets);
    return(genes);
})
names(migsaTCGAresIEG) <- colnames(TCGA_EGS)[-1];

stopifnot(all(lapply(migsaTCGAresIEG, length) == nrow(migsaTCGAres)));
# all rownames are equal
stopifnot(all(apply(do.call(cbind, lapply(migsaTCGAresIEG, names)), 1, function(x) length(unique(x)))==1));

migsaTCGAresIEG <- do.call(cbind, lapply(migsaTCGAresIEG, function(actTCGAGenes) {
    unlist(lapply(actTCGAGenes, paste, collapse=" "));
}))

# genesInSets function returns as rownames GS_Name ++ '_' ++ id . So be carefull with this replacement
rownames(migsaTCGAresIEG) <- sub('.*?_', '', sub('.*?_', '', rownames(migsaTCGAresIEG)));
migsaTCGAresIEG <- migsaTCGAresIEG[rownames(aux),colnames(aux)];
stopifnot(all(rownames(aux) == rownames(migsaTCGAresIEG)));
## Supplementary File 6 IEG
invisible(toCsv(migsaTCGAresIEG, 'MigsaTCGA_IEG.csv'));

# add CCET_IEG to these genes
CCET_IEG <- CCETs$IEG;
colnames(CCET_IEG) <- paste0('C_', colnames(CCET_IEG));
CCET_IEG <- cbind(id=rownames(CCET_IEG), CCET_IEG);
migsaTCGAresIEG <- cbind(id=rownames(migsaTCGAresIEG), migsaTCGAresIEG);

dim(CCET_IEG); dim(migsaTCGAresIEG);
# [1] 1728    7
# [1] 4694   19
CCET_Wtcga_IEG <- as.matrix(merge(CCET_IEG, migsaTCGAresIEG, all=!F, by='id'));
CCET_Wtcga_IEG[is.na(CCET_Wtcga_IEG)] <- '';
dim(CCET_Wtcga_IEG);
# [1] 4885   25

gs_ids <- CCET_Wtcga_IEG[,1];
CCET_Wtcga_IEG <- CCET_Wtcga_IEG[,-1];
rownames(CCET_Wtcga_IEG) <- gs_ids;


# lets intersect for each contrast Microarray and RNA-Seq as TET
contrasts <- unique(sub('.*?_', '', colnames(CCET_Wtcga_EGS)));
techs <- c('CCET', 'PET', 'TET');
migsaEGS <- do.call(cbind, lapply(contrasts, function(actCtrst) {
#     actCtrst <- 'BaB';
    actMatrix <- CCET_Wtcga_EGS[, sub('.*?_', '', colnames(CCET_Wtcga_EGS)) == actCtrst];
    newActMatrix <- actMatrix[, paste0(c('C_', 'I_'), actCtrst)];
    newActMatrix <- cbind(newActMatrix, 
        rowSums(actMatrix[, paste0(c('M_', 'R_'), actCtrst)], na.rm=!F) == 2);
    colnames(newActMatrix) <- paste(techs, actCtrst, sep='_');
    stopifnot(ncol(newActMatrix) == 3);
    return(newActMatrix);
}))

# now get the genes
migsaIEG <- do.call(cbind, lapply(contrasts, function(actCtrst) {
#     actCtrst <- 'BaB';
    actMatrix <- CCET_Wtcga_IEG[, sub('.*?_', '', colnames(CCET_Wtcga_IEG)) == actCtrst];
    tetMatrix <- actMatrix[, paste0(c('M_', 'R_'), actCtrst)];
    actGenes <- unlist(lapply(seq_len(nrow(tetMatrix)), function(i) {
        actGenes <- strsplit(as.character(tetMatrix[i,]), ' ');
        paste(intersect(
                actGenes[[1]],
                actGenes[[2]]),
            collapse=' ');
    }))
    
    newActMatrix <- cbind(
        actMatrix[, paste0(c('C_', 'I_'), actCtrst)], 
        actGenes);

    colnames(newActMatrix) <- paste(techs, actCtrst, sep='_');
    stopifnot(ncol(newActMatrix) == 3);
    return(newActMatrix);
}))


stopifnot(all(dim(migsaEGS) == dim(migsaIEG)));
stopifnot(all(colnames(migsaEGS) == colnames(migsaIEG)));
stopifnot(all(rownames(migsaEGS) %in% rownames(migsaIEG)));

# give the same order for IEG
migsaIEG <- migsaIEG[rownames(migsaEGS),];
stopifnot(all(rownames(migsaEGS) == rownames(migsaIEG)));

migsaIEG <- migsaIEG[rowSums(migsaEGS, na.rm=!F) > 0,];
migsaEGS <- migsaEGS[rowSums(migsaEGS, na.rm=!F) > 0,];
migsaEGS[is.na(migsaEGS)] <- FALSE;

## Supplementary File 7 EGS
invisible(toCsv(migsaEGS, 'Migsa_EGS.csv'));


## Supplementary File 7 IEG
invisible(toCsv(migsaIEG, 'Migsa_IEG.csv'));


# lets do a table for each subtype
FTP_Wtcga_EGS <- do.call(cbind, lapply(subtypes, function(actStype) {
#     actStype <- 'A'; otherStype <- 'H'; acttech <- 'CCET';
    res <- do.call(cbind, lapply(techs, function(acttech) {
        actCols <- unlist(lapply(setdiff(subtypes, actStype), function(otherStype) {
            actCtrst <- paste0(acttech, '_', actStype, otherStype);
            if (length(grep(actCtrst, colnames(migsaEGS))) == 0)
                actCtrst <- paste0(acttech, '_', otherStype, actStype);
            return(actCtrst);
        }))
        rowSums(migsaEGS[, actCols], na.rm=!F) == 3;
    }));
    colnames(res) <- paste(techs, actStype, sep='_');
    return(res);
}))

# now get the genes
FTP_Wtcga_IEG <- do.call(cbind, lapply(subtypes, function(actStype) {
#     actStype <- 'A'; otherStype <- 'H'; acttech <- 'CCET';
    res <- do.call(cbind, lapply(techs, function(acttech) {
        actCols <- unlist(lapply(setdiff(subtypes, actStype), function(otherStype) {
            actCtrst <- paste0(acttech, '_', actStype, otherStype);
            if (length(grep(actCtrst, colnames(migsaIEG))) == 0)
                actCtrst <- paste0(acttech, '_', otherStype, actStype);
            return(actCtrst);
        }))
        stopifnot(length(actCols) == 3);
        actIEG <- migsaIEG[, actCols];
        actGenes <- unlist(lapply(seq_len(nrow(actIEG)), function(i) {
            actGenes <- strsplit(as.character(actIEG[i,]), ' ');
            paste(intersect(intersect(
                    actGenes[[1]],
                    actGenes[[2]]),
                    actGenes[[3]]),
                collapse=' ');
        }))
    }));
    colnames(res) <- paste(techs, actStype, sep='_');
    return(res);
}))
rownames(FTP_Wtcga_IEG) <- rownames(migsaIEG);

stopifnot(all(dim(FTP_Wtcga_EGS) == dim(FTP_Wtcga_IEG)));
stopifnot(all(colnames(FTP_Wtcga_EGS) == colnames(FTP_Wtcga_IEG)));
stopifnot(all(rownames(FTP_Wtcga_EGS) == rownames(FTP_Wtcga_IEG)));

FTP_Wtcga_IEG <- FTP_Wtcga_IEG[rowSums(FTP_Wtcga_EGS, na.rm=!F) > 0,];
FTP_Wtcga_EGS <- FTP_Wtcga_EGS[rowSums(FTP_Wtcga_EGS, na.rm=!F) > 0,];

## Supplementary File 8 EGS
invisible(toCsv(FTP_Wtcga_EGS, 'FTP_Wtcga_EGS.csv'));


## Supplementary File 8 IEG
invisible(toCsv(FTP_Wtcga_IEG, 'FTP_Wtcga_IEG.csv'));


library('VennDiagram');
stypestable <- data.frame(short=subtypes,
                        long=c('Basal-like', 'Luminal B', 'Her2-enriched', 'Luminal A'),
                        color=c('#00BFC4', '#7CAE00', '#C77CFF', '#F8766D'));
rownames(stypestable) <- stypestable$short;

## Supplementary Figs. 12-15
pdf("FTP_Wtcga_EGS_venns.pdf");
byStype <- do.call(rbind, lapply(subtypes, function(actStype) {
#     actStype <- 'A'; otherStype <- 'H'; acttech <- 'CCET';
    actCols <- paste(techs, actStype, sep='_');
    actMatrix <- FTP_Wtcga_EGS[,actCols];
    actMatrix[is.na(actMatrix)] <- F;
    stopifnot(ncol(actMatrix) == 3);
    
    forVenn <- lapply(colnames(actMatrix), function(acttech) {
        rownames(actMatrix)[actMatrix[, acttech]];
    })
    names(forVenn) <- sub('CCET', 'FTP', sub('_.*', '', colnames(actMatrix)));
    tmp <- venn.diagram(forVenn, filename=NULL, fill=rainbow(3), main=stypestable[actStype, 'long']);
    grid.newpage();
    grid.draw(tmp);
    
    colnames(actMatrix) <- sub('_.*', '', colnames(actMatrix));
    
    c(
        sum(actMatrix[,'CCET'] == 1 & rowSums(actMatrix) == 1),
        sum(actMatrix[,'PET'] == 1 & rowSums(actMatrix) == 1),
        sum(actMatrix[,'TET'] == 1 & rowSums(actMatrix) == 1),
        
        sum(rowSums(actMatrix[,c('CCET','PET')]) == 2 & rowSums(actMatrix) == 2),
        sum(rowSums(actMatrix[,c('CCET','TET')]) == 2 & rowSums(actMatrix) == 2),
        sum(rowSums(actMatrix[,c('PET','TET')]) == 2 & rowSums(actMatrix) == 2),
        
        sum(rowSums(actMatrix) == 3),
        
        sum(rowSums(actMatrix) > 0));
}))
dev.off();

byStype <- rbind(byStype, colSums(byStype));
rownames(byStype) <- c(subtypes, 'total');
colnames(byStype) <- c('FTP', 'PET', 'TET', 'FTP&PET', 'FTP&TET', 'PET&TET', 'FTP&PET&TET', 'total');
## Table 2
byStype; # View(byStype);
#        FTP PET TET  FTP&PET  FTP&TET PET&TET  FTP&PET&TET total
# Ba     109  46  21       12       19       0           15   222
# B       37   7  25        0        7       0            0    76
# H       42  20   1        0        2       0            0    65
# A      234   2  15        3      186       0           14   454
# total  422  75  62       15      214       0           29   817

# and now let's plot some GO trees

colnames(FTP_Wtcga_EGS);
#  [1] "CCET_Ba" "PET_Ba"  "TET_Ba"  "CCET_B"  "PET_B"   "TET_B"   "CCET_H" 
#  [8] "PET_H"   "TET_H"   "CCET_A"  "PET_A"   "TET_A"

library('GO.db');
library('graph');
treeFun <- function(dataMatrix) {
    actStype <- unique(sub('.*?_', '', colnames(dataMatrix)));
    stopifnot(length(actStype) == 1);
    
    colnames(dataMatrix) <- sub('_.*', '', colnames(dataMatrix));
    
    dataMatrix <- dataMatrix[rowSums(dataMatrix, na.rm=!F) > 0,];
    
    treeInfo <- matrix(!F, ncol=3, nrow=nrow(dataMatrix));
    rownames(treeInfo) <- rownames(dataMatrix);
    colnames(treeInfo) <- c("Enriched", "Important", "Color");
    treeInfo[,"Important"] <- F;
    treeInfo <- data.frame(treeInfo);
    
    vec <- rep("", nrow(treeInfo));
    aux <- do.call(rbind, lapply(colnames(dataMatrix), function(actCol) {
        vec[dataMatrix[,colnames(dataMatrix) == actCol] == 1] <- actCol;
        return(vec);
    }))
    aux2 <- apply(aux, 2, paste, collapse="x");
    table(aux2);
    
    aux2["CCETxx" == aux2] <- "firebrick";
    aux2["xPETx" == aux2] <- "green";
    aux2["xxTET" == aux2] <- "dodgerblue2";
    
    aux2["CCETxPETx" == aux2] <- "yellow";
    aux2["CCETxxTET" == aux2] <- "darkorchid1";
#     aux2["xPETxTET" == aux2] <- "cyan";
    stopifnot(all(!"xPETxTET" == aux2));
    
    aux2["CCETxPETxTET" == aux2] <- "grey75";
    
    treeInfo$Color <- aux2;
    
    legend <- matrix(c('FTP', 'PET', 'TET', 'FTP & PET', 'FTP & TET', 'In all of them', 'firebrick', 'green', 'dodgerblue2', 'yellow', 'darkorchid1', 'grey75'), nrow=2, byrow=!F);
    
    res <- lapply(c("BP", "CC", "MF"), function(ont) {
        plotInfo <- treeInfo[Ontology(rownames(treeInfo)) == ont, ];
        goTree(plotInfo, ont, legend, legendPos="bottomright");
    });
    return(res);
}

## Supplementary Figs. 16-27
# actStype <- subtypes[[1]];
trees <- lapply(subtypes, function(actStype) {
    print(actStype);
    actData <- FTP_Wtcga_EGS[,grep(paste0('_', actStype, '$'), colnames(FTP_Wtcga_EGS))];
    
    # for BP remove some gene sets, if not the tree is too big.
    # from BP we are removing gene sets of depth less or equal to 3
#     actData <- actData[Ontology(rownames(actData)) != 'BP' | getHeights(rownames(actData)) > 3,];
    
    pdf(paste0('FTP_Wtcga_trees_', actStype, '.pdf'), height=14, width=14);
    actRes <- treeFun(actData);
    dev.off();
    return(actRes);
})
names(trees) <- subtypes;

## lets get for each subtype, its PET exclusive terms
petExclusives <- lapply(subtypes, function(actStype) {
    actData <- FTP_Wtcga_EGS[,grep(paste0('_', actStype, '$'), colnames(FTP_Wtcga_EGS))];
    
    rownames(actData)[rowSums(actData, na.rm=!F) == 1 & actData[, paste0('PET_', actStype)]];
}); names(petExclusives) <-  subtypes;

# cat(paste(lapply(names(petExclusives), function(actStype) {
#     stype <- sub('\\<Ba\\>','Basal-like', sub('A', 'Luminal A', sub('H', 'Her2-enriched', sub('\\<B\\>', 'Luminal B', actStype))));
#     paste(stype, unlist(lapply(petExclusives[[actStype]], Term)), sep=' | ', collapse='\n');
# }), collapse='\n'));


## lets get for each subtype, its FTPP exclusive terms
ftppExclusives <- lapply(subtypes, function(actStype) {
    actData <- FTP_Wtcga_EGS[,grep(paste0('_', actStype, '$'), colnames(FTP_Wtcga_EGS))];
    
    rownames(actData)[rowSums(actData, na.rm=!F) == 3];
}); names(ftppExclusives) <-  subtypes;

# cat(paste(lapply(names(ftppExclusives), function(actStype) {
#     stype <- sub('\\<Ba\\>','Basal-like', sub('A', 'Luminal A', sub('H', 'Her2-enriched', sub('\\<B\\>', 'Luminal B', actStype))));
#     paste(stype, unlist(lapply(ftppExclusives[[actStype]], Term)), sep=' | ', collapse='\n');
# }), collapse='\n'));


## lets get for each subtype, its PET exclusive genes
petExclusiveGenes <- lapply(subtypes, function(actStype) {
    actData <- FTP_Wtcga_IEG[,grep(paste0('_', actStype, '$'), colnames(FTP_Wtcga_IEG))];
    stopifnot(ncol(actData) == 3);
    actData <- actData[petExclusives[[actStype]],, drop=F]; # only analyze the pet exclusive gene sets
    actPetGenes <- actData[, grepl(paste0('PET_', actStype), colnames(actData)), drop=F];
    actOtherGenes <- actData[, !grepl(paste0('PET_', actStype), colnames(actData))];
    
    actGenes <- unlist(lapply(seq_len(nrow(actPetGenes)), function(i) {
        aux <- strsplit(as.character(actOtherGenes[i,]), ' ');
        setdiff(
            strsplit(as.character(actPetGenes[i,]), ' ')[[1]],
            union(aux[[1]], aux[[2]])
        )
    }))
    return(unique(actGenes));
}); names(petExclusiveGenes) <-  subtypes;

# cat(paste(lapply(names(petExclusiveGenes), function(actStype) {
#     stype <- sub('\\<Ba\\>','Basal-like', sub('A', 'Luminal A', sub('H', 'Her2-enriched', sub('\\<B\\>', 'Luminal B', actStype))));
#     paste(stype, petExclusiveGenes[[actStype]], sep=' | ', collapse='\n');
# }), collapse='\n'));


## lets get for each subtype, its FTPP exclusive genes
ftppExclusiveGenes <- lapply(subtypes, function(actStype) {
    actData <- FTP_Wtcga_IEG[,grep(paste0('_', actStype, '$'), colnames(FTP_Wtcga_IEG))];
    stopifnot(ncol(actData) == 3);
    
    actGenes <- unlist(lapply(seq_len(nrow(actData)), function(i) {
        stypeGenes <- strsplit(actData[i,], ' ');
        intersect(intersect(
                stypeGenes[[1]],
                stypeGenes[[2]]),
                stypeGenes[[3]])
    }))
    return(unique(actGenes));
}); names(ftppExclusiveGenes) <-  subtypes;

# cat(paste(lapply(names(ftppExclusiveGenes), function(actStype) {
#     stype <- sub('\\<Ba\\>','Basal-like', sub('A', 'Luminal A', sub('H', 'Her2-enriched', sub('\\<B\\>', 'Luminal B', actStype))));
#     paste(stype, ftppExclusiveGenes[[actStype]], sep=' | ', collapse='\n');
# }), collapse='\n'));
