# alias R=R3.5.0;
# R
# stopifnot(R.Version()$version.string ==
#     'R Under development (unstable) (2017-12-13 r73907)');
setwd('~/MIGSAdata/'); # or any directory where Supplementary Data was downloaded to

library('GO.db');
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


# Download HKres.xz from
# https://figshare.com/articles/HKres_xz/6741689
migsaHKres <- readRDS('HKres.xz');

# set enrichment cutoff of 0.01
migsaHKres <- setEnrCutoff(migsaHKres, 0.01);

dim(migsaHKres);
# [1] 22302   121
migsaHKres <- migsaHKres[ rowSums(migsaHKres[,-(1:3)], na.rm=!F) > 0, ];
dim(migsaHKres); # enriched in at least one experiment
# [1] 9488  121

# get each experiment name, and contrasted subtypes
aux <- colnames(migsaHKres)[-(1:3)];
aux <- strsplit(gsub('LumA', 'A', gsub('LumB', 'B', gsub('Basal', 'Ba', gsub('Her2', 'H', aux)))), '_');
aux <- lapply(aux, function(x) unlist(strsplit(x, 'Vs')));

categories <- list();
categories$dataset <- factor(unlist(lapply(aux, function(x) x[[1]])));
categories$contrast <- factor(unlist(lapply(aux, function(x) paste(x[[2]], x[[3]], sep=''))));
categories$cond1 <- unlist(lapply(aux, function(x) x[[2]]));
categories$cond2 <- unlist(lapply(aux, function(x) x[[3]]));

# lines to get nice colors in heatmap
aux <- unique(c(categories$cond1, categories$cond2));
c1 <- setdiff(aux, unique(categories$cond1));
c2 <- setdiff(aux, unique(categories$cond2));
categories$cond1[which(categories$contrast==paste0(c2,c1))[[1]]] <- c1;
categories$cond2[which(categories$contrast==paste0(c2,c1))[[1]]] <- c2;

subtypes <- unique(unlist(categories[3:4]));
categories$cond1 <- factor(categories$cond1);
categories$cond2 <- factor(categories$cond2);

## Fig. 2
pdf('migsaHKresHeatmap.pdf');
aa <- migsaHeatmap(migsaHKres, categories=categories, categLabels=F)
dev.off();
## Supplementary File 1 EGSs
invisible(toCsv(aa$data, 'MigsaHk_EGS.csv'));

# lets get important enriching genes for each experiment
library('org.Hs.eg.db');
migsaHKresIEG <- lapply(colnames(aa$data), function(actCol) {
    print(actCol);
    actIdx <- grep(actCol, colnames(migsaHKres));
    actHKres <- migsaHKres[,c(1:3, actIdx)];
    genesinsets <- genesInSets(actHKres);
    genes <- lapply(1:nrow(genesinsets), function(i) {
        actGset <- genesinsets[i,];
        gsetGenes <- names(actGset)[actGset == 1];
        
        res <- list();
        if (length(gsetGenes) > 0) {
            res <- suppressMessages(select(org.Hs.eg.db, gsetGenes, 'SYMBOL', 'ENTREZID'));
            stopifnot(nrow(res) == length(gsetGenes));
            res <- res$SYMBOL;
        }
        return(res);
    })
    names(genes) <- rownames(genesinsets);
    return(genes);
})
names(migsaHKresIEG) <- colnames(aa$data);

stopifnot(all(lapply(migsaHKresIEG, length) == nrow(migsaHKres)));
# check all rownames are equal
stopifnot(all(apply(do.call(cbind, lapply(migsaHKresIEG, names)), 1, function(x) length(unique(x)))==1));

migsaHKresIEG <- do.call(cbind, lapply(migsaHKresIEG, function(actHKGenes) {
    unlist(lapply(actHKGenes, paste, collapse=' '));
}))

# genesInSets function returns as rownames GS_Name ++ '_' ++ id . So be carefull with this replacement
rownames(migsaHKresIEG) <- sub('.*?_', '', sub('.*?_', '', rownames(migsaHKresIEG)));
migsaHKresIEG <- migsaHKresIEG[rownames(aa$data),];
stopifnot(all(rownames(aa$data) == rownames(migsaHKresIEG)));
## Supplementary File 1 IEGs
invisible(toCsv(migsaHKresIEG, 'MigsaHk_IEG.csv'));


# heatmap per each subtype. if > 50% of enriched datasets
# library('VennDiagram');

stypestable <- data.frame(short=subtypes,
                        long=c('Basal-like', 'Luminal B', 'Her2-enriched', 'Luminal A'),
                        color=c('#00BFC4', '#7CAE00', '#C77CFF', '#F8766D'));
rownames(stypestable) <- stypestable$short;
percEnr <- 0; # 0 to 100
consPerc <- 0.5; # 0 to 1
## Supplementary Figs. 1-4
pdf('migsaHKresHeatmap_perSubtype.pdf');
stypeGSetsFull <- lapply(subtypes, function(actStype) {
#     actStype <- 'Ba';
    actCols <- which(categories[[3]] == actStype | categories[[4]] == actStype);
    actMigsaRes <- migsaHKres[, c(1:3, 3+actCols)];
    actCategories <- lapply(categories, function(x) x[actCols])[-2];
    actConds <- lapply(actCategories[2:3], function(x) {
        x <- as.character(x);
        x[x==actStype] <- '';
        return(x);
    })
    actCategories <- list(actCategories[[1]], 
        factor(paste0(actConds[[1]], actConds[[2]]), levels=levels(actCategories$cond2)));
    
    enrFilter <- percEnr*(ncol(actMigsaRes)-3)/100;
    
    migsaHeatmap(actMigsaRes, categories=actCategories, categLabels=F, enrFilter=enrFilter);
    
    actMigsaRes <- actMigsaRes[rowSums(actMigsaRes[,-(1:3)], na.rm=!F) >= enrFilter,];
    otherStypes <- setdiff(subtypes, actStype);
    enriched <- lapply(otherStypes, function(otherStype) {
#         otherStype <- 'A';
        otherCols <- actCategories[[2]] == otherStype;
        actMigsaRes[
            rowSums(actMigsaRes[, 3+which(otherCols)], na.rm=!F) > 
            consPerc*sum(otherCols), 'id'
        , drop=!F];
    });
    names(enriched) <- otherStypes;
    enriched <- enriched[order(names(enriched))]; # to have same colors as in heatmap
    
    enriched2Venn <- enriched;
    names(enriched2Venn) <- stypestable[names(enriched), 'long'];
    cols <- as.character(stypestable[names(enriched), 'color']);
    
#     tmp <- venn.diagram(enriched2Venn, filename=NULL, fill=cols, main=stypestable[actStype, 'long']);
#     grid.newpage();
#     grid.draw(tmp);
    
    return(enriched); # return total intersection
});
dev.off();
names(stypeGSetsFull) <- subtypes;

# lets put it stypeGSetsFull in matrix, by contrast
gsets <- unique(unlist(stypeGSetsFull));
aux <- rep(0, length(gsets)); names(aux) <- gsets;
contrasts <- combn(subtypes, 2);
CCET_EGS <- apply(contrasts, 2, function(actCtrst) {
    aux[unlist(stypeGSetsFull[[ actCtrst[1] ]][[ actCtrst[2] ]])] <- 1;
    return(aux);
})
colnames(CCET_EGS) <- apply(contrasts, 2, paste, collapse=' Vs ');
# lets save it as csv
## Supplementary File 2 EGSs
invisible(toCsv(CCET_EGS, 'CCET_EGS.csv'));

# lets analyze the genes that are mostly enriching these gene sets
library('org.Hs.eg.db');
specAllMigsaRes <- migsaHKres[migsaHKres$id %in% rownames(CCET_EGS),];
stopifnot(all(rownames(CCET_EGS) %in% specAllMigsaRes$id));
contrasts <- as.character(unique(categories$contrast));

CCET_IEG <- lapply(contrasts, function(actCtrst) {
#     actCtrst <- 'BaH';
    print(actCtrst);
    actCols <- which(actCtrst == as.character(categories$contrast));
    
    actMigsaRes <- specAllMigsaRes[, c(1:3, 3+actCols)];
    genesinsets <- genesInSets(actMigsaRes);
    genes <- lapply(1:nrow(genesinsets), function(i) {
        actGset <- sort(genesinsets[i,], decreasing=!F); # most important first
        gsetGenes <- names(actGset)[actGset > (ncol(actMigsaRes)-3)*consPerc];
        
        res <- list();
        if (length(gsetGenes) > 0) {
            res <- suppressMessages(select(org.Hs.eg.db, gsetGenes, 'SYMBOL', 'ENTREZID'));
            stopifnot(nrow(res) == length(gsetGenes));
            res <- res$SYMBOL;
        }
        return(res);
    })
    names(genes) <- rownames(genesinsets);
    return(genes);
})
names(CCET_IEG) <- contrasts;

stopifnot(all(lapply(CCET_IEG, length) == nrow(specAllMigsaRes)));
# all rownames are equal
stopifnot(all(apply(do.call(cbind, lapply(CCET_IEG, names)), 1, function(x) length(unique(x)))==1));

CCET_IEG <- do.call(cbind, lapply(CCET_IEG, function(actCtrstGenes) {
    unlist(lapply(actCtrstGenes, paste, collapse=' '));
}))

# genesInSets function returns as rownames GS_Name ++ '_' ++ id . So be carefull with this replacement
rownames(CCET_IEG) <- sub('.*?_', '', sub('.*?_', '', rownames(CCET_IEG)));
CCET_IEG <- CCET_IEG[rownames(CCET_EGS),];
stopifnot(all(rownames(CCET_EGS) == rownames(CCET_IEG)));
## Supplementary File 2 IEGs
invisible(toCsv(CCET_IEG, 'CCET_IEG.csv'));

# File also given in Suppl Data
# saveRDS(list(EGs=stypeGSetsFull, IEGs=CCET_IEG), file='CCETs.xz', compress='xz');
CCETs <- readRDS('CCETs.xz');
stypeGSetsFull <- CCETs$EGs;
CCET_IEG <- CCETs$IEGs;

stypeGSets <- lapply(stypeGSetsFull, function(enriched) Reduce(intersect, enriched));
names(stypeGSets) <- subtypes;

# cat(paste(lapply(names(stypeGSets), function(actStype) {
#     stype <- sub('\\<Ba\\>','Basal-like', sub('A', 'Luminal A', sub('H', 'Her2-enriched', sub('\\<B\\>', 'Luminal B', actStype))));
#     paste(stype, unlist(lapply(stypeGSets[[actStype]], Term)), sep=' | ', collapse='\n');
# }), collapse='\n'));

# lets put it in matrix
gsets <- unique(unlist(stypeGSets));
aux <- rep(0, length(gsets)); names(aux) <- gsets;
FTP_EGS <- do.call(cbind, lapply(subtypes, function(actStype) {
    aux[unlist(stypeGSets[actStype])] <- 1;
    return(aux);
}))
colnames(FTP_EGS) <- subtypes;

## Table 1
invisible(lapply(subtypes, function(actStype) {
#     actStype <- 'A';
    otherStypes <- setdiff(subtypes, actStype);
    res <- do.call(c, lapply(otherStypes, function(otherStype) {
#         otherStype <- otherStypes[[1]];
        length(stypeGSetsFull[[actStype]][[otherStype]]);
    }));
    
    names(res) <- otherStypes;
    res <- c(res, sum(FTP_EGS[, actStype]));
    print(res);
}));
#   B   H   A     
# 488 437 768 155 
#  Ba   H   A     
# 488 184 983  44 
#  Ba   B   A     
# 437 184 707  44 
#  Ba   B   H     
# 768 983 707 437


# pairwise intersections
invisible(lapply(subtypes, function(actStype) {
#     actStype <- 'A';
    otherStypes <- setdiff(subtypes, actStype);
    res <- do.call(c, lapply(otherStypes, function(otherStype) {
#         otherStype <- otherStypes[[1]];
        sum(rowSums(FTP_EGS[,c(actStype, otherStype)], na.rm=!F) == 2 & 
            rowSums(FTP_EGS, na.rm=!F) == 2
        );
    }));
    
    names(res) <- otherStypes;
    print(res);
}));
#  B  H  A 
#  8 10  7 
# Ba  H  A 
#  8  0  8 
# Ba  B  A 
# 10  0  0 
# Ba  B  H 
#  7  8  0

table(rowSums(FTP_EGS)); # exclusive and consensus gene sets
#   1   2   4 
# 506  33  27
colSums(FTP_EGS[rowSums(FTP_EGS) == 1,]);
#  Ba   B   H   A 
# 103   1   7 395

# lets save it as csv
## Supplementary File 3 EGSs
invisible(toCsv(FTP_EGS, 'FTP_EGS.csv'));

# lets analyze the genes that are mostly enriching these gene sets
FTP_IEG <- do.call(cbind, lapply(subtypes, function(actStype) {
    otherStypes <- setdiff(subtypes, actStype);
    contrasts <- c(paste0(otherStypes, actStype), paste0(actStype, otherStypes));
    contrasts <- contrasts[contrasts %in% colnames(CCET_IEG)];
    actIEG <- CCET_IEG[,contrasts];
    stopifnot(ncol(actIEG) == 3);
    
    actGenes <- unlist(lapply(seq_len(nrow(actIEG)), function(i) {
        stypeGenes <- strsplit(actIEG[i,], ' ');
        paste(intersect(intersect(
                stypeGenes[[1]],
                stypeGenes[[2]]),
                stypeGenes[[3]]),
            collapse=' ');
    }))
    
    names(actGenes) <- rownames(actIEG);
    return(actGenes);
}))
colnames(FTP_IEG) <- subtypes;

FTP_IEG <- FTP_IEG[rownames(FTP_EGS),];
stopifnot(all(rownames(FTP_EGS) == rownames(FTP_IEG)));
## Supplementary File 3 IEGs
invisible(toCsv(FTP_IEG, 'FTP_IEG.csv'));

# plot the three big GO tree with all subtypes
treeInfo <- matrix(!F, ncol=3, nrow=nrow(FTP_EGS));
rownames(treeInfo) <- rownames(FTP_EGS);
colnames(treeInfo) <- c('Enriched', 'Important', 'Color');
treeInfo[,'Important'] <- F;
treeInfo <- data.frame(treeInfo);

vec <- rep('', nrow(treeInfo));
aux <- do.call(rbind, lapply(colnames(FTP_EGS), function(actStype) {
    vec[FTP_EGS[,colnames(FTP_EGS) == actStype] == 1] <- actStype;
    return(vec);
}))
aux2 <- apply(aux, 2, paste, collapse='x'); table(aux2)
# BaxBxHxA   BaxBxx   BaxxHx    Baxxx   BaxxxA     xBxx    xBxxA     xxHx     xxxA 
#       27        8       10      103        7        1        8        7      395

aux2['xxxA' == aux2] <- '#f8766d';
aux2['Baxxx' == aux2] <- '#00bfc4';
aux2['xBxx' == aux2] <- '#7cae00';
aux2['xxHx' == aux2] <- '#c77cff';

aux2['xBxxA' == aux2] <- '#7f0000';
aux2['BaxxxA' == aux2] <- '#ff0000';
aux2['BaxBxx' == aux2] <- '#ffff00';
aux2['BaxxHx' == aux2] <- '#619cff';

aux2['BaxBxHxA' == aux2] <- 'grey75';

stopifnot(all(rowSums(FTP_EGS) != 3));

treeInfo$Color <- aux2;

legend <- matrix(
    c('Luminal A', 'Basal-like', 'Luminal B', 'Her2-enriched',
    'Lum A & Lum B', 'Lum A & Basal', 'Basal & Lum B', 'Basal & Her2-enriched',
    'In all of them',
    '#f8766d', '#00bfc4', '#7cae00', '#c77cff',
    '#7f0000', '#ff0000', '#ffff00', '#619cff',
    'grey75'),
    nrow=2, byrow=!F);

## Supplementary Figs. 5-7
pdf('FTP_trees.pdf', height=14, width=14);
invisible(lapply(c('BP', 'CC', 'MF'), function(ont) {
    plotInfo <- treeInfo[Ontology(rownames(treeInfo)) == ont, ];
    goTree(plotInfo, ont, legend, legendPos='bottomright');
}));
dev.off();
