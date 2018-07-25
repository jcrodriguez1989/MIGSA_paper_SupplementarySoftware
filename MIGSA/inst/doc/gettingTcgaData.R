### R code from vignette source 'gettingTcgaData.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: basalSubjects
###################################################
library(MIGSAdata);
data(tcgaMAdata);

names(tcgaMAdata$subtypes)[ tcgaMAdata$subtypes == "Basal" ];


###################################################
### code chunk number 2: lumaSubjects
###################################################
library(MIGSAdata);
data(tcgaMAdata);

names(tcgaMAdata$subtypes)[ tcgaMAdata$subtypes == "LumA" ];


###################################################
### code chunk number 3: tcgabiolinksCode (eval = FALSE)
###################################################
## ## Not run:
## 
## library(TCGAbiolinks);
## 
## R.Version()$version.string;
## # [1] "R version 3.2.3 (2015-12-10)"
## packageVersion("TCGAbiolinks");
## # [1] ‘1.0.10’
## 
## query <- TCGAquery(tumor="BRCA");
## matSamples <- TCGAquery_integrate(query);
## 
## # subjects in both microarray and RNAseq data
## matSamples["AgilentG4502A_07_3", "IlluminaHiSeq_RNASeq"];
## # [1] 495
## 
## # we filter only microarray data
## geneExprSubjects <- TCGAquery(tumor="BRCA", platform="AgilentG4502A_07_3", 
##     level=3);
## 
## # we filter only RNAseq data
## rnaSeqSubjects <- TCGAquery(tumor="BRCA", platform="IlluminaHiSeq_RNASeq", 
##     level=3);
## 
## geneExprbarcodes <- geneExprSubjects$barcode;
## geneExprbarcodes <- strsplit(geneExprbarcodes, ",");
## geneExprbarcodes <- Reduce(union, geneExprbarcodes);
## 
## rnaSeqbarcodes <- rnaSeqSubjects$barcode;
## rnaSeqbarcodes <- strsplit(rnaSeqbarcodes, ",");
## rnaSeqbarcodes <- Reduce(union, rnaSeqbarcodes);
## 
## commonSubjects <- intersect(geneExprbarcodes, rnaSeqbarcodes);
## rm(geneExprbarcodes); rm(rnaSeqbarcodes);
## length(commonSubjects);
## # [1] 547
## 
## # we filter microarray and RNAseq data (but just common subjects)
## geneExprSubjects <- TCGAquery(tumor="BRCA", platform="AgilentG4502A_07_3", 
##     samples=commonSubjects, level=3);
## rnaSeqSubjects <- TCGAquery(tumor="BRCA", platform="IlluminaHiSeq_RNASeq", 
##     samples=commonSubjects, level=3);
## 
## 
## #### this lines are the ones which are not working any more (TCGAdownload)
## # TCGAdownload(geneExprSubjects, path="geneExpr/", samples=commonSubjects);
## # TCGAdownload(rnaSeqSubjects,   path="rnaSeq/",   samples=commonSubjects, 
## #     type="gene.quantification");
## 
## ## However, we can provide you necessary files to skip the TCGAdownload step.
## 
## ## type is any of:
## # RNASeq:             exon.quantification
## #                     spljxn.quantification
## #                     gene.quantification
## # genome_wide_snp_6:  hg18.seg
## #                     hg19.seg,nocnv_hg18.seg
## #                     nocnv_hg19.seg
## 
## geneExpr <- TCGAprepare(geneExprSubjects, dir="geneExpr/");
## rnaSeq <- TCGAprepare(rnaSeqSubjects, dir="rnaSeq/", 
##     type="gene.quantification");
## 
## library(SummarizedExperiment);
## 
## assays(geneExpr);
## # names(1): raw_counts
## 
## # It would be a better way of conversion
## geneExpr <- head(assay(geneExpr, "raw_counts"), n=nrow(geneExpr));
## 
## assays(rnaSeq);
## # names(3): raw_counts median_length_normalized RPKM
## rnaSeq_raw <- head(assay(rnaSeq, "raw_counts"), n=nrow(rnaSeq));
## rnaSeq_medianNorm <- head(assay(rnaSeq, "median_length_normalized"), 
##     n=nrow(rnaSeq));
## rnaSeq_rpkm <- head(assay(rnaSeq, "RPKM"), n=nrow(rnaSeq));
## 
## ## checking if we have the same subjects in every experiment
## stopifnot(all(colnames(geneExpr) %in% colnames(rnaSeq_raw)));
## stopifnot(all(colnames(rnaSeq_raw) %in% colnames(rnaSeq_medianNorm)));
## stopifnot(all(colnames(rnaSeq_medianNorm) %in% colnames(rnaSeq_rpkm)));
## stopifnot(all(colnames(rnaSeq_rpkm) %in% colnames(geneExpr)));
## 
## mapping <- do.call(rbind, strsplit(rownames(rnaSeq_raw), "|", fixed=!F));
## colnames(mapping) <- c("Symbol", "Entrez");
## 
## #### Now let's get subjects subtypes
## 
## library(genefu);
## 
## rnaSeq <- rnaSeq_rpkm;
## rm(rnaSeq_rpkm);
## 
## ## Also request this file!
## pam50Annot <- read.csv("pam50_annotation.txt",sep="\t");
## 
## library(limma);
## dim(geneExpr);
## geneExpr <- avereps(geneExpr);
## dim(geneExpr);
## 
## rownames(rnaSeq) <- mapping[, "Symbol" ];
## dim(rnaSeq);
## rnaSeq <- rnaSeq[ mapping[, "Symbol" ] != "?" , ];
## dim(rnaSeq);
## rnaSeq <- avereps(rnaSeq);
## dim(rnaSeq);
## 
## geneExpr <- geneExpr[as.character(pam50Annot$GeneName),, drop=F];
## dim(geneExpr);
## rnaSeq <- rnaSeq[as.character(pam50Annot$GeneName),, drop=F];
## dim(rnaSeq);
## rnaSeq <- log(rnaSeq);
## 
## pam50Annot <- pam50Annot[,c("GeneName", "EntrezGene")];
## colnames(pam50Annot) <- c("probe", "EntrezGene.ID");
## pam50Annot$probe <- as.character(pam50Annot$probe);
## 
## ## get subtypes
## dataset <- apply(geneExpr, 1, as.numeric);
## rownames(dataset) <- colnames(geneExpr);
## subtypesGeneExpr <- intrinsic.cluster.predict(sbt.model=pam50.scale, 
##     data=dataset, annot=pam50Annot, do.mapping=!F, do.prediction.strength=!F, 
##     verbose=!F);
## 
## ## get subtypes
## dataset <- apply(rnaSeq, 1, as.numeric);
## rownames(dataset) <- colnames(rnaSeq);
## subtypesRnaSeq <- intrinsic.cluster.predict(sbt.model=pam50.scale, 
##     data=dataset, annot=pam50Annot, do.mapping=!F, do.prediction.strength=!F, 
##     verbose=!F);
## 
## table(subtypesGeneExpr$subtype);
## #  Basal   Her2   LumA   LumB Normal 
## #    101     77    150    157     62
## 
## table(subtypesRnaSeq$subtype);
## #  Basal   Her2   LumA   LumB Normal
## #    101     81    165    137     63
## 
## subtypesGeneExpr <- subtypesGeneExpr$subtype;
## subtypesRnaSeq <- subtypesRnaSeq$subtype[names(subtypesGeneExpr)];
## 
## ## how many subjects got the same subtype between microarray and RNAseq data
## concSubtypes <- table(subtypesGeneExpr, subtypesRnaSeq);
## concSubtypes;
## #          Basal Her2 LumA LumB Normal
## #   Basal     95    2    1    2      1
## #   Her2       0   72    0    4      1
## #   LumA       1    0  142    4      3
## #   LumB       3    7   19  127      1
## #   Normal     2    0    3    0     57
## sum(diag(concSubtypes)) / sum(concSubtypes);
## # [1] 0.9012797 # 90% of concordant subjects
## 
## stopifnot(all(names(subtypesGeneExpr) == names(subtypesRnaSeq)));
## 
## ## I am going to use the subjects that got the same classification in both
## subtypes <- subtypesGeneExpr[subtypesGeneExpr == subtypesRnaSeq];
## length(subtypes);
## # [1] 493
## 
## #### Now just translate GeneSymbols to EntrezGene IDs
## 
## ## Also request this file!
## annotAgi <- read.csv("AgilentG4502A_07_3.csv", sep="|");
## geneExprSymbol <- rownames(geneExpr);
## # we first search into Agilent annotation file
## geneExprEntrez <- annotAgi[ match(geneExprSymbol, annotAgi[, "Symbol"]), 
##     "Entrez" ];
## sum(is.na(geneExprEntrez));
## # [1] 796
## # then we look into the mapping given by RNASeq TCGA data
## geneExprEntrez[ is.na(geneExprEntrez) ] <- mapping[ match(geneExprSymbol[ 
##     is.na(geneExprEntrez) ], mapping[, "Symbol"]), "Entrez" ];
## sum(is.na(geneExprEntrez));
## # [1] 772
## 
## geneExpr <- geneExpr[ !is.na(geneExprEntrez), ];
## rownames(geneExpr) <- geneExprEntrez[ !is.na(geneExprEntrez) ];
## dim(geneExpr);
## geneExpr <- avereps(geneExpr);
## dim(geneExpr);
## 
## rownames(rnaSeq) <- do.call(rbind, strsplit(rownames(rnaSeq), "|", 
##     fixed=!F))[,2];
## dim(rnaSeq);
## rnaSeq <- avereps(rnaSeq);
## dim(rnaSeq);
## 
## load("rnaSeq_raw.RData");
## rownames(rnaSeq_raw) <- do.call(rbind, strsplit(rownames(rnaSeq_raw), "|", 
##     fixed=!F))[,2];
## dim(rnaSeq_raw);
## rnaSeq_raw <- avereps(rnaSeq_raw);
## dim(rnaSeq_raw);
## 
## #### And keep only Basal and Luminal A subjects
## rnaSeq_raw <- rnaSeq_raw[, names(subtypes)[subtypes %in% c("Basal", "LumA")] ];
## geneExpr <- geneExpr[, names(subtypes)[subtypes %in% c("Basal", "LumA")] ];
## 
## subtypes <- subtypes[subtypes %in% c("Basal", "LumA")];
## 
## ## And these are the two data objects used.
## tcgaRNAseqData <- list(rnaSeq=rnaSeq_raw, subtypes=subtypes);
## tcgaMAdata <- list(geneExpr=geneExpr, subtypes=subtypes);
## ## End(Not run)


###################################################
### code chunk number 4: Session Info
###################################################
sessionInfo()


