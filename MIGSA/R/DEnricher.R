setGeneric(name="DEnricher", def=function(params, M, fit_options, gene_sets) {
    standardGeneric("DEnricher")
})

#'@importClassesFrom GSEABase GeneSetCollection
#'@importFrom futile.logger flog.info
#'@importFrom GSEABase GeneSetCollection geneIds geneIds<- setIdentifier
#'setIdentifier<-
#'@include FitOptions.R
#'@include Genesets.R
#'@include IGSAinput.R
#'@include SEAparams.R
#'@include SEAres.R
# params <- seaParams(igsaInput); M <- expr_data; 
# gene_sets <- merged_gene_sets
setMethod(
    f="DEnricher",
    signature=c("SEAparams", "ExprData", "FitOptions", "GeneSetCollection"),
    definition=function(params, M, fit_options, gene_sets) {
        if (length(de_genes(params)) == 1 && is.na(de_genes(params)[[1]])) {
            # if there are not set the DE genes then calculate them
            dif <- igsaGetDEGenes(params, M, fit_options);
        } else {
            # otherwise we use them
            dif <- de_genes(params);
            flog.info(paste("DE genes", length(dif)));
        }
        
        br <- br(params);
        # get all the genes present in the gene sets
        allGenes <- unique(unlist(asList(gene_sets)));
        
        if (length(br) > 1) {
            # if its a user provided br, then filter the genes that are in the
            # gene sets (statistical importance)
#             br <- intersect(br, allGenes);
            flog.info(paste("Using user provided BR:", length(br), "genes."));
            if (length(br) < 2) {
                stop("No genes in br after intersecting with experiment genes");
            }
        } else if (br == "bri") {
            # bri is the genome, however we use gene sets genes as genome
            # (statistical importance)
            br <- allGenes;
            flog.info(paste("Using BRI:", length(br), "genes."));
        } else if (br == "briii") {
            # briii are the reliably detected genes by the experiment
            br <- unique(rownames(M));
            flog.info(paste("Using BRIII:", length(br), "genes."));
        } else {
            stop("Incorrect br option.");
        }
        
        # DE genes must be in the br (statistical importance)
        dif <- intersect(dif, br);
        
        # so we filter gene sets that dont have genes in this br
        validGSets <- lapply(gene_sets, function(actGs) {
            validGenes <- intersect(geneIds(actGs), br);
            if (length(validGenes) == 0) {
                return(NULL);
            }
            actGsId <- setIdentifier(actGs);
            geneIds(actGs) <- validGenes;
            setIdentifier(actGs) <- actGsId;
            return(actGs);
        });
        validGSets <- validGSets[ !unlist(lapply(validGSets, is.null)) ];
        gene_sets <- GeneSetCollection(validGSets);
        
        test <- test(params);
        
        # and run dEnricher
        seaRes <- runDEnricher(dif, gene_sets, br, test=test);
        
        seaRes <- SEAres(gene_sets_res=seaRes);
        
        return(seaRes);
    }
)

setGeneric(name="runDEnricher",
    def=function(dif, genesets, br, test) {
    standardGeneric("runDEnricher")
})

#'@importClassesFrom GSEABase GeneSetCollection
#'@importFrom BiocParallel bplapply bpparam
#'@importFrom futile.logger flog.info
#'@importFrom GSEABase geneIds setIdentifier setName
#'@importFrom stats fisher.test pbinom phyper
#'@include GenesetRes.R
#'@include Genesets.R
# genesets <- gene_sets; test="FisherTest"
setMethod(
    f="runDEnricher",
    signature=c("character", "GeneSetCollection", "character", "character"),
    definition=function(dif, genesets, br, test) {
        # most of this function is copied from dEnricher function
        
        GeneID <- dif;
        genes.group <- GeneID[!is.na(GeneID)];
        
        gs <- genesets;
        nSet <- length(gs)
        
        # defining the three tests
        doFisherTest <- function(genes.group, genes.term, genes.universe) {
            genes.hit <- intersect(genes.group, genes.term)
            X <- length(genes.hit)
            K <- length(genes.group)
            M <- length(genes.term)
            N <- length(genes.universe)
            cTab <- matrix(c(X, K - X, M - X, N - M - K + X), nrow = 2, 
                dimnames = list(c("anno", "notAnno"),
                                c("group", "notGroup")))
            
            if (any(cTab < 0)) return(1);
            p.value <- ifelse(all(cTab == 0), 1, stats::fisher.test(cTab, 
                                            alternative = "greater")$p.value)
            return(p.value)
        }
        doHypergeoTest <- function(genes.group, genes.term, genes.universe) {
            genes.hit <- intersect(genes.group, genes.term)
            X <- length(genes.hit)
            K <- length(genes.group)
            M <- length(genes.term)
            N <- length(genes.universe)
            x <- X
            m <- M
            n <- N - M
            k <- K
            if (any(c(x, m, n, k) < 0)) return(1);
            p.value <- ifelse(m == 0 || k == 0, 1, stats::phyper(x, 
                                    m, n, k, lower.tail=FALSE, log.p=FALSE))
            return(p.value)
        }
        doBinomialTest <- function(genes.group, genes.term, genes.universe) {
            genes.hit <- intersect(genes.group, genes.term)
            X <- length(genes.hit)
            K <- length(genes.group)
            M <- length(genes.term)
            N <- length(genes.universe)
            p.value <- ifelse(K == 0 || M == 0 || N == 0, 1, stats::pbinom(X, 
                                    K, M/N, lower.tail=FALSE, log.p=FALSE))
            if (is.na(p.value)) return(1);
            return(p.value)
        }
        
        genes.universe <- br;
        genes.group <- intersect(genes.universe, genes.group)
        if (length(genes.group) == 0) {
            warning(paste("There are no differentialy expressed genes",
                "being used."));
        }
        
        flog.info(paste("Running SEA at cores:", bpparam()$workers));
        # I am passing MIGSA's not exported functions to bplapply to avoid
        # SnowParam environment errors
        seaRes <- bplapply(gs, function(actGset, GenesetRes) {
        #     seaRes <- lapply(gs@gene_sets, function(actGset) {
            genes.term <- unique(unlist(geneIds(actGset)))
            p.value <- switch(test,
                FisherTest   = doFisherTest(genes.group, genes.term,
                                                            genes.universe),
                HypergeoTest = doHypergeoTest(genes.group, genes.term,
                                                            genes.universe),
                BinomialTest = doBinomialTest(genes.group, genes.term,
                                                            genes.universe))
            
            # important genes are the DE that are in the term
            impGenes <- intersect(genes.group, genes.term);
            
            # We put name as setIdentifier and id as setName because we 
            # disagree with GSEABase's usage as they use Names as unique, and 
            # identifiers not.
            actRes <- GenesetRes(id=setName(actGset), 
                                name=setIdentifier(actGset), pvalue=p.value, 
                                genes=geneIds(actGset),
                                enriching_genes=impGenes);
            return(actRes);
        }, GenesetRes=GenesetRes)
        #     })

        return(seaRes);
    }
)
