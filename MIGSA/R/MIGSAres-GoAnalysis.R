# to avoid R CMD check errors we set them as NULL
is_GO = NULL;

#'Explore Gene Ontology gene sets in MIGSAres
#'
#'\code{migsaGoTree} plots the GO tree/s present in migsaRes.
#'\code{getHeights} returns the heights of given a list of ids (GO IDs).
#'
#'@param migsaRes MIGSAres object. It must contain at least one GO gene set.
#'@param ids character vector indicating the queried GO ids.
#'@param minHeight logical indicating if the minimum or maximum height must be 
#'calculated. If it is FALSE then the longest path to the root is calculated, 
#'otherwise, the shortest path.
#'@param categories vector. Each experiment category, it will print different 
#'node color for each. Can have NAs. Must have length equal to number of 
#'experiments, i.e. length(categories) == ncol(migsaRes)-3
#'@param categColors character. Color for each category. Must have the same 
#'length as the number of different categories.
#'@param ont character. One of "BP", "CC" or "MF". Selected ontology to plot.
#'@param legendPos . Parameter passed to legend function.
#'@param treeInfo . Data.frame with GO ids as rownames, and three columns:
#'Enriched (logical), Important (logical) and Color (character). If Enriched is
#'true then it will be ploted, if Important is true then it will have another 
#'shape, Color has the color name to use.
#'@param legends . Matrix with two columns, each col is a pair (lengend, color).
#'@param ... not in use.
#'
#'@return If migsaGoTree: A list with the used data to plot. If getHeights: A 
#'list with each term height.
#'
#'@docType methods
#'@name MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'
#'@examples
#'## Lets load breast cancer results.
#'data(bcMigsaRes);
#'
#'###### migsaGoTree
#'## Get the first 40 Gene Ontology gene sets results from CC.
#'goRes <- bcMigsaRes[bcMigsaRes$GS_Name == "CC",];
#'fst40goRes <- goRes[1:40,];
#'
#'## And lets plot the results GO trees.
#'\dontrun{
#'aux <- migsaGoTree(fst40goRes);
#'}
#'
#'###### getHeights
#'## Get the first 40 Gene Ontology gene sets IDs.
#'goIds <- bcMigsaRes[bcMigsaRes$GS_Name %in% c("BP", "CC", "MF"), "id"];
#'fst40goIds <- goIds[1:40,];
#'
#'\dontrun{
#'## And lets get the heights in the GO tree structure.
#'getHeights(fst40goIds);
#'}
#'
setGeneric(name="MIGSAres-GOanalysis", def=function(migsaRes) {
    standardGeneric("MIGSAres-GOanalysis")
})

#'@name migsaGoTree
#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases migsaGoTree,MIGSAres-method
#'@exportMethod migsaGoTree
#'
setGeneric(name="migsaGoTree", def=function(migsaRes, ...) {
    standardGeneric("migsaGoTree")
})

#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases migsaGoTree,MIGSAres-method
#'
#'@importFrom AnnotationDbi Ontology
#'@importFrom futile.logger flog.info
#'@importFrom graphics legend
#'@importFrom grDevices col2rgb rgb
#'@importFrom utils combn
#'@include MIGSAres.R
#'@include GoAnalysis.R
#'
setMethod(
    f="migsaGoTree",
    signature=c("MIGSAres"),
    definition=function(migsaRes, categories=rep(NA, ncol(migsaRes)-3),
        categColors="red", ont="BP", legendPos="topleft") {
        stopifnot(validObject(migsaRes));
        stopifnot(length(categories) == ncol(migsaRes)-3);
        stopifnot(length(unique(categories))
            -any(is.na(categories))+all(is.na(categories))
            == length(unique(categColors)));
        stopifnot(ont %in% c("BP", "CC", "MF"));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        if (any(!is.na(categories))) {
            migsaRes <- migsaRes[,c(1:3, 3+which(!is.na(categories)))];
            categories <- categories[!is.na(categories)];
        } else {
            categories[is.na(categories)] <- "Enriched";
        }
        
        # get gene set names and if they are GO
        goGsets <- unique(as.data.frame(migsaResAll(migsaRes)[,
                                            list(gene_set_name, is_GO)]));
        
        if (!any(goGsets[ goGsets$gene_set_name == ont, "is_GO" ])) {
            stop("No Gene Ontology gene set to plot.");
        }
        
        # and filter the MIGSAres object with only the GO gene sets
        actRes <- migsaRes[ migsaRes[, "GS_Name", drop=!FALSE] == ont, ];
        
        # get the gene sets and the number of enriched experiments
        plotRes <- data.frame(id=actRes[,1], gs_name=actRes[,3]);
        
        # get the real ontology, I could use GS_Name, but this is more trustable
        plotRes <- cbind(plotRes,
                        ont=Ontology(as.character(plotRes$id)));
        
        actRes <- actRes[!is.na(plotRes$ont),];
        plotRes <- plotRes[!is.na(plotRes$ont),];
        ontsPresent <- unique(as.character(plotRes$ont));
        
        # for each condition, the times it gets enriched
        enrichedPerc <- do.call(cbind, by(t(actRes[,-(1:3)]), 
            categories, colMeans, na.rm=TRUE));
        enrichedPerc[is.nan(enrichedPerc)] <- 0;
        
        uniqColors <- lapply(categColors, function(x) 
            t(col2rgb(x)/255));
        
        plotRes$colors <- apply(enrichedPerc, 1, function(actRow) {
            colors <- unlist(lapply(seq_along(actRow), function(i) {
                # we put as alpha the percentage of enrichment
                rgb(uniqColors[[i]], alpha=actRow[[i]]);
            }))
            rgb(t(mixColors(colors)));
        })
        
        out <- data.frame(matrix(!FALSE, nrow=nrow(plotRes), ncol=3));
        colnames(out) <- c("Enriched", "Important", "Color");
        rownames(out) <- plotRes$id;
        
        out$Important <- FALSE;
        out$Color <- plotRes$colors;
        
        # lets create the legend
        categories <- colnames(enrichedPerc);
        legends <- do.call(cbind, lapply(seq_along(categories), function(i) {
            apply(combn(categories, i), 2, function(x) {
                aux <- mixColors(categColors[categories %in% x]);
                return(c(paste(x, collapse=" "), 
                    rgb(t(aux))));
            })
        }))
        
        acttree <- goTree(out, ont, legends, legendPos);
        
        return(invisible(acttree));
    }
)

#'@name goTree
#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases goTree,MIGSAres-method
#'@exportMethod goTree
#'
setGeneric(name="goTree", def=function(treeInfo, ...) {
    standardGeneric("goTree")
})

#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases goTree,MIGSAres-method
#'
#'@importFrom graphics legend
#'@include GoAnalysis.R
#'
setMethod(
    f="goTree",
    signature=c("data.frame"),
    definition=function(treeInfo, ont="BP", legends=NA, legendPos="topleft") {
        stopifnot(ont %in% c("BP", "CC", "MF"));
        stopifnot(colnames(treeInfo) == c("Enriched", "Important", "Color"));
        stopifnot(unlist(lapply(seq_len(ncol(treeInfo)), function(i) {
            class(treeInfo[,i])
        })) == c("logical", "logical", "character"))
        
        treeInfo <- treeInfo[Ontology(rownames(treeInfo)) == ont, ];
        treeData <- list(treeInfo);
        names(treeData) <- ont;
        
        graph <- createGoGraph(treeData);
        acttree <- list(graph=graph, gotree=treeData);
        
        # select grid depending on number of ontologies to plot
        acttree$graph <- plotGoTree(acttree, ont);
        
        if (!all(is.na(legends)))
            legend(legendPos, legends[1,], fill=legends[2,], cex=.75);
        return(invisible(acttree));
    }
)

#'@name getHeights
#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases getHeights,character-method
#'@exportMethod getHeights
#'
setGeneric(name="getHeights", def=function(ids, ...) {
    standardGeneric("getHeights")
})

#'@inheritParams MIGSAres-GOanalysis
#'@rdname MIGSAres-GoAnalysis
#'@aliases getHeights,character-method
#'
#'@importFrom GO.db GOBPPARENTS GOCCPARENTS GOMFPARENTS
#'@include GoAnalysis.R
#'
setMethod(
    f="getHeights",
    signature=c("character"),
    definition=function(ids, minHeight=TRUE) {
        allParents <- as.list(GOBPPARENTS);
        allParents <- c(allParents, as.list(GOMFPARENTS));
        allParents <- c(allParents, as.list(GOCCPARENTS));
        
        ids <- gsub(" ", "", ids);
        
        # for each GO id find its height
        result <- lapply(ids, function(x) {
            getHeight(x, minHeight, allParents)
        });
        
        return(unlist(result));
    }
)
