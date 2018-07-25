# # library(RColorBrewer);

# percentages: two column dataframe with the percentage of enriched datasets 
# for each method (SEA, mGSZ)
# nColors: the number of opacities used for each one of the 3 colors
# giveColors <- function(percentages, nColors) {
#   nColors <- ceiling(nColors);
#   
#   SEA_colors <- brewer.pal(max(3,nColors), "Blues");
#   mGSZ_colors <- brewer.pal(max(3,nColors), "Reds");
#   both_colors <- brewer.pal(max(3,nColors), "Greens");
#   colors <- rep("white", nrow(percentages));
#   names(colors) <- rownames(percentages);
#   
#   both_rows <- apply(percentages>0, 1, function(x) (x[1] && x[2]));
#   mGSZ_rows <- apply(percentages>0, 1, function(x) (!x[1] && x[2]));
#   SEA_rows <- apply(percentages>0, 1, function(x) (x[1] && !x[2]));
#   
#   colors[both_rows] <- both_colors[floor(rowSums(percentages)[both_rows]/2)];
#   colors[mGSZ_rows] <- mGSZ_colors[ percentages[mGSZ_rows, "mGSZ"] ];
#   colors[SEA_rows] <- SEA_colors[ percentages[SEA_rows, "SEA"] ];
#   
#   return(colors);
# }

# convert2percentages <- function(enrichedCounts, percentages, nDatasets) {
#   percentages <- sort(percentages);
#   minEnricheds <- percentages*nDatasets/100;
#   result <- enrichedCounts;
#   result[ enrichedCounts < minEnricheds[1] ] <- 0;
#   
#   for(perc in (1:length(percentages))) {
#     result[ enrichedCounts >= minEnricheds[perc] ] <- perc;
#   };
#   return(result);
# }

# createtree <- function(data2plot, percentages, importantGO=NULL) {
#   dataResult <- lapply(names(data2plot), function(dataSetName) {
#     dataSet <- data2plot[[ dataSetName ]];
#     
#     methodResult <- lapply(names(dataSet), function(methodName) {
#       method <- dataSet[[ methodName ]];
#       Reduce(rbind, lapply(names(method), function(ontology) {
#         expResult <- method[[ ontology ]];
#         return(cbind(expResult[,c("NAME", "Enriched")], ont=ontology));
#       }));
#     })
#     dataSetResult <- Reduce(function(...) merge(...,by=c("NAME", "ont"), 
#             all=!F), methodResult);
#     colnames(dataSetResult)[3:4] <- paste(dataSetName, names(dataSet),
#             sep="_");
#     return(dataSetResult);
#   });
#   dataResult <- Reduce(function(...) merge(...,by=c("NAME", "ont"), all=!F), 
#                             dataResult);
#   
#   
#   ontologies <- c("BP", "CC", "MF");
#   treeData <- lapply(ontologies, function(actualOnt) {
#     actualResults <- subset(dataResult, ont==actualOnt, select=-2);
#     rownames(actualResults) <- actualResults$NAME;
#     actualResults <- actualResults[,-1];
#     
#     #   actualResults <- actualResults[
#                     rowSums(actualResults, na.rm=!F) > 0, ];
#     actualResults[is.na(actualResults)] <- F;
#     
#     actualSEAResults <- actualResults[, grep("SEA", 
#                 colnames(actualResults))];
#     actualSEAResults <- rowSums(actualSEAResults);
#     
#     actualmGSZResults <- actualResults[, grep("mGSZ", 
#                 colnames(actualResults))];
#     actualmGSZResults <- rowSums(actualmGSZResults);
#     
#     actualResults <- cbind(SEA=actualSEAResults, mGSZ=actualmGSZResults);
#     actualResults <- convert2percentages(actualResults, percentages, 
#                                         length(data2plot));
#     
#     out <- data.frame(matrix(NA, nrow=nrow(actualResults), ncol=3));
#     colnames(out) <- c("Enriched", "Important", "Color");
#     rownames(out) <- rownames(actualResults);
#     
#     
#     out$Enriched <- rowSums(actualResults) > 0;
#     out$Important <- rownames(out) %in% importantGO;
#     out$Color <- giveColors(actualResults, length(percentages));
#     
#     return(out);
#   });
#   names(treeData) <- ontologies;
#   return(treeData);
# }

# terms: list of GO Ids to plot
# setGeneric(name="createTreeFromList", def=function(terms, ...) {
#     standardGeneric("createTreeFromList")
# })

# internal function
# #'@importFrom AnnotationDbi Ontology
# #'@importFrom futile.logger flog.info
# #'@importFrom GO.db GOTERM
# setMethod(
#     f="createTreeFromList",
#     signature=c("character"),
#     definition=function(terms, color="red") {
#         stopifnot(length(terms) > 0);
#         
#         terms <- cbind(terms, ont=Ontology(GOTERM)[terms]);
#         terms <- terms[ !is.na(rownames(terms)), ];
#         
#         ontsPresent <- names(table(terms[,"ont"]));
#         flog.info(ontsPresent);
#         
#         treeData <- lapply(ontsPresent, function(actOnt) {
#             actualResults <- terms[ terms[,"ont"] == actOnt, ];
#             out <- data.frame(matrix(!FALSE, nrow=nrow(actualResults),
#                                 ncol=3));
#             colnames(out) <- c("Enriched", "Important", "Color");
#             rownames(out) <- rownames(actualResults);
#             
#             out$Important <- FALSE;
#             out$Color <- color;
#             
#             return(out);
#         });
#         names(treeData) <- ontsPresent;
# 
#         graph <- createGoGraph(treeData);
#         return(list(graph=graph, gotree=treeData));
#     }
# )

setGeneric(name="createGoGraph", def=function(ontologies) {
    standardGeneric("createGoGraph")
})

#'@importFrom AnnotationDbi Ontology
#'@importFrom GO.db GOBPPARENTS GOCCPARENTS GOMFPARENTS
#'@importFrom GOstats GOGraph
#'@importFrom graph nodes removeNode
#'@importFrom stats na.omit
setMethod(
    f="createGoGraph",
    signature=c("list"),
    definition=function(ontologies) {
        GOparents <- c(BP=GOBPPARENTS, CC=GOCCPARENTS, MF=GOMFPARENTS)
        
        # For each ontology create the GO graph
        GOgraph <- lapply(ontologies, function(ontology, GOparents) {
            graph <- GOGraph(
                rownames(ontology)[ontology$Enriched | ontology$Important],
                GOparents[[as.character(na.omit(unique(
                                        Ontology(row.names(ontology)))))]]
            )
            
            # Remove "all" node if present
            if("all" %in% nodes(graph)){
                graph <- removeNode("all", graph);
            }
            
            # invert the graph (if not it plots it 180 degrees)
            graph <- nodeLevel(invertGraph(graph));
            return(graph);
        }, GOparents);
        
        return(GOgraph)
    }
)

setGeneric(name="invertGraph", def=function(graph) {
    standardGeneric("invertGraph")
})

#'@importClassesFrom graph graphNEL
#'@importFrom graph addEdge edges nodes
#'@importFrom RBGL isomorphism
setMethod(
    f="invertGraph",
    signature=c("graphNEL"),
    definition=function(graph) {
        invEdges <- do.call(rbind, 
            lapply(names(edges(graph)), function(nodeto) {
            # for each edge, get its starting node
            do.call(rbind, lapply(edges(graph)[[nodeto]], function(nodefrom) {
                # for each destination, return the other way round
                return(c(nodefrom, nodeto));
            }));
        }));
        fromNodes <- unique(invEdges[,1]);
        
        # return it as a list
        invEdges <- lapply(fromNodes, function(nodefrom) {
            invEdges[ invEdges[,1] == nodefrom, 2];
        });
        names(invEdges) <- lapply(invEdges, function(x) names(x)[[1]]);
        
        # dont give the edges to graphNEL, so it creates also the edges of the 
        # nodes with no children
        ginv <- new("graphNEL", edgemode="directed", nodes=nodes(graph));
        newEdges <- edges(ginv);
        # over write edges that we inverted
        newEdges[ names(invEdges) ] <- invEdges;
        
        # create the new graph with the complete edges structure
        ginv <- new("graphNEL", edgemode="directed", nodes=nodes(graph),
            edgeL=newEdges);
        
        # control
        stopifnot(isomorphism(graph,ginv)$isomorphism);
        return(ginv)
    }
)

setGeneric(name="nodeLevel", def=function(graph, ...) {
    standardGeneric("nodeLevel")
})

#'@importClassesFrom graph graphNEL
#'@importFrom graph degree
#'@importFrom Rgraphviz 'nodeData<-' 'nodeDataDefaults<-'
setMethod(
    f="nodeLevel",
    signature=c("graphNEL"),
    definition=function(graph, root=NULL) {
        # default attribute "-1"
        nodeDataDefaults(graph, "level") <- -1;
        
        # inicialization
        if (is.null(root)) {
            root <- names(which(degree(graph)$inDegree == 0));
            if (length(root) > 1){
                warning("Unmatching GO ids found")
                root <- root[which(degree(graph)$outDegree[root] ==
                                    max(degree(graph)$outDegree[root]))];
            }
        }
        nodeData(graph, root, "level") <- 0;
        return(graph)
    }
)

# doGOtree <- function(data2plot, percentages=c(50, 75, 100), 
#     importantGO=NULL) {
#   gotree <- createtree(data2plot, percentages, importantGO);
#   flog.info(checkEnrichment(gotree));
#   
#   graph <- createGoGraph(gotree);
#   
#   return(list(graph=graph, gotree=gotree));
# }

setGeneric(name="checkEnrichment", def=function(ontologies) {
    standardGeneric("checkEnrichment")
})

setMethod(
    f="checkEnrichment",
    signature=c("list"),
    definition=function(ontologies) {
        # Check for enrichment presence
        do.call(c, lapply(ontologies, function(onto) {
            sum(onto$Enriched);
        }));
    }
)

setGeneric(name="plotGoTree", def=function(graph, ont) {
    standardGeneric("plotGoTree")
})

setMethod(
    f="plotGoTree",
    signature=c("list", "character"),
    definition=function(graph, ont) {
        gotree <- graph$gotree;
        GOgraph <- graph$graph;
        
        # name.height: node height
        # name.width: node width
        # name.fontsize: node name size
        # name.lwd: vert spessor
        # this is for pretty printing, well the best I could do :P
        # vertSpessor: connectors and outter circles
        vertSpessor <- 50/checkEnrichment(gotree)[[ont]];
        fontSize <- max(50, checkEnrichment(gotree)[[ont]] / 16);
        size <- c(name.legend.pt=8, name.legend.text=3, name.height=1.9,
                name.width=0.9, name.fontsize.full=11, name.fontsize=fontSize,
                name.lwd=vertSpessor);
        
        plotGo(GOgraph[[ont]], gotree[[ont]], lroot=ont, size);
    }
)

setGeneric(name="plotGo", def=function(graph, ontology, ...) {
    standardGeneric("plotGo")
})

#'@importFrom AnnotationDbi Term
#'@importFrom GO.db GOTERM
#'@importClassesFrom graph graphNEL
#'@importFrom graph 'edgeRenderInfo<-' 'nodeRenderInfo<-' 'graphRenderInfo<-'
#'@importFrom Rgraphviz layoutGraph nodeData renderGraph
setMethod(
    f="plotGo",
    signature=c("graphNEL", "data.frame"),
    definition=function(graph, ontology, lroot="MF",
        size=c(name.lwd=3, name.height=1.9, name.width=0.9,
        name.fontsize=30)) {
        
        # only nodes to plot
        plotNodes <- rownames(ontology)[ontology$Enriched | 
                                        ontology$Important];
        
        # change the id for the Name
        terms <- sapply(plotNodes, function(id) {
            if(!is.null(GOTERM[[id]])) {
                Term(GOTERM[[id]]);
            } else {
                id; # database version issue!!!!
            }
        })
        terms <- gsub("\nor"," or",
                    gsub("cell\n","cell ",
                        gsub("\nof"," of",
                            gsub(" ","\n", terms))));
        
        # add GO root name
        root <- unlist(nodeData(graph, attr="level"))[
                    unlist(nodeData(graph,attr="level")) == 0];
        root[1] <- lroot; 
        terms <- c(terms,root);
        
        # color look up
        color <- ontology$Color;
        names(color) <- rownames(ontology);
        if (!(names(root[1]) %in% names(color))) {
            color <- c("white", color);
            names(color)[1] <- names(root[1]);
        }
        color <- color[names(terms)];
        
        # if no letters are needed, just a circule for the nodes
        stopifnot(size[["name.fontsize.full"]] > 0)
        
        # size
        width <- sapply(terms, function(term) {
            max(nchar((unlist(strsplit(term, split="\n")))))
        })*size["name.width"];
        
        height <- sapply(terms,function(term) {
            length(unlist(strsplit(term, split="\n")))
        })*size["name.height"];
        
        # not in use right now, but if we have important nodes then change the 
        # shape
        importantGO <- rownames(ontology[ontology$Important,]);
        shapes <- ifelse(names(terms) %in% importantGO, "box", "ellipse");
#         shape: The shape of the node. Current acceptable values are circle, 
#         rectangle, rect, box and ellipse. The circle shape is the default. 
#         Note that box, rect and rectangle all correspond to the same actual 
#         shape
        names(width) <- names(height) <- names(shapes) <- names(terms);
        
        nodeAttrs <- list(width=width, height=height, shape=shapes);
        
        node <- list();
        node$node <- list(fixedsize=FALSE);
        graph <- layoutGraph(graph, nodeAttrs=nodeAttrs, attrs=node);
        edgeRenderInfo(graph) <- list(lwd=size[["name.lwd"]], 
                                        arrowhead="none");
# other possible ways of encoding info
# color: Basic drawing color for the node, corresponding to the outside edge 
# of the node.
# fontcolor: Color used for text. This defaults to black.
# fontname: Font used for text. The default of this is Times Roman.
        # others
        nodeRenderInfo(graph) <- list(label=terms, fill=color,
                                    fontsize=size[["name.fontsize"]],
                                    lwd=size[["name.lwd"]]/2);
        
        graphRenderInfo(graph) <- list(main=lroot);
        # rendering at last
        graph <- renderGraph(graph);
        
        return(graph)
    }
)

setGeneric(name="mixColors", def=function(colorList) {
    standardGeneric("mixColors")
})

setMethod(
    f="mixColors",
    signature=c("character"),
    definition=function(colorList) {
    #     if (length(colorList) == 1) {
    #         return(colorList);
    #     }
    #     colMeans(do.call(rbind, lapply(colorList, function(actColor) {
    #         newColor <- t(col2rgb(actColor, alpha=TRUE)/255);
    #         newColor <- newColor * newColor[, "alpha"];
    #         return(newColor[,1:3])
    #     })))
        newColors <- colSums(do.call(rbind, 
            lapply(colorList, function(actColor) {
                newColor <- t(col2rgb(actColor, alpha=TRUE)/255);
                newColor <- newColor * (newColor[, "alpha"] != 0);
                return(newColor[,1:3])
        })))
        unlist(lapply(newColors, min, 1));
    }
)

setGeneric(name="getHeight", def=function(term, minHeight, allParents) {
    standardGeneric("getHeight")
})

setMethod(
    f="getHeight",
    signature=c("character", "logical", "list"),
    definition=function(term, minHeight, allParents) {
        # bp, mf and cc ids
        fstTerms <- c("GO:0008150", "GO:0003674", "GO:0005575");
        actualHeight <- 1;
        actualParents <- allParents[[term]];
        
        if (term %in% fstTerms) {
            return(0);
        } else if (is.null(actualParents)) {
            return(NA);
        }
        
        if (minHeight) {
            # while we dont get to the root node sum 1 to the result.
            # remember we are starting from our node, finding its parents 
            # pseudo recursively
            while ((!any(fstTerms %in% actualParents)) && 
                    (actualParents != "all")) {
                actualHeight <- actualHeight +1;
                actualParents <- unique(unlist(allParents[actualParents]));
            }
        } else {
            while (!all(actualParents %in% as.list(c("all", fstTerms)))) {
                actualHeight <- actualHeight +1;
                actualParents <- unique(unlist(allParents[actualParents]));
            }
        }
        
        return(actualHeight);
    }
)
