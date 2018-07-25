# to avoid R CMD check errors we set them as NULL
Cols = Rows = number = gene = geneset = Count = Value = x = y = NULL;

#'MIGSAres plots
#'
#'\code{genesHeatmap} plots a heatmap with the number of experiments in which 
#'each gene contributed to enrich each gene set.
#'\code{genesBarplot} generates a barplot of the number of gene sets in which 
#'each gene contributed to enrich. Each gene set counts 1 regardless if it was 
#'enriched in many experiments. x-axis each gene, y-axis number of gene sets 
#'in which it contributed to enrich.
#'\code{migsaHeatmap} plots the enrichment heatmap of the MIGSAres object.
#'\code{geneSetBarplot} generates a barplot of the number of experiments in 
#'which each gene set was enriched. x-axis each gene set, y-axis times it was 
#'enriched (0 to #experiments).
#'
#'@param migsaRes MIGSAres object.
#'@param categories list. List of character vectors, each vector must have the
#'same length as the number of experiments (ncol(migsaRes)-3). It will plot 
#'for each category/column a color representing its category.
#'@param categLabels logical. Indicates if labels should be plotted for each 
#'category.
#'@param colPal vector. Character vector of two colors, first value will 
#'represent FALSE/1 on heatmap, and second value TRUE/0.
#'@param enrFilter numeric. Keep gene sets enriched in at least enrFilter 
#'experiments.
#'@param gsFilter numeric. Keep genes enriched in at least gsFilter gene sets.
#'@param expFilter numeric. Keep experiments which enriched at least expFilter 
#'gene sets.
#'@param col.dist character. Distance algorithm to be used in columns, passed 
#'to vegdist function. If migsaRes has cutoff then default is jaccard, else, 
#'default is euclidean.
#'@param row.dist character. Distance algorithm to be used in rows, passed to 
#'vegdist function. If migsaRes has cutoff then default is jaccard, else, 
#'default is euclidean.
#'@param remove0Rows logical. Whether remove gene sets that are not enriched in 
#'any experiment.
#'@param breaks numeric. If migsaRes does not have cutoff then break P-values 
#'in breaks intervals.
#'@param dendrogram character. Wheter to plot "row", "col", "none" or "both"
#'dendrograms.
#'@param layout matrix. The layout used for heatmaps.
#'@param ... not in use.
#'
#'@return In heatmap functions: A list returned by heatmap.2 function 
#'(plotted data). In other functions: A ggplot object used as graphic.
#'
#'@docType methods
#'@name MIGSAres-plots
#'@rdname MIGSAres-plots
#'@seealso \code{\link{setEnrCutoff}}
#'
#'@examples
#'data(migsaRes);
#'## For this example lets work with the first 50 results
#'migsaRes <- migsaRes[1:50,];
#'
#'#### genesHeatmap
#'## First lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets check what genes contributed to enrich the highest number of gene 
#'## sets (in more than one gene set).
#'genesHeatmap(migsaResWCoff, gsFilter=1);
#'
#'## Moreover we can keep gene sets which where enriched in more than 
#'## enrFilter experiments. To do this, we can use the enrFilter parameter.
#'
#'#### genesBarplot
#'## Lets set a cutoff of 0.01
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.01);
#'
#'## Lets check what genes contributed to enrich the highest number of gene 
#'## sets (in more than one gene set).
#'genesBarplot(migsaResWCoff, gsFilter=1);
#'
#'## Moreover we can keep gene sets which where enriched in more than 
#'## enrFilter experiments. To do this, we can use the enrFilter parameter.
#'
#'#### migsaHeatmap
#'## Lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets visually check enriched gene sets shared between experiments.
#'migsaHeatmap(migsaResWCoff);
#'
#'#### geneSetBarplot
#'## Lets set a cutoff of 0.1
#'migsaResWCoff <- setEnrCutoff(migsaRes, 0.1);
#'
#'## Lets check in how many experiments each gene set was enriched (in more 
#'## than one experiment).
#'geneSetBarplot(migsaResWCoff, enrFilter=1);
#'
setGeneric(name="MIGSAres-plots", def=function(migsaRes) {
    standardGeneric("MIGSAres-plots")
})

#'@name genesHeatmap
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesHeatmap,MIGSAres-method
#'@exportMethod genesHeatmap
#'
setGeneric(name="genesHeatmap", def=function(migsaRes, ...) {
    standardGeneric("genesHeatmap")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesHeatmap,MIGSAres-method
#'
#'@importFrom reshape2 melt
#'@include MIGSAres-genesManipulation.R
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(f="genesHeatmap",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, gsFilter=0, 
        colPal=c("white", "red"), col.dist=NA, row.dist=NA, layout=NA,
        dendrogram="row") {
        stopifnot(dendrogram %in% c("none", "row", "col", "both"));
        ## todo: agregar funcionalidad a dendrogram
        ## todo: ver que el dendrograma no pega bien con las celdas
        stopifnot(length(colPal) == 2);
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=TRUE) >
                                    enrFilter, ];
        plotGenes <- genesInSets(actRes);
        
        # keep genes enriched in more than gsFilter gene sets
        plotGenes <- plotGenes[, colSums(plotGenes > 0) > gsFilter];
        plotGenes <- plotGenes[ rowSums(plotGenes) > 0 ,];
        
        if (is.na(col.dist)) {
            col.dist="jaccard";
        }
        if (is.na(row.dist)) {
            row.dist="jaccard";
        }
        plotInfo <- melt(plotGenes);
        colnames(plotInfo) <- c("geneset", "gene", "Count");
        
        ddc.s <- suppressWarnings(
            vegdist(t(plotGenes), col.dist, na.rm=TRUE));
        ddc.s[is.na(ddc.s)] <- 0;
        hc.s <- hclust(ddc.s, method="average");
        colInd <- hc.s$order;
        
        # jaccard clustering per row (gene set)
        ddr.s <- suppressWarnings(
            vegdist(plotGenes, row.dist, na.rm=TRUE));
        ddr.s[is.na(ddr.s)] <- 0;
        hr.s <- hclust(ddr.s, method="average");
        rowInd <- hr.s$order;
        
        # reordering the plot by distances
        plotInfo$gene <- factor(as.character(plotInfo$gene), 
            levels=colnames(plotGenes)[colInd]);
        plotInfo$geneset <- factor(as.character(plotInfo$geneset), 
            levels=rownames(plotGenes)[rowInd]);
        
        p_heatm <- ggplot(plotInfo, aes(gene, geneset));
        p_heatm <- p_heatm + geom_tile(aes(fill=Count));
        p_heatm <- p_heatm + xlab(NULL) + ylab(NULL);
        p_heatm <- p_heatm + scale_y_discrete(breaks=NULL);
        p_heatm <- p_heatm + scale_fill_gradient(low=colPal[[2]], 
            high=colPal[[1]]);
        p_heatm <- p_heatm + theme(legend.position="bottom", 
                            legend.key.size=unit(8, "points"), 
                            legend.text=element_text(size=unit(8, "points")), 
                            plot.margin=unit(c(0,0,0,0),"mm"),
                            axis.text.x=element_text(angle = 45, hjust=1));
        p_heatm <- p_heatm + guides(fill=guide_legend(nrow=2));
        
        
        allPlots <- list();
        # now the x dendro
        if (dendrogram %in% c("col", "both")) {
            p_x_dendro <- ggdendrogram(hr.s, labels=FALSE);
            p_x_dendro <- p_x_dendro + xlab(NULL) + ylab(NULL)
            p_x_dendro <- p_x_dendro + theme(axis.text.x=element_blank(),
                                    axis.text.y=element_blank(), 
                                    plot.margin=unit(c(0,0,0,0),"mm"));
            p_x_dendro <- p_x_dendro + coord_flip();
            p_x_dendro <- p_x_dendro + scale_y_continuous(expand=c(0.005,0), 
                trans="reverse");
            p_x_dendro <- p_x_dendro + scale_x_continuous(expand=
                c(1/nrow(plotGenes)/2,0));
            allPlots <- list(p_x_dendro);
        }
        
        # now the y dendro
        if (dendrogram %in% c("row", "both")) {
            p_y_dendro <- ggdendrogram(hc.s, labels=FALSE);
            p_y_dendro <- p_y_dendro + scale_y_continuous(expand=c(0.005,0));
            p_y_dendro <- p_y_dendro + scale_x_continuous(expand=
                c(1/ncol(plotGenes)/2,0));
            p_y_dendro <- p_y_dendro + xlab(NULL) + ylab(NULL)
            p_y_dendro <- p_y_dendro + theme(axis.text.x=element_blank(),
                                    axis.text.y=element_blank(), 
                                    plot.margin=unit(c(0,0,0,0),"mm"));
            allPlots[[length(allPlots)+1]] <- p_y_dendro;
        }
        
        allPlots[[length(allPlots)+1]] <- p_heatm;
        
        if (is.na(layout)) {
            if (dendrogram == "none") {
                layout <- matrix(1, nrow=1, ncol=1); # heatmap
            } else if (dendrogram == "col") {
                layout <- cbind(
                    matrix(c(rep(1, 182), rep(0, 18)), ncol=2, byrow=TRUE),
                    matrix(c(rep(2, 8*100)), ncol=8, byrow=TRUE));
            } else if (dendrogram == "row") {
                layout <- c(
                        rep(1, 2),
                        rep(2, 8)) # heatmap
                layout <- matrix(layout, nrow=10, ncol=1);
            } else if (dendrogram == "both") {
                layout <- rbind(
                    matrix(c(rep(0, 2*10), rep(2, 8*10)), ncol=10),
                    cbind(
                        matrix(c(rep(1, 182), rep(0, 18)), ncol=2, byrow=TRUE),
                        matrix(c(rep(3, 8*100)), ncol=8, byrow=TRUE)));
            } else {
                stop("bad dendrogram option provided.");
            }
        }
        
        multiplot(plotlist=allPlots, layout=layout);
        
        res <- list(plot=allPlots, plotLayout=layout, rowInd=rowInd,
            colInd=colInd, data=plotGenes[rowInd, colInd], rowDendrogram=hr.s, 
            colDendrogram=hc.s);
        
        return(invisible(res));
    }
)

#'@name genesBarplot
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesBarplot,MIGSAres-method
#'@exportMethod genesBarplot
#'
setGeneric(name="genesBarplot", def=function(migsaRes, ...) {
    standardGeneric("genesBarplot")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases genesBarplot,MIGSAres-method
#'
#'@importFrom ggplot2 aes geom_bar ggplot
#'@include MIGSAres-genesManipulation.R
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="genesBarplot",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, gsFilter=0) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=TRUE) >
                                    enrFilter, ];
        plotGenes <- genesInSets(actRes);
        
        # keep genes enriched in more than gsFilter gene sets
        plotGenes <- plotGenes[, colSums(plotGenes > 0) > gsFilter];
        
        # count 1 disregard if the gene enriched in more than 1 experiment
        plotGenes <- data.frame(id=colnames(plotGenes),
                                number=colSums(plotGenes > 0));
        
        ## todo: add symbol to plotGenes
        # plotGenes$symbol <- entrez2symbol(plotGenes$id);
        
        # reordering bars
        plotGenes$id <- factor(plotGenes$id,
                                levels=plotGenes$id[order(plotGenes$number,
                                                        decreasing=TRUE)]);
        plotGenes <- plotGenes[order(plotGenes$number, decreasing=TRUE),];
        
        p <- ggplot(plotGenes);
        p <- p + geom_bar(aes(id, y=number), stat="identity");
        print(p);
        
        return(p);
    }
)

#'@name migsaHeatmap
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases migsaHeatmap,MIGSAres-method
#'@exportMethod migsaHeatmap
#'
setGeneric(name="migsaHeatmap", def=function(migsaRes, ...) {
    standardGeneric("migsaHeatmap")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases migsaHeatmap,MIGSAres-method
#'
#'@importFrom ggdendro ggdendrogram
#'@importFrom ggplot2 aes element_blank element_text geom_text geom_tile ggplot
#'@importFrom ggplot2 guide_legend guides scale_fill_gradient scale_fill_manual
#'@importFrom ggplot2 scale_x_continuous scale_y_continuous scale_y_discrete 
#'@importFrom ggplot2 theme unit xlab ylab coord_flip
#'@importFrom stats as.dendrogram hclust
#'@importFrom vegan vegdist
#'@include MIGSAmultiplot.R
#'@include MIGSAres.R
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="migsaHeatmap",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0, expFilter=0, categories=list(), 
    categLabels=TRUE, colPal=c("white", "red"), col.dist=NA, row.dist=NA, 
    layout=NA, remove0Rows=TRUE, breaks=NA) {
        stopifnot(length(colPal) == 2);
        stopifnot(all(unlist(lapply(categories, length)) == ncol(migsaRes)-3));
        
        cOff <- enrCutoff(migsaRes);
        if (is.na(cOff)) {
            plotMigsaRes <- as.matrix(migsaRes[,-(1:3)]);
            rownames(plotMigsaRes) <- migsaRes$id;
            if (!is.na(breaks)) {
                limits <- seq(0, 1, length.out=breaks+1);
                limits[length(limits)] <- 1.1;
                limits <- cbind(limits[1:(length(limits)-1)],
                                limits[-1]);
                
                aux <- plotMigsaRes;
                aux[!is.na(aux)] <- unlist(lapply(aux, function(x) 
                    limits[ which(limits[,1] <= x & x < limits[,2]), 1 ]));
                plotMigsaRes <- aux;
            }
        } else {
            if (remove0Rows) {
                migsaRes <- migsaRes[rowSums(migsaRes[,-(1:3)], 
                    na.rm=TRUE) > 0,];
            }
            # terms filtering by enrFilter and expFilter
            keepRows <- which(rowSums(migsaRes[,-(1:3)], na.rm=TRUE) >= 
                enrFilter);
            keepCols <- which(colSums(migsaRes[,-(1:3)], na.rm=TRUE) >= 
                expFilter);
            migsaRes <- migsaRes[ keepRows, c(1:3, 3+keepCols) ];
            
            flog.debug(paste("In migsaHeatmap, after filtering, dim=",
                            dim(migsaRes)));

            plotMigsaRes <- as.matrix(migsaRes[,-(1:3)]);
            rownames(plotMigsaRes) <- migsaRes$id;
        }
        
        res <- migsaHeatmap(plotMigsaRes, categories=categories, 
            categLabels=categLabels, colPal=colPal, col.dist=col.dist, 
            row.dist=row.dist, layout=layout);
        return(invisible(res));
    }
)

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases migsaHeatmap,matrix-method
#'
#'@importFrom ggdendro ggdendrogram
#'@importFrom ggplot2 aes element_blank element_text geom_text geom_tile ggplot
#'@importFrom ggplot2 guide_legend guides scale_fill_gradient scale_fill_manual
#'@importFrom ggplot2 scale_x_continuous scale_y_continuous scale_y_discrete 
#'@importFrom ggplot2 theme unit xlab ylab coord_flip
#'@importFrom grDevices hcl
#'@importFrom stats as.dendrogram hclust
#'@importFrom vegan vegdist
#'@include MIGSAmultiplot.R
#'@include MIGSAres.R
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="migsaHeatmap",
    signature=c("matrix"),
    definition=function(migsaRes, categories=list(), categLabels=TRUE, 
    colPal=c("white", "red"), col.dist=NA, row.dist=NA, layout=NA) {
        stopifnot(length(colPal) == 2);
        stopifnot(all(unlist(lapply(categories, length)) == ncol(migsaRes)));
        stopifnot((!is.null(rownames(migsaRes))) & 
                            (!is.null(colnames(migsaRes))));
        
        plotInfo <- melt(migsaRes);
        colnames(plotInfo) <- c("Rows", "Cols", "Value");
        isNumeric <- is.numeric(migsaRes[1,1]);
        if (is.na(col.dist)) {
            col.dist <- ifelse(isNumeric, "euclidean", "jaccard");
        }
        if (is.na(row.dist)) {
            row.dist <- ifelse(isNumeric, "euclidean", "jaccard");
        }
        
        # jaccard clustering per col
        ddc.s <- suppressWarnings(vegdist(t(migsaRes), col.dist, na.rm=TRUE));
        ddc.s[is.na(ddc.s)] <- 0;
        hc.s <- hclust(ddc.s, method="average");
        colInd <- hc.s$order;
        
        # jaccard clustering per row
        ddr.s <- suppressWarnings(vegdist(migsaRes, row.dist, na.rm=TRUE));
        ddr.s[is.na(ddr.s)] <- 0;
        hr.s <- hclust(ddr.s, method="average");
        rowInd <- hr.s$order;
        
        # reordering the plot by distances
        plotInfo$Cols <- factor(plotInfo$Cols, 
            levels=colnames(migsaRes)[colInd]);
        plotInfo$Rows <- factor(plotInfo$Rows, 
            levels=rownames(migsaRes)[rowInd]);
        
        p_heatm <- ggplot(plotInfo, aes(Cols, Rows));
        p_heatm <- p_heatm + geom_tile(aes(fill=Value));
        p_heatm <- p_heatm + xlab(NULL) + ylab(NULL);
        p_heatm <- p_heatm + scale_y_discrete(breaks=NULL);
        
        if (isNumeric) {
            p_heatm <- p_heatm + scale_fill_gradient(low=colPal[[2]], 
                high=colPal[[1]]);
        } else {
            p_heatm <- p_heatm + scale_fill_manual(values=colPal)
        }
        
        p_heatm <- p_heatm + theme(legend.position="bottom", 
                            legend.key.size=unit(8, "points"), 
                            legend.text=element_text(size=unit(8, "points")), 
                            plot.margin=unit(c(0,0,0,0),"mm"),
                            axis.text.x=element_text(angle = 45, hjust=1));
        p_heatm <- p_heatm + guides(fill=guide_legend(nrow=2));
        
        # now the dendro
        p_dendro <- ggdendrogram(hc.s, labels=TRUE);
        p_dendro <- p_dendro + scale_y_continuous(expand=c(0.005,0));
        p_dendro <- p_dendro + scale_x_continuous(expand=
            c(1/(ncol(migsaRes)-3)/2,0));
        p_dendro <- p_dendro + xlab(NULL) + ylab(NULL)
        p_dendro <- p_dendro + theme(axis.text.x=element_blank(),
                                axis.text.y=element_blank(), 
                                plot.margin=unit(c(0,0,0,0),"mm"));
        
        # now each category plot
        categPlots <- lapply(categories, function(actCategory) {
            actCategory <- actCategory[colInd];
            actCategory <- data.frame(x=factor(1:length(actCategory)), 
                            y=factor(1), cat=actCategory);
            
            p_actCat <- ggplot(actCategory, aes(x, y));
            p_actCat <- p_actCat + geom_tile(aes(fill=cat));
            
            if (categLabels) {
                p_actCat <- p_actCat + geom_text(aes(label=cat));
            }
            p_actCat <- p_actCat + xlab(NULL) + ylab(NULL);
            p_actCat <- p_actCat + scale_y_discrete(breaks=NULL);
            p_actCat <- p_actCat + theme(legend.position="none", 
                                plot.margin=unit(c(0,0,0,0),"mm"), 
                                axis.title.x=element_blank(), 
                                axis.text.x=element_blank(), 
                                axis.ticks.x=element_blank());
            p_actCat <- p_actCat + guides(fill=guide_legend(nrow=1));
            
            n <- length(levels(actCategory$cat));
            actColors <- hcl(h=seq(15, 375, length=n+1), l=65, c=100)[1:n];
            p_actCat <- p_actCat + 
                scale_fill_manual(values=actColors, drop=FALSE)
            
            return(p_actCat);
        })
        
        allPlots <- list();
        allPlots[[1]] <- p_dendro;
        allPlots[seq_len(length(categPlots))+1] <- categPlots;
        allPlots[[length(allPlots)+1]] <- p_heatm;
        
        if (is.na(layout)) {
            layout <- c(
                    rep(1, 10),
                    rep(2:(length(categPlots)+1), 
                        each=floor(10 / length(categPlots))),
                    rep(length(allPlots), 80)) # heatmap
            layout <- c(rep(1, 100-length(layout)), layout);
            layout <- matrix(layout, nrow=100, ncol=1);
        }
        
        multiplot(plotlist=allPlots, layout=layout);
        
        res <- list(plot=allPlots, plotLayout=layout, rowInd=rowInd, 
            colInd=colInd, data=migsaRes[ rowInd, colInd ], 
            rowDendrogram=hr.s, colDendrogram=hc.s);
        return(invisible(res));
    }
)

#'@name geneSetBarplot
#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases geneSetBarplot,MIGSAres-method
#'@exportMethod geneSetBarplot
#'
setGeneric(name="geneSetBarplot", def=function(migsaRes, ...) {
    standardGeneric("geneSetBarplot")
})

#'@inheritParams MIGSAres-plots
#'@rdname MIGSAres-plots
#'@aliases geneSetBarplot,MIGSAres-method
#'
#'@importFrom ggplot2 aes geom_bar ggplot
#'@include MIGSAres-setEnrCutoff.R
#'
setMethod(
    f="geneSetBarplot",
    signature=c("MIGSAres"),
    definition=function(migsaRes, enrFilter=0) {
        stopifnot(validObject(migsaRes));
        
        # we must have a cutoff for this function
        migsaRes <- setDefaultEnrCutoff(migsaRes);
        
        # keep gene sets enriched in more than enrFilter experiments
        actRes <- migsaRes[ rowSums(migsaRes[,-(1:3)], na.rm=TRUE) >
                                enrFilter, ];
        plotRes <- data.frame(actRes[,1:3],
                                number=rowSums(actRes[,-(1:3)], na.rm=TRUE));
        
        # reordering bars
        plotRes$id <- factor(plotRes$id,
                            levels=plotRes$id[order(plotRes$number,
                                                    decreasing=TRUE)]);
        plotRes <- plotRes[order(plotRes$number, decreasing=TRUE),];
        
        p <- ggplot(plotRes);
        p <- p + geom_bar(aes(id, y=number), stat="identity");
        print(p);
        
        return(p);
    }
)
