#'List and download gene sets from enrichr database
#'
#'\code{enrichrGeneSets} lists the database names present at enrichr.
#'\code{downloadEnrichrGeneSets} creates a list of GeneSetCollection 
#'objects downloading the specified ones from enrichr website
#' (http://amp.pharm.mssm.edu /Enrichr/).
#'
#'@param pattern character indicating a pattern to filter the database names.
#'@param geneSetNames list of characters with the names of the gene sets to 
#'download. Must be listed at \code{\link{enrichrGeneSets}}.
#'@param deleteMultipleEntrez logical indicating if multiple Entrez IDs should 
#'be deleted or repeated.
#'Note: Enrichr uses Gene Symbol, if org.Hs.eg translates it into several 
#'Entrez IDs then if deleteMultipleEntrez == FALSE all Entrez are removed else 
#'all are included.
#'@param ... not in use.
#'
#'@return enrichrGeneSets: character with present database names.
#'downloadEnrichrGeneSets: list of GeneSetCollection objects (Genes are with 
#'their EntrezGene ID).
#'
#'@docType methods
#'@name Genesets-enrichr
#'@rdname Genesets-enrichr
#'@seealso \code{\link{as.Genesets}}
#'@seealso \code{\link{geneSetsFromFile}}
#'@seealso \code{\link{loadGo}}
#'
#'@examples
#'## Lets list all the gene sets that can be downloaded from Enichr website.
#'enrichrGeneSets();
#'
#'\dontrun{
#'## Now lets list only the gene sets that have BioCarta in their names 
#'## (different BioCarta versions).
#'enrichrGeneSets("BioCarta");
#'}
#'
#'## And lets download the latest BioCarta gene sets from Enrichr.
#'## Make sure you use the same names as listed with enrichrGeneSets() .
#'biocartaGSs <- downloadEnrichrGeneSets(c("BioCarta_2015"));
#'
setGeneric(name="Genesets-enrichr", def=function(pattern) {
    standardGeneric("Genesets-enrichr")
})

#'@name enrichrGeneSets
#'@inheritParams Genesets-enrichr
#'@rdname Genesets-enrichr
#'@aliases enrichrGeneSets,character-method
#'@exportMethod enrichrGeneSets
#'
setGeneric(name="enrichrGeneSets", def=function(pattern=".*") {
    standardGeneric("enrichrGeneSets")
})

#'@inheritParams Genesets-enrichr
#'@rdname Genesets-enrichr
#'@aliases enrichrGeneSets,character-method
#'
#'@importFrom Biobase testBioCConnection
#'@importFrom RJSONIO fromJSON
#'
setMethod(
    f="enrichrGeneSets",
    signature=character(),
    definition=function(pattern=".*") {
        if (!testBioCConnection()) stop("You must have internet connection.");
        
        # enrichr url
        datasetStatisticsUrl <-
                    "http://amp.pharm.mssm.edu/Enrichr/datasetStatistics";
        
        # donwload genesets list (basic info)
        datasetStatistics <- do.call(rbind,
                        fromJSON(datasetStatisticsUrl)$statistics);
        datasetStatistics <- apply(datasetStatistics, 2, unlist);
        datasetStatistics <- data.frame(datasetStatistics);
        
        # giving some format
        datasetStatistics[,2:4] <- apply(datasetStatistics[,2:4], 2,
                                                                as.numeric);
        datasetStatistics[,c(1,5)] <- apply(datasetStatistics[,c(1,5)], 2,
                                                                as.character);
        
        # filtering with pattern (by default returns all)
        datasetStatistics <- datasetStatistics[
            grep(pattern, datasetStatistics[,1], ignore.case=!FALSE), ];
        
        # order by gene set name
        datasetStatistics <- datasetStatistics[order(datasetStatistics[,1]),];
        
        return(datasetStatistics);
    }
)

#'@name downloadEnrichrGeneSets
#'@inheritParams Genesets-enrichr
#'@rdname Genesets-enrichr
#'@aliases downloadEnrichrGeneSets,character-method
#'@exportMethod downloadEnrichrGeneSets
#'
setGeneric(name="downloadEnrichrGeneSets", def=function(geneSetNames, ...) {
    standardGeneric("downloadEnrichrGeneSets")
})

#'@inheritParams Genesets-enrichr
#'@rdname Genesets-enrichr
#'@aliases downloadEnrichrGeneSets,ANY-method
#'
#'@importFrom AnnotationDbi as.list
#'@importFrom Biobase testBioCConnection
#'@importFrom futile.logger flog.info
#'@importFrom GSEABase GeneSet GeneSetCollection GOCollection NullCollection
#'@importFrom org.Hs.eg.db org.Hs.egALIAS2EG
#'
setMethod(
    f="downloadEnrichrGeneSets",
    signature=c("character"),
    definition=function(geneSetNames,
    deleteMultipleEntrez=!FALSE) {
        # enrichr url
        downloadUrlFst <-
"http://amp.pharm.mssm.edu/Enrichr/geneSetLibrary?mode=text&libraryName=";
        
        # filtering invalid gene sets
        if (length(geneSetNames) < 1) {
            stop("No gene sets found to download from Enrichr");
        }
        
        if (!testBioCConnection()) stop("You must have internet connection.");
        
        # to translate Symbols to Entrez
        symbol2entrez <- as.list(org.Hs.egALIAS2EG);
        
        libraries <- lapply(geneSetNames, function(libName) {
            # if name starts with GO then is_GO true
            is_GO <- grepl("^GO_", libName);
            collection <- NullCollection();
            if (is_GO)
                collection <- GOCollection();
            
            actUrl <- paste(downloadUrlFst, libName, sep="");
            # downloading gene set
            tmp <- try({ readLines(actUrl); });
            if (inherits(tmp, 'try-error')) return(NA);
            
            flog.info(paste("Downloaded", libName));
            tmp <- strsplit(tmp, "\t");
            flog.info("Converting Symbol to Entrez");
            geneSets <- lapply(tmp, function(actGS) {
                actId <- actGS[[1]];
                actName <- actGS[[2]];
                entrez <- lapply(actGS[2:length(actGS)], function(symbol) {
                    # trying to convert from Symbol to Entrez
                    translates <- symbol2entrez[[symbol]];
                    
                    if (is.null(translates)) {
                        # some times Enrichr returns gene symbols followed by a 
                        # ',1.0' , this is some kind of fixing it.
                        translates <- symbol2entrez[[gsub(",1.0", "", symbol)]];
                    }
                    
                    # if multiple Entrez returned then delete them or not
                    if (length(translates) > 1 && deleteMultipleEntrez) {
                        translates <- NULL;
                    }
                    
                    return(translates);
                })
                entrez <- Reduce(union, entrez);
                
                res <- NA;
                
                # if returned any gene, then create the Geneset
                if (length(entrez) > 0) {
                    res <- GeneSet(entrez, setIdentifier=actName,
                        setName=actId, geneIdType=EntrezIdentifier(),
                        collectionType=collection);
                }
                
                return(res);
            })
            
            # filter invalid gene sets
            geneSets <- geneSets[!is.na(geneSets)];
            res <- NA;
            
            if (length(geneSets) > 0) {
                res <- GeneSetCollection(geneSets);
            }
            
            return(res);
        })
        names(libraries) <- geneSetNames;
        libraries <- libraries[!is.na(libraries)];
        
        return(libraries);
    }
)
