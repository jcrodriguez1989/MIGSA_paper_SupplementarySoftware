#'Creates a GeneSetCollection object using the Gene Ontology data base
#'
#'\code{loadGo} creates a GeneSetCollection object from the data present at 
#'the Gene Ontology data base (org.Hs.eg.db R package).
#'
#'@param ontology character indicating which ontology must be loaded. 
#'Must be one of BP, MF or CC.
#'
#'@return A GeneSetCollection object  (Genes are with their EntrezGene ID).
#'
#'@docType methods
#'@name loadGo
#'@rdname Genesets-loadGo
#'@seealso \code{\link{as.Genesets}}
#'@seealso \code{\link{Genesets-enrichr}}
#'@seealso \code{\link{geneSetsFromFile}}
#'
#'@exportMethod loadGo
setGeneric(name="loadGo", def=function(ontology) {
    standardGeneric("loadGo")
})

#'@rdname Genesets-loadGo
#'@inheritParams loadGo
#'@aliases loadGo,character-method
#'
#'@importFrom AnnotationDbi as.list mappedkeys Ontology Term
#'@importFrom GSEABase GOCollection GeneSet EntrezIdentifier GeneSetCollection
#'@importFrom org.Hs.eg.db org.Hs.egGO2ALLEGS
#'@examples
#'## Lets load the Cellular Components gene sets from the Gene Ontology.
#'ccGSets <- loadGo("CC");
#'
setMethod(
    f="loadGo",
    signature=c("character"),
    definition=function(ontology) {
        ontology <- toupper(ontology);
        if (!ontology %in% c("BP", "CC", "MF")) {
            stop("Ontology must be one of: 'BP', 'CC', 'MF'.");
        }
        
        # download GO gene sets
        # org.Hs.egGO is an R object that provides mappings between entrez 
        # gene identifiers and the GO identifiers that they are directly 
        # associated with
        # go <- org.Hs.egGO;
        
        # org.Hs.egGO2ALLEGS is an R object that provides mappings between 
        # a given GO identifier and all of the Entrez Gene identifiers 
        # annotated at that GO term OR TO ONE OF ITS CHILD NODES
        go <- org.Hs.egGO2ALLEGS;
        goIds <- mappedkeys(go);
        
        # filter the desired ontology
        ontologies <- Ontology(goIds);
        go <- go[!is.na(ontologies) & ontologies == ontology];
        go <- as.list(go);
        go <- lapply(go, unique);
        
        goCollection <- GOCollection(ontology=ontology,
            evidenceCode=as.character(NA));
        
        terms <- Term(names(go));
        gsets <- lapply(seq_along(go), function(i) {
            gset <- GeneSet(go[[i]],
                setIdentifier=terms[[i]],
                setName=names(go)[[i]],
                geneIdType=EntrezIdentifier(), 
                collectionType=goCollection);
            
            return(gset);
        });
        
        gsetsColl <- GeneSetCollection(gsets);
        return(gsetsColl);
    }
)
