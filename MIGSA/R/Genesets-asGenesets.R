#'Creates a GeneSetCollection from a list
#'
#'\code{as.Genesets} creates a GeneSetCollection object from the data present 
#'in a list. Each element will parse to a GeneSet. For each list element, its 
#'name will be the GeneSet setName, and the content are the genes.
#'
#'@param x list of character vectors which are the genes corresponding to each 
#'GeneSet. The list must have names (unique).
#'@param is_GO logical indicating if this gene sets are from the Gene Ontology.
#'If true, then each GeneSet setName must be a GO id.
#'@param ... not in use.
#'
#'@return A GeneSetCollection object.
#'
#'@docType methods
#'@name as.Genesets
#'@rdname Genesets-as.Genesets
#'@seealso \code{\link{Genesets-enrichr}}
#'@seealso \code{\link{geneSetsFromFile}}
#'@seealso \code{\link{loadGo}}
#'
#'@exportMethod as.Genesets
setGeneric("as.Genesets", def=function(x, ...) {
    standardGeneric("as.Genesets")
})

#'@rdname Genesets-as.Genesets
#'@inheritParams as.Genesets
#'@aliases as.Genesets,list-method
#'
#'@importFrom GSEABase GeneSet GeneSetCollection GOCollection NullCollection
#'@examples
#'## Lets create a list with three manually created gene sets and load it as a 
#'## GeneSetCollection object.
#'myGs1 <- as.character(1:10);
#'myGs2 <- as.character(15:21);
#'myGs3 <- as.character(25:30);
#'myGssList <- list(myGs1, myGs2, myGs3);
#'names(myGssList) <- c("myGs1", "myGs2", "myGs3");
#'myGss <- as.Genesets(myGssList);
#'
setMethod(
    f="as.Genesets",
    signature=c("list"),
    definition=function(x, is_GO=FALSE) {
        if (length(unique(names(x))) != length(x)) {
            stop("The list must have all unique names.");
        }
        
        collection <- NullCollection();
        if (is_GO)
            collection <- GOCollection();
        
        act_gss <- lapply(names(x), function(act_id) {
            # these are the genes
            act_genes <- unique(x[[act_id]]);
            
            # if there is any gene then create the Geneset
            if (length(act_genes) == 0 ||
                (length(act_genes) == 1 && act_genes == "")) {
                act_gs <- NULL;
            } else {
                act_gs <- GeneSet(act_genes, setName=act_id, 
                    collectionType=collection, setIdentifier="");
            }
            return(act_gs);
        })
        # delete invalid genesets
        act_gss <- act_gss[!unlist(lapply(act_gss, is.null))];
        
        if (length(act_gss) < 1) {
            stop("Could not create any Geneset object,
                    maybe they where empty");
        }
        res <- GeneSetCollection(act_gss);
        return(res)
    }
)
