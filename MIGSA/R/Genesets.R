#'@importClassesFrom GSEABase GeneSetCollection
#'@importFrom GSEABase setName
setGeneric(name="asList", def=function(object) {
    standardGeneric("asList")
})

setMethod(f="asList",
    signature=c("GeneSetCollection"),
    definition=function(object) {
        # translate a Genesets object to a list. Why? mGSZ uses lists.
        to <- lapply(object, geneIds);
        names(to) <- lapply(object, setName);
        
        return(to);
    }
)
