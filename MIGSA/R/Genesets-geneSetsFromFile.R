#'Creates a GeneSetCollection object from a file
#'
#'\code{geneSetsFromFile} creates a GeneSetCollection object from the data 
#'present in a file. The file must be a tab separated values file (tsv). Each 
#'line will parse to a GeneSet. First field will be the GeneSet setName, the 
#'second the setIdentifier and the remaining are the genes.
#'
#'@param filePath character with the path of the file to parse.
#'@param is_GO logical indicating if this gene sets are from the Gene Ontology.
#'If true, then each gene GeneSet setName must be a GO id.
#'@param ... not in use.
#'
#'@return A GeneSetCollection object.
#'
#'@docType methods
#'@name geneSetsFromFile
#'@rdname Genesets-geneSetsFromFile
#'@seealso \code{\link{as.Genesets}}
#'@seealso \code{\link{Genesets-enrichr}}
#'@seealso \code{\link{loadGo}}
#'
#'@exportMethod geneSetsFromFile
#'@examples
#'## Create some fake gene sets in a data.frame to save them in disk and then
#'## load them (10 gene sets with 20 genes each (it is not neccesary that they
#'## have the same number of genes).
#'gsets <- data.frame(
#'     IDs=paste("set", 1:10), 
#'     Names=rep("", 10), 
#'     matrix(paste("gene", 1:(10*20)), nrow=10));
#'
#'## And save this file as a tab separated file.
#'geneSetsFile <- paste(tempdir(), "/fakeGsets.tsv", sep="");
#'write.table(gsets, file=geneSetsFile, sep="\t",
#'     col.names=FALSE, row.names=FALSE, quote=FALSE);
#'
#'## Now lets load this tsv file as a GeneSetCollection object.
#'myGsets <- geneSetsFromFile(geneSetsFile);
#'
#'## And lets delete this tsv file (so we dont have garbage in our disk).
#'unlink(geneSetsFile);
#'
setGeneric(name="geneSetsFromFile", def=function(filePath, ... ) {
    standardGeneric("geneSetsFromFile")
})

#'@rdname Genesets-geneSetsFromFile
#'@inheritParams geneSetsFromFile
#'@aliases geneSetsFromFile,character-method
#'
#'@importFrom futile.logger flog.error
setMethod(
    f="geneSetsFromFile",
    signature=c("character"),
    definition=function(filePath, is_GO=FALSE) {
        stopifnot(length(filePath) == 1);
        
        # check if file exists, its not a directory, its readable
        if (!file.exists(filePath) || dir.exists(filePath)
                            || (file.access(filePath, 4) == -1)) {
            flog.error("Gene sets file path error, check if it is readable.");
        }
        
        collection <- NullCollection();
        if (is_GO)
            collection <- GOCollection();
        
        tmp <- readLines(filePath);
        gsets <- lapply(tmp, function(actLine) {
            # for each line, fst item is id, snd is name, rest are genes
            t <- strsplit(actLine,'\t')[[1]]
            
            actId <- t[[1]];
            actName <- t[[2]];
            if (actName == "")
                actName <- actId;
            
            actGS <- GeneSet(unique(t[3:length(t)]), setIdentifier=actName,
                    setName=actId, collectionType=collection);
            return(actGS);
        })
        
        return(GeneSetCollection(gsets));
    }
)
