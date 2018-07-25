#'MIGSAres exploratory functions
#'
#'Several R base overwritten functions to manipulate a MIGSAres object as a 
#'data.frame way.
#'NOTE: When subsetting a MIGSAres object, if it does not have the id, GS_Name 
#'and (at least) one experiment columns, then it wont be a MIGSAres object, 
#'i.e., migsaRes[,c("id","igsaInput1")] is no longer a MIGSAres object.
#'
#'@param x MIGSAres object.
#'@param y MIGSAres object.
#'@param object MIGSAres object.
#'@param name as used in \code{$}.
#'@param n as used in \code{\link[utils]{head}} and \code{tail}.
#'@param i as used in \code{[}.
#'@param j as used in \code{[}.
#'@param drop as used in \code{[} (default: FALSE).
#'
#'@return Desired object.
#'
#'@docType methods
#'@name MIGSAres-common
#'@rdname MIGSAres-common
#'
#'@include MIGSAres-class.R
#'@aliases dim,MIGSAres-method
#'@exportMethod dim
#'@examples
#'data(migsaRes);
#'## As we ran MIGSA for two experiments and 200 gene sets, it must have 200 
#'## rows, and five columns (id, Name, GS_Name, and the experiments names).
#'dim(migsaRes);
#'
#'## migsaRes shown as data.frame has these column names: id, Name, GS_Name, 
#'## and the experiments names. As we ran two experiments, names igsaInput1 
#'## and igsaInput2, we can use $ in these ways:
#'head(migsaRes$id);
#'table(migsaRes$Name);
#'table(migsaRes$GS_Name);
#'head(migsaRes$igsaInput1);
#'head(migsaRes$igsaInput2);
#'
#'colnames(migsaRes);
#'
#'head(migsaRes);
#'
#'## Or see the first 10
#'head(migsaRes, n=10);
#'
#'tail(migsaRes);
#'
#'## Or see the last 10
#'tail(migsaRes, n=10);
#'
#'## migsaRes shown as data.frame has these column names: id, Name, GS_Name, 
#'## and the experiments names. As we ran two experiments, names igsaInput1 
#'## and igsaInput2, we can use [ in these ways:
#'
#'## Lets get the first 5 rows and 4 columns (the result is a MIGSAres object).
#'migsaRes[1:5, 1:4];
#'class(migsaRes[1:5, 1:4]);
#'
#'## Lets get the experiments results. Note that this is not any more a 
#'## MIGSAres object.
#'migsaRes[, c("igsaInput1", "igsaInput2")];
#'class(migsaRes[, c("igsaInput1", "igsaInput2")]);
#'
#'migsaRes;
#'
#'migsaResDFrame <- as.data.frame(migsaRes);
#'head(migsaResDFrame);
#'
#'migsaRes1 <- migsaRes[,1:4];
#'migsaRes2 <- migsaRes[,c(1:3,5)];
#'migsaResMerged <- merge(migsaRes1, migsaRes2);
#'
setMethod("dim",
    signature=c("MIGSAres"),
    function(x) {
        stopifnot(validObject(x));
        return(dim(x@migsa_res_summary));
    }
)

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases $,MIGSAres-method
#'@include MIGSAres.R
#'@exportMethod $
#'
setMethod("$",
    signature=c("MIGSAres"),
    function(x, name) {
        stopifnot(validObject(x));
        res <- get_summary(x);
        
        return(res[,name]);
    }
)

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases colnames,MIGSAres
#'@exportMethod colnames
#'
setMethod("colnames",
    signature=c("MIGSAres"),
    function(x) {
        stopifnot(validObject(x));
        return(colnames(x@migsa_res_summary));
    }
)

## todo: define colnames(migsaRes) <-

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases head,MIGSAres
#'@exportMethod head
#'
setMethod("head",
    signature=c("MIGSAres"),
    function(x, n=6L) {
        stopifnot(validObject(x));
        n <- min(nrow(x@migsa_res_summary), n);
        return(x[1:n,, drop=FALSE]);
    }
)

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases tail,MIGSAres
#'@exportMethod tail
#'
setMethod("tail",
    signature=c("MIGSAres"),
    function(x, n=6L) {
        stopifnot(validObject(x));
        n <- min(nrow(x@migsa_res_summary), n);
        return(x[(nrow(x)-n+1):nrow(x),, drop=FALSE]);
    }
)

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases [,MIGSAres,ANY,ANY,ANY-method
#'@include MIGSAres.R
#'@exportMethod "["
#'
setMethod("[",
    signature=c("MIGSAres"),
    function(x, i, j, drop=FALSE) {
        stopifnot(validObject(x));
        # subset from the results summary
        data_frame <- x@migsa_res_summary[i, j, drop=drop];
        
        if (is.null(dim(data_frame)) || ncol(data_frame) == 0 || 
            nrow(data_frame) == 0) {
            # if the structure has broken then return a non MIGSAres object
            res <- get_summary(x)[i, j, drop=drop];
            return(res);
        }
        
        # present experiments in the subsetted object
        act_methods <- setdiff(colnames(data_frame),
                                    c("id", "GS_Name", "Name"));
        data_frame_all <- x@migsa_res_all;
        
        # if id and GS_Name are removed then object is not valid, so we return
        # just the data.frame
        if (!all(c("id", "GS_Name") %in% colnames(data_frame)) | 
                length(act_methods) == 0) {
            res <- get_summary(x)[i, j, drop=drop];
            return(res);
        }
        
        # filter the removed experiments from the all results data frame
        data_frame_all <- data_frame_all[
            data_frame_all$experiment_name %in%  act_methods &
            data_frame_all$id %in% data_frame$id &
            data_frame_all$gene_set_name %in% data_frame$GS_Name , ];
        
        x@migsa_res_all <- data_frame_all;
        x@migsa_res_summary <- data_frame;
        
        return(x);
    }
)

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases show,MIGSAres
#'@include MIGSAres.R
#'@exportMethod show
#'
setMethod("show",
    signature=c("MIGSAres"),
    function(object) {
        stopifnot(validObject(object));
        res <- get_summary(object);
        
        print(res);
    }
)

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases as.data.frame,MIGSAres
#'@include MIGSAres.R
#'@exportMethod as.data.frame
#'
setMethod(
    f="as.data.frame",
    signature=c("MIGSAres"),
    definition=function(x) {
        res <- get_summary(x);
        return(res)
    }
)

#'@inheritParams MIGSAres exploratory functions
#'@rdname MIGSAres-common
#'@aliases merge,MIGSAres,MIGSAres
#'@include MIGSAres.R
#'@exportMethod merge
#'
setMethod(
    f="merge",
    signature=c("MIGSAres", "MIGSAres"),
    definition=function(x, y) {
        x_all <- x@migsa_res_all;
        y_all <- y@migsa_res_all;
        
        merged_all <- rbind(x_all, y_all);
        
        # checking that we dont have experiment names repeated between x and y
        if (nrow(unique(merged_all[,list(experiment_name,gene_set_name, id)])) 
            != nrow(merged_all[,list(experiment_name,gene_set_name, id)])) {
            stop("Merging: Experiment names can not be the same");
        }
        newGenesRank <- unique(c(x@genes_rank, y@genes_rank));
        
        .Object <- MIGSAres(merged_all, 
            genes_rank=newGenesRank);
        return(.Object);
    }
)
