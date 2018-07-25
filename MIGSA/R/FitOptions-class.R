#'FitOptions S4 class implementation in R
#' 
#'This S4 class contains the parameters to provide for model fitting.
#'If the vector of samples is provided (must be two different, e.g. 
#'c("C1", "C1", "C2")) then it will contrast C1 vs. C2. If not, it should be 
#'provided with a data.frame x, the formula and the contrast, it will create 
#'the model matrix using x as data, and the formula.
#'
#'@param x There are two options for x:
#'\itemize{
#'\item It can be a character vector containing the two conditions (length 
#'must be the same as the number of subjects to use).
#'\item It can be a data.frame used as data by 
#'\code{\link[stats]{model.matrix}}.
#'}
#'@param formula (only used if x is data.frame) used by 
#'\code{\link[stats]{model.matrix}}.
#'@param contrast (only used if x is data.frame) the contrast to test.
#'@param ... not in use.
#'
#'@return FitOptions object.
#'
#'@docType methods
#'@name FitOptions-class
#'@rdname FitOptions-class
#'
#'@examples
#'## Supose we have 15 subjects, the first 8 from Condition1 and the last 7 
#'## from Condition2, lets create the corresponding FitOptions object to test
#'## Condition1 vs. Condition2.
#'l <- c(rep("Condition1", 8), rep("Condition2", 7));
#'fit_options <- FitOptions(l);
#'
#'## Otherwise if we have the data and formula for model.matrix function and 
#'## the desired contrast, we can create the FitOptions object as:
#'myData <- data.frame(cond=c(rep("Condition1", 8), rep("Condition2", 7)));
#'myFormula <- ~cond - 1;
#'myContrast <- c(-1, 1);
#'fit_options <- FitOptions(myData, myFormula, myContrast);
#'
#'@importFrom futile.logger flog.error
#'@importFrom stats model.matrix
#'@exportClass FitOptions
#'
setClass(
    Class="FitOptions",
    slots=c(
        col_data="data.frame",
        formula="formula",
        contrast="numeric",
        design_matrix="matrix"
    ),
    prototype=list(
    ),
    validity=function(object) {
        design_matrix <- object@design_matrix;
        formula <- object@formula;
        col_data <- object@col_data;
        contrast <- object@contrast;
        
#         design_matrix_error <- all.equal(design_matrix, model.matrix(formula,
#                                             data=col_data));
#         contst_names <- make.names(
#                     colnames(model.matrix(formula, data=col_data)));
        contrast_ok <- length(contrast) == ncol(design_matrix);
        if (!contrast_ok) {
            flog.error("Contrast length must be equal to design_matrix
                                                        columns number.");
        }
        
        return(contrast_ok);
    }
)

#'@inheritParams FitOptions-class
#'@rdname FitOptions-class
#'@export FitOptions
#'
FitOptions <- function(x, ...) {
    UseMethod("FitOptions", x);
}

#'@inheritParams FitOptions-class
#'@rdname FitOptions-class
#'@aliases FitOptions.default
#'@export FitOptions.default
#'@importFrom stats model.matrix
#'
FitOptions.default <- function(x, ...) {
    # checking that exactly two conditions are present
    if (length(x) < 2) {
        stop("More than two labels required.");
    }
    if (length(unique(x)) != 2) {
        stop("Exactly two possible conditions required.");
    }
    
    # creating the model from the two conditions
    act_col_data <- data.frame(cond=factor(x));
    act_formula  <- ~cond-1;
    
    act_contrast <- c(-1, 1);
    .Object <- FitOptions(x=act_col_data, formula=act_formula,
                            contrast=act_contrast);
    return(.Object);
}

#'@inheritParams FitOptions-class
#'@rdname FitOptions-class
#'@aliases FitOptions.data.frame
#'@export FitOptions.data.frame
#'@importFrom stats model.matrix
#'
FitOptions.data.frame <- function(x, formula, contrast, ...) {
    stopifnot(is(x, "data.frame"));
    stopifnot(is(formula, "formula"));
    stopifnot(is(contrast, "numeric"));
    
    # completing the model
    act_design <- model.matrix(formula, data=x);
    .Object <- new("FitOptions", col_data=x, formula=formula, 
                    contrast=contrast, design_matrix=act_design);
    return(.Object);
}
