% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IGSAinput-getterSetters.R
\docType{methods}
\name{IGSAinput-getterSetters}
\alias{IGSAinput-getterSetters}
\alias{name}
\alias{name,IGSAinput-method}
\alias{name<-}
\alias{name<-,IGSAinput-method}
\alias{fitOptions}
\alias{fitOptions,IGSAinput-method}
\alias{fitOptions<-}
\alias{fitOptions<-,IGSAinput-method}
\alias{exprData}
\alias{exprData,IGSAinput-method}
\alias{exprData<-}
\alias{exprData<-,IGSAinput-method}
\alias{geneSetsList}
\alias{geneSetsList,IGSAinput-method}
\alias{geneSetsList<-}
\alias{geneSetsList<-,IGSAinput-method}
\alias{gseaParams}
\alias{gseaParams,IGSAinput-method}
\alias{gseaParams<-}
\alias{gseaParams<-,IGSAinput-method}
\alias{seaParams}
\alias{seaParams,IGSAinput-method}
\alias{seaParams<-}
\alias{seaParams<-,IGSAinput-method}
\title{Accessors for IGSAinput class}
\usage{
"IGSAinput-getterSetters"(object)

name(object)

\S4method{name}{IGSAinput}(object)

name(object) <- value

\S4method{name}{IGSAinput}(object) <- value

fitOptions(object)

\S4method{fitOptions}{IGSAinput}(object)

fitOptions(object) <- value

\S4method{fitOptions}{IGSAinput}(object) <- value

exprData(object)

\S4method{exprData}{IGSAinput}(object)

exprData(object) <- value

\S4method{exprData}{IGSAinput}(object) <- value

geneSetsList(object)

\S4method{geneSetsList}{IGSAinput}(object)

geneSetsList(object) <- value

\S4method{geneSetsList}{IGSAinput}(object) <- value

gseaParams(object)

\S4method{gseaParams}{IGSAinput}(object)

gseaParams(object) <- value

\S4method{gseaParams}{IGSAinput}(object) <- value

seaParams(object)

\S4method{seaParams}{IGSAinput}(object)

seaParams(object) <- value

\S4method{seaParams}{IGSAinput}(object) <- value
}
\arguments{
\item{object}{IGSAinput object.}

\item{value}{value to replace in the slot.}
}
\value{
Modified IGSAinput object or desired slot.
}
\description{
Getters and setters functions to access IGSAinput object slots.
}
\examples{
## Lets create a basic IGSAinput object.
## First create a expression matrix.
maData <- matrix(rnorm(10000),ncol=4);
rownames(maData) <- 1:nrow(maData); # It must have rownames (gene names).
maExprData <- new("MAList",list(M=maData));

## Now lets create the FitOptions object.
myFOpts <- FitOptions(c("Cond1", "Cond1", "Cond2", "Cond2"));

## And now we can create our IGSAinput ready for MIGSA.
igsaInput <- IGSAinput(name="igsaInput", expr_data=maExprData, 
fit_options=myFOpts);

## Lets get igsaInput values, and modify its name.
name(igsaInput);
name(igsaInput) <- "newName";
fitOptions(igsaInput); 
exprData(igsaInput);
geneSetsList(igsaInput);
gseaParams(igsaInput);
seaParams(igsaInput);

}
\seealso{
\code{\link{IGSAinput-class}}
}
