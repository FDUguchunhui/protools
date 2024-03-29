% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{SpC_List}
\alias{SpC_List}
\title{SpC list constructor}
\usage{
SpC_List(
  df,
  annotation,
  NA_substitution = NULL,
  proteins_filter = NULL,
  replicates_remove = NULL,
  replicates_keep = NULL
)
}
\arguments{
\item{df}{a data.frame with proteins accession number as row names, replicate (sample) identifier as column names, with
numeric data in each cell}

\item{annotation}{an annotation table provide extra information about the each replicate in the df variable,
such as disease type, subtype, and biology processing approach, etc. The rowname of annotation should be replicate identifiers that
is also used in the column names of df. Caution: the rownames of annotation are supposed to match with colnames of df.}

\item{NA_substitution}{how NA value in the df should be replace, can be NULL or a scalar value. If it is NULL
then all NA value will be kept as NA}

\item{proteins_filter}{A String vector, proteins that you want to keep in the data and should be compatible with rownames in df.}

\item{replicates_remove}{A String vector, replicates that you want to remove
from the data (corresponding replicates will be remove both from df and annotation),}
}
\value{
a list with two elements: matrix and annotation
\itemize{
\item matrix - the data matrix that is converted from the dataframe df using provided parameters
\item annotation - the original annotation passed into the function after removing rownames in replicates_remove
}
}
\description{
Create SpC_list object for conviently conducting downstream analysis
}
