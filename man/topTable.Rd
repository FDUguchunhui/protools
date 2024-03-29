% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statistical_testing.R
\name{topTable}
\alias{topTable}
\title{Table of the Top Differentially Expressed tags}
\usage{
topTable(res, logFC, PAjusted)
}
\arguments{
\item{res}{summary table protein expression analysis including logFC, PValue, and PAdjusted}

\item{fold_change}{absolute log2 fold change threshold}

\item{PAdjusted}{filter tags that less than the adjusted p-value used}
}
\description{
filter the top differentially expressed tags based on adjustment p-value
and fold-change
}
