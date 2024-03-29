% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{nsaf}
\alias{nsaf}
\title{Calculate Normalized Spectral Abundance Factor (NASF)}
\usage{
nsaf(x, protein_length, per_count, na_fill = NULL)
}
\arguments{
\item{x}{a matrix of raw counts}

\item{protein_length}{a data.frame that contains information about protein length for
each protein represented by protein id in each row}
}
\value{
the NASF normalized matrix in the same shape as x
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
\deqn{
NSAF_j = \frac{SpC_j/proteinLength_j}{denominator} \\
denomiator = \sum{over all protein i} SpC_i/proteinLength_i
}
}
\references{
https://github.com/moldach/proteomics-spectralCount-normalization/blob/master/nsaf.R

McIlwain S, Mathews M, Bereman MS, Rubel EW, MacCoss MJ, Noble WS.
Estimating relative abundances of proteins from shotgun proteomics data.
BMC Bioinformatics. 2012 Nov 19;13:308. doi: 10.1186/1471-2105-13-308. PMID: 23164367; PMCID: PMC3599300.
}
