% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/truncate_sequences.length.R
\name{truncate_sequences.length}
\alias{truncate_sequences.length}
\title{Truncate DNA Sequences to Specified Length}
\usage{
truncate_sequences.length(sequences, length, quality_scores)
}
\arguments{
\item{sequences}{A character vector of DNA sequences to truncate.}

\item{length}{Numeric. The length to truncate DNA sequences to.}

\item{quality_scores}{An optional character vector of DNA sequence quality scores. If supplied, these will be truncated to their corresponding truncated DNA sequences.}
}
\value{
If quality scores are not provided, then a character vector of truncated DNA sequences is returned. If quality scores are provided, then a list containing two elements is returned. The first element is a character vector of truncated DNA sequences, and the second element is a character vector of quality scores which have been truncated to their corresponding truncated DNA sequences.
}
\description{
Truncates DNA sequences to a specified length.
}
\examples{
truncate_sequences.length(sequences=c("ATATAGCGCG","TGCCGATATA","ATCTATCACCGC"),
                          length=5,
                          quality_scores=c("989!.C;F@\"","A((#-#;,2F","HD8I/+67=1>?"))
}
\seealso{
\code{\link[=truncate_sequences.quality_score]{truncate_sequences.quality_score}} for truncating DNA sequences by Phred quality score. \cr \cr
\code{\link[=truncate_sequences.probability]{truncate_sequences.probability}} for truncating DNA sequences by cumulative probability that all bases were called correctly. \cr \cr
\code{\link[=truncate_and_merge_pairs]{truncate_and_merge_pairs}} for truncating and merging forward and reverse DNA sequence reads.
}
