% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{build_10x_atac}
\alias{build_10x_atac}
\title{Read in multiple 10X ATAC samples with Signac and merge with common peak set}
\usage{
build_10x_atac(manifest = NULL, project = NULL)
}
\arguments{
\item{manifest}{The file manifest from build_manifest or a list of named lists where each entry is a sample with peaks.bed, singlecell.csv and fragments}

\item{project}{Name to give merged sample set}
}
\description{
Takes file manifest from build_manifest
following preffered merge procedure in Signac: https://satijalab.org/signac/articles/merging.html#merge-fragment-files-1
}
\details{
"Fail to connect" or "Network is unreachable" errors are likely the result of a curl command reaching for the annotation
but unable to connect to the server, possibly becuase of firewall issues. Try running the annotation block on a laptop
outside of a firewall and then import the result back into the main notebook.
}
\examples{
\dontrun{
srt = build_10X_atac(manifest = m, project = "A6440")
}

}
