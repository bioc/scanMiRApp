% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scanMiRAnno.R
\docType{class}
\name{ScanMiRAnno-class}
\alias{ScanMiRAnno-class}
\alias{ScanMiRAnno}
\title{ScanMiRAnno}
\usage{
ScanMiRAnno(
  species = NULL,
  genome = NULL,
  ensdb = NULL,
  models = NULL,
  scan = NULL,
  aggregated = NULL,
  version = NULL,
  addDBs = list(),
  ...
)
}
\arguments{
\item{species}{The species/build acronym for automatic construction; if
omitted, `genome` and `ensdb` should be given. Current possible values are:
GRCh38, GRCm38, GRCm39, Rnor_6.}

\item{genome}{A \link[BSgenome]{BSgenome-class}, or a
\code{\link[rtracklayer]{TwoBitFile}}}

\item{ensdb}{An \link[ensembldb]{EnsDb-class} (or a
\link[GenomicFeatures]{TxDb-class}) object}

\item{models}{An optional KdModelList}

\item{scan}{An optional full scan (IndexedFst or GRanges)}

\item{aggregated}{An optional per-transcript aggregation (IndexedFst or
data.frame)}

\item{version}{optional ensembl version}

\item{addDBs}{A named list of additional tx-miRNA databases, each of which
should be a data.frame with the columns 'transcript', 'miRNA', and 'score'.}

\item{...}{Arguments passed to `AnnotationHub`}
}
\value{
A `ScanMiRAnno` object
}
\description{
ScanMiRAnno
}
\examples{
anno <- ScanMiRAnno(species="fake")
anno
}
