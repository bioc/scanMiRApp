#' runFullScan
#'
#' Runs a full miRNA scan on all protein-coding transcripts (or UTRs) of an
#' annotation.
#'
#' @param annotation A \code{\link{ScanMiRAnno}} object
#' @param mods An optional `KdModelList` (defaults to the one in `annotation`)
#' @param annoFilter An optional `AnnotationFilter` or `AnnotationFilterList` to
#' filter the set of transcripts to be extracted
#' @param extract Which parts of the transcripts to extract. For `UTRonly`
#' (default) only the 3' UTR regions are extracted, `withORF` additionally
#' extracts the coding regions, and `exons` extracts all exons
#' @param shadow The size of the ribosomal shadow at the UTR starts
#' @param cores The number of threads to use. Alternatively accepts a
#' \code{\link[BiocParallel]{BiocParallelParam-class}}, as for instance produced
#' by \code{\link[BiocParallel]{MulticoreParam}}.
#' @param maxLogKd The maximum log_kd of sites to report
#' @param save.path Optional, the path to which to save the results
#' @param ... Arguments passed to `scanMiR::findSeedMatches`
#'
#' @return A `GRanges` object
#'
#' @export
#' @import Biostrings scanMiR
#' @importFrom BiocParallel SerialParam MulticoreParam bpnworkers
#' @importFrom GenomicFeatures extractTranscriptSeqs threeUTRsByTranscript
#' exonsBy cdsBy
#' @importFrom GenomicRanges strand
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom stats setNames
#' @importFrom AnnotationFilter AnnotationFilterList SeqNameFilter
#' @examples
#' anno <- ScanMiRAnno("fake")
#' m <- runFullScan( annotation=anno )
#' m
runFullScan <- function(annotation, mods=NULL, annoFilter = NULL,
                        extract=c("UTRonly", "withORF", "exons"),
                        onlyCanonical=TRUE, shadow=15, cores=1,
                        maxLogKd=c(-1,-1.5), save.path=NULL, ...){

  stopifnot(is(annotation, "ScanMiRAnno"))
  if(is.null(mods)) mods <- annotation$models
  stopifnot(is(mods,"KdModelList"))

  if(is.numeric(cores) && length(cores)==1){
    cores <- as.integer(cores)
    if(cores>1){
      BP <- MulticoreParam(cores, progress=TRUE)
    }else{
      BP <- SerialParam()
    }
  }else if(is(cores,"BiocParallelParam")){
    BP <- cores
  }else{
    warning("`cores` argument of unknown type -- will be ignored.")
    BP <- SerialParam()
  }

  message("Loading annotation")
  genome <- annotation$genome
  ensdb <- annotation$ensdb
  if(is(ensdb,"EnsDb")) seqlevelsStyle(genome) <- "Ensembl"

  # restrict to canonical chromosomes
  canonical_chroms <- seqlevels(genome)[!grepl('_', seqlevels(genome))]
  filt <- SeqNameFilter(canonical_chroms)

  if(!is.null(annoFilter)){
    if(is(annoFilter, "AnnotationFilterList") ||
       is(annoFilter, "AnnotationFilter"))
      filt <- AnnotationFilterList(filt, annoFilter)
    else
      stop("filter must be either `AnnotationFilter` or `AnnotationFilterList`")
    if(!is(ensdb,"EnsDb"))
      warning("`annoFilter` is ignored when the annotation is not in ",
              "`EnsDb` format", immediate.=TRUE)
  }
  message("Extracting transcripts")
  seqs <- getTranscriptSequence(tx=NULL, annotation, annoFilter=annoFilter,
                                extract=extract)

  message("Scanning with ", bpnworkers(BP), " thread(s)")
  m <- findSeedMatches(seqs, mods, shadow=shadow, maxLogKd=maxLogKd,
                       onlyCanonical=onlyCanonical, BP=BP, verbose=FALSE, ...)

  md <- metadata(anno$ensdb)
  md <- setNames(md$value,md$name)
  if(!("genome_build" %in% names(md)))
    md[["genome_build"]] <- md[["Genome"]]
  if(!("ensembl_version" %in% names(md))) md[["ensembl_version"]] <- NA

  if(is(m, "GRanges")) {
    metadata(m)$genome_build <- md[["genome_build"]]
    metadata(m)$ensembl_version <- md[["ensembl_version"]]
  } else {
    attr(m, "ensembl_version") <- md[["ensembl_version"]]
    attr(m, "genome_build") <- md[["genome_build"]]
  }
  if(is.null(save.path)) return(m)
  saveRDS(m, file=save.path)
  rm(m)
  gc()
  message("Saved in: ", save.path)
}
