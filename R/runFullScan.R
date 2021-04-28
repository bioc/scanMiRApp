#' runFullScan
#'
#' Runs a full miRNA scan on all protein-coding transcripts (or UTRs) of an
#' annotation.
#' 
#' @param annotation A \code{\link{ScanMiRAnno}} object
#' @param mods An optional `KdModelList` (defaults to the one in `annotation`)
#' @param UTRonly Whether to scan only UTRs
#' @param shadow The size of the ribosomal shadow at the UTR starts
#' @param cores The number of threads to use
#' @param maxLogKd The maximum log_kd of sites to report
#' @param save.path Optional, the path to which to save the results
#' @param ... Arguments passed to `scanMiR::findSeedMatches`
#'
#' @export
#' @importFrom BiocParallel SerialParam MulticoreParam
#' @importFrom GenomicFeatures extractTranscriptSeqs threeUTRsByTranscript cdsBy
#' @import Biostrings scanMiR
#' @importFrom S4Vectors metadata metadata<-
#' @examples 
#' # not run
#' # anno <- ScanMiRAnno("Rnor_6")
#' # seq <- runFullScan( annotation=anno )
runFullScan <- function(annotation, mods=NULL, UTRonly=TRUE, onlyCanonical=TRUE,
                        shadow=15, cores=1, maxLogKd=c(-1,-1.5), 
                        save.path=NULL, ...){
  message("Loading annotation")
  stopifnot(is(annotation, "ScanMiRAnno"))
  if(is.null(mods)) mods <- annotation$models
  stopifnot(is(mods,"KdModelList"))
  genome <- annotation$genome
  ensdb <- annotation$ensdb
  seqlevelsStyle(genome) <- "Ensembl"
  
  # restrict to canonical chromosomes
  canonical_chroms <- seqlevels(genome)[!grepl('_', seqlevels(genome))]
  filt <- SeqNameFilter(canonical_chroms)
  
  message("Extracting transcripts")
  grl_UTR <- suppressWarnings(threeUTRsByTranscript(ensdb, filter=filt))
  seqs <- extractTranscriptSeqs(genome, grl_UTR)
  utr.len <- lengths(seqs)
  names(utr.len) <- names(seqs)
  if(!UTRonly){
    grl_ORF <- cdsBy(ensdb, by="tx", filter=filt)
    seqs_ORF <- extractTranscriptSeqs(genome, grl_ORF)
    tx_info <- data.frame(strand=unlist(unique(strand(grl_ORF))))
    orf.len <- lengths(seqs_ORF)
    names(orf.len) <- names(grl_ORF)
    tx_info$ORF.length <- orf.len[row.names(tx_info)]
    seqs_ORF[names(seqs)] <- xscat(seqs_ORF[names(seqs)],seqs)
    seqs <- seqs_ORF
    rm(seqs_ORF)
    mcols(seqs)$ORF.length <- orf.len[names(seqs)]
  }else{
    tx_info <- data.frame(strand=unlist(unique(strand(grl_UTR))))
  }
  tx_info$UTR.length <- utr.len[row.names(tx_info)]
  
  message("Scanning with ", cores, " core(s)")
  if(cores>1){
    BP <- MulticoreParam(cores, progress=TRUE)
  }else{
    BP <- SerialParam()
  }
  m <- findSeedMatches(seqs, mods, shadow=shadow, maxLogKd=maxLogKd, 
                       onlyCanonical=onlyCanonical, BP=BP, ...)
  
  md <- metadata(anno$ensdb)
  md <- setNames(md$value,md$name)

  if(is(m, "GRanges")) {
    metadata(m)$tx_info <- tx_info
    metadata(m)$genome_build <- md[["genome_build"]]
    metadata(m)$ensembl_version <- md[["ensembl_version"]]
  } else {
    attr(m, "tx_info") <- tx_info
    attr(m, "ensembl_version") <- md[["ensembl_version"]]
    attr(m, "genome_build") <- md[["genome_build"]]
  }
  if(is.null(save.path)) return(m)
  saveRDS(m, file=save.path)
  rm(m)
  gc()
  message("Saved in: ", save.path)
}
