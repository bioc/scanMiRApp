#' runFullScan
#' 
#' @param species 
#' @param mods 
#' @param UTRonly 
#' @param shadow 
#' @param cores 
#' @param maxLogKd 
#' @param save.path 
#' @param ... 
#'
#' @export
runFullScan <- function(species, mods=NULL, UTRonly=TRUE, shadow=15, cores=8, maxLogKd=c(-0.3,-0.3), save.path=NULL, ...){
  message("Loading annotation")
  suppressPackageStartupMessages({
    library(ensembldb)
    library(AnnotationHub)
    library(BSgenome)
    library(BiocParallel)
  })
  ah <- AnnotationHub()
  species <- match.arg(species, c("mmu","hsa","rno"))
  if(species=="hsa"){
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    if(is.null(mods)) mods <- readRDS(file = "/mnt/schratt/miRNA_KD/Data_Output/mods_hsa_comp.rds")
    ahid <- rev(query(ah, c("EnsDb", "Homo sapiens"))$ah_id)[1]
  }else if(species=="mmu"){
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(is.null(mods)) mods <- readRDS(file = "/mnt/schratt/miRNA_KD/Data_Output/mods_mmu_comp.rds")
    ahid <- rev(query(ah, c("EnsDb", "Mus musculus"))$ah_id)[2]
  }else if(species=="rno"){
    genome <- BSgenome.Rnorvegicus.UCSC.rn6::BSgenome.Rnorvegicus.UCSC.rn6
    if(is.null(mods)) mods <- readRDS(file = "/mnt/schratt/miRNA_KD/Data_Output/mods_rno_comp.rds")
    ahid <- rev(query(ah, c("EnsDb", "Rattus norvegicus"))$ah_id)[1]
  }
  ensdb <- ah[[ahid]]
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
  
  message("Scanning with ", cores, " cores")
  if(cores>1){
    BP <- MulticoreParam(cores, progress=TRUE)
  }else{
    BP <- SerialParam()
  }
  m <- findSeedMatches(seqs, mods, shadow=shadow, maxLogKd=maxLogKd, BP=BP, ...)
  
  if(is(m, "GRanges")) {
    metadata(m)$tx_info <- tx_info
    metadata(m)$ah_id <- ahid
  } else {
    attr(m, "tx_info") <- tx_info
    attr(m, "ah_id") <- ahid
  }
  if(is.null(save.path)) return(m)
  else saveRDS(m, file=save.path)
  rm(m)
  gc()
  message("Saved in: ", save.path)
}
