#' getTranscriptSequence
#'
#' Utility wrapper to extracts the sequence of a given transcript (UTR or
#' CDS+UTR).
#'
#' @param tx The ensembl ID of the transcript(s)
#' @param annotation A \code{\link{ScanMiRAnno}} object.
#' @param annoFilter An optional `AnnotationFilter` or `AnnotationFilterList` to
#' further filter the set of transcripts to be extracted
#' @param extract Which parts of the transcripts to extract. For `UTRonly`
#' (default) only the 3' UTR regions are extracted, `withORF` additionally
#' extracts the coding regions, and `exons` extracts all exons
#' @param ... Passed to \code{\link{AnnotationHub}}
#'
#' @return A \code{\link{DNAStringSet}}.
#' @export
#'
#' @importFrom GenomicFeatures exonsBy extractTranscriptSeqs cdsBy
#' threeUTRsByTranscript seqlevels<-
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom ensembldb metadata seqlevelsStyle seqlevelsStyle<-
#' @importFrom AnnotationFilter AnnotationFilterList AnnotationFilter
#' @importFrom stats setNames
#' @import Biostrings
#' @examples
#' anno <- ScanMiRAnno("fake")
#' seq <- getTranscriptSequence( tx="ENSTFAKE0000056456", annotation=anno )
getTranscriptSequence <- function(tx=NULL, annotation, annoFilter=NULL,
                                  extract=c("UTRonly", "withORF", "exons"),...){
  if(!is.null(tx) && is(annotation$ensdb,"EnsDb"))
    tx <- gsub("\\.[0-9]+$","",as.character(tx))
  if(!is.null(annoFilter)){
    if(!is(annotation$ensdb,"EnsDb")){
      warning("`annoFilter` is ignored when the annotation is not in ",
              "`EnsDb` format")
      annoFilter <- NULL
    }else{
      if(!is(annoFilter, "AnnotationFilterList") &&
         !is(annoFilter, "AnnotationFilter"))
        stop("filter must be either `AnnotationFilter` or `AnnotationFilterList`")
      if(!is.null(tx))
        annoFilter <- AnnotationFilterList(annoFilter, ~tx_id %in% tx)
    }
  }
  if(is.null(annoFilter)){
    if(!is.null(tx)){
      annoFilter <- AnnotationFilter(~tx_id %in% tx)
    } else {
      annoFilter <- AnnotationFilterList()
    }
  }
  ensdb <- annotation$ensdb
  extract <- match.arg(extract)
  if(is(ensdb,"EnsDb")){
    if(extract=="exons") {
      gr <- exonsBy(ensdb, by="tx", filter=annoFilter)
    } else {
      gr <- suppressWarnings(threeUTRsByTranscript(ensdb, filter=annoFilter))
    }
  }else{
    if(is.null(tx)){
      if(extract=="exons") {
        gr <- exonsBy(ensdb, by="tx", use.names=TRUE)
      } else {
        gr <- suppressWarnings(threeUTRsByTranscript(ensdb, use.names=TRUE))
      }
    }else{
      if(extract=="exons") {
        gr <- .getExonsFromTxDb(tx, ensdb)
      }else{
        gr <- .get3UTRFromTxDb(tx, ensdb)
      }      
    }
  }
  tx_strand <- unlist(unique(strand(gr)))
  genome <- annotation$genome
  if(is(ensdb,"EnsDb")) seqlevelsStyle(genome) <- "Ensembl"
  seqs <- DNAStringSet()
  if(length(gr)==0){
    if(extract=="withORF"){
      seqs <- DNAStringSet()
    } else {
      message("Nothing found!")
      return(DNAStringSet())
    }
  } else {
    gr <- gr[seqnames(gr) %in% seqlevels(genome)]
    seqs <- extractTranscriptSeqs(genome, gr)
  }
  if(extract=="withORF"){
    if(is(ensdb,"EnsDb")){
      grl_ORF <- suppressWarnings(cdsBy(annotation$ensdb, by="tx",
                                        filter=annoFilter))
    }else{
      if(is.null(tx)){
        grl_ORF <- suppressWarnings(cdsBy(annotation$ensdb, by="tx",
                                          use.names=TRUE))
      }else{
        grl_ORF <- .getExonsFromTxDb(tx, ensdb, fn=GenomicFeatures::cds)
      }
    }
    tx_strand <- unlist(unique(strand(grl_ORF)))
    if(length(grl_ORF)==0) {
      if(length(seqs) ==0) message("Nothing found!")
      return(seqs)
    }
    seqs_ORF <- extractTranscriptSeqs(genome, grl_ORF)
    orf.len <- setNames(lengths(seqs_ORF), names(seqs_ORF))
    seqs_ORF[names(seqs)] <- xscat(seqs_ORF[names(seqs)],seqs)
    seqs <- seqs_ORF
    rm(seqs_ORF)
    mcols(seqs)$ORF.length <- orf.len[names(seqs)]
  }
  mcols(seqs)$strand <- tx_strand[names(seqs)]
  seqs
}


#' @importFrom GenomicFeatures exons cds
.getExonsFromTxDb <- function(tx, db, fn=GenomicFeatures::exons){
  a <- fn(db,
          filter = list(tx_name = tx),
          columns = list("tx_name"))
  a$tx_name <- a$tx_name[which(a$tx_name %in% tx)]
  names(tx) <- tx <- unique(unlist(a$tx_name,use.names=FALSE))
  a <- GRangesList(lapply(tx, FUN=function(x){
    y <- a[any(a$tx_name==x)]
    y$tx_name <- x
    sort(y, decreasing=as.character(unlist(strand(y)))[1]=="-")
  }))
  names(a) <- tx
  a
}

.get3UTRFromTxDb <- function(tx, db){
  te <- .getExonsFromTxDb(tx, db, fn=GenomicFeatures::exons)
  if(length(te)==0) return(NULL)
  tcds <- .getExonsFromTxDb(tx, db, fn=GenomicFeatures::cds)
  i <- intersect(names(tcds[lengths(tcds)>0]), names(te[lengths(te)>0]))
  if(length(i)==0) return(NULL)
  te <- GRangesList(mapply(te=te[i], tcds=tcds[i], FUN=function(te,tcds){
    te <- setdiff(te,tcds)
    if(as.character(unlist(strand(te),use.names=FALSE)[1])=="-"){
      return(te[end(te)<=min(start(tcds))])
    }
    te[start(te)>=max(end(tcds))]
  }))
  te
}


#' plotSitesOnUTR
#'
#' Wrapper function with minimal arguments to plot scanMiR-Binding
#' sites on 3'UTRs of specified transcripts. The red dashed line indicates the
#' background threshhold is indicated, the lightblue dashed line shows the
#' average 8mer dissociation rate of the given miRNA
#'
#' @param tx An ensembl TranscriptID
#' @param annotation A \code{\link{ScanMiRAnno}} object.
#' @param miRNA A miRNA name in the mirbase format (eg. "hsa-miR-485-5p"), a
#' `KdModel`, or a miRNA sequence or target seed.
#' @param label_6mers Logical whether to label 6mer sites in the plot
#' @param label_notes Logical whether to label special sites in the plot (as
#'   TDMD or Slicing)
#' @param verbose Logical; whether to print updates on the processing
#' @param ... Any further arguments passed to
#' \code{\link[scanMiR]{findSeedMatches}}
#'
#' @return Returns a ggplot.
#' @importFrom ggplot2 ggplot geom_hline geom_point geom_text labs xlim aes
#' theme_light
#' @import scanMiR
#' @export
#' @examples
#' anno <- ScanMiRAnno("fake")
#' plotSitesOnUTR( tx="ENSTFAKE0000056456", annotation=anno,
#'                 miRNA="hsa-miR-155-5p" )
plotSitesOnUTR <- function(tx, annotation, miRNA = NULL, label_6mers = FALSE, 
                           label_notes = FALSE, verbose = TRUE, ...){
  stopifnot(is(annotation, "ScanMiRAnno"))
  if (verbose) message("Prepare miRNA model")
  mods <- annotation$models
  if((is(miRNA, "character") && length(miRNA) == 1 && 
     nchar(gsub("[ACGTU]", "", miRNA)) == 0) || is(miRNA, "KdModel")){
    mods <- miRNA
    if(!(is(miRNA,"KdModel"))){
      title <- paste0(tx, " - ", miRNA)
    }
  } else if(miRNA %in% names(mods)){
    mods <- mods[[miRNA]]
  } else if(length(w <- grep(miRNA, names(mods), ignore.case=TRUE)) == 1){
    mods <- mods[[w]]
  } else {
    stop("The specified microRNA is not listed for this species. Please check", 
         "\nthe spelling (eg. 'hsa-mir-485-5p') and the organism")
  }
  if (verbose) 
    message("Get Transcript Sequence")
  Seq <- getTranscriptSequence(tx = tx, annotation = annotation)
  if (verbose) 
    message("Scan")
  m <- findSeedMatches(seqs = Seq, seeds = mods, shadow = 15L, 
                       keepMatchSeq = TRUE, p3.extra = TRUE, ret = "data.frame", 
                       verbose = FALSE, ...)
  if (isFALSE(label_6mers)){ 
    m$type <- ifelse(grepl("6mer", m$type), "", as.character(m$type))}
  if(is(miRNA, "KdModel")){
    m$logKd <- m$log_kd/1000
    m$type <- ifelse(m$type == "non-canonical", "", as.character(m$type))
    m <- as.data.frame(m[m$logKd < -1, ])
    mer8 <- get8merRange(mods)/-1000
    max8 <- max(mer8)
    title <- paste0(tx, " - ",mods$name)
    p <- ggplot(m, aes(x = start, y = -logKd)) + 
      geom_hline(yintercept = 1,linetype = "dashed", color = "red", size = 1) + 
      geom_point(size = 2) + geom_text(label = m$type, nudge_y = -max8/50) + 
      labs(x = "sequence length", y = "-logKd", title = title) +  
      theme_light() + xlim(0, width(Seq)) + expand_limits(y=c(0,max8)) + 
      geom_rect(
               xmin = -Inf, xmax = Inf, 
               ymin = min(mer8), ymax = max8,  fill = "chartreuse3", alpha=.3)
    
    if (label_notes){ 
      p <- p + geom_text(data = m[m$note != "-", ], label = aes(note), 
                         nudge_y = max8/50)
    }
  }else{
    m <- m[m$type != "",]
    p <- ggplot(m, aes(x = start, y = type)) + 
      geom_point(size = 2) + 
      labs(x = "sequence length", y = "site type", title = title) + xlim(0, width(Seq)) + 
      theme_light()
  }
  p
}


