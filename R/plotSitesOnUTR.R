#' plotSitesOnUTR
#'
#' Wrapper function with minimal arguments to plot scanMiR-Binding
#' sites on 3'UTRs of specified transcripts. The red dashed line indicates the
#' background threshhold is indicated, the lightblue dashed line shows the
#' average 8mer dissociation rate of the given miRNA
#'
#'
#' @param species Either Mouse (="mmu"), Rat (="rno") or Human ("hsa")
#' @param transcriptID An ensembl TranscriptID
#' @param miRNA A miRNA name in the mirbase format (eg. "hsa-miR-485-5p")
#' @param label_6mers Logical whether to label 6mer sites in the plot
#' @param label_notes Logical whether to label special sites in the plot (as
#'   TDMD or Slicing)
#' @param verbose Logical; whether to print updates on the processing
#'
#' @return Returns a ggplot.
#' @import ggplot2
#' @export
plotSitesOnUTR <- function(species=NULL, transcriptID=NULL, miRNA=NULL, 
                           label_6mers=FALSE, label_notes=FALSE, verbose=TRUE){
  
  # Prepare everything & scan
  if(verbose) message("Prepare miRNA model")
  species <- match.arg(species, c("hsa","rno","mmu"))
  mods <- getKdModels(species = species, NULL)
  if(miRNA %in% names(mods)){
    mods <- mods[[miRNA]]
  }else if(length(w <- grep(miRNA,names(mods),ignore.case = TRUE))==1){
    mods <- mods[[w]]
  }else{
    stop("The specified microRNA is not listed for this species. Please check",
         "\nthe spelling (eg. 'hsa-mir-485-5p') and the organism")
  }
  if(verbose) message("Get Transcript Sequence")
  Seq <- getTranscriptSequence(tx = transcriptID ,species = species,UTRonly = TRUE)
  if(verbose) message("Scan")
  m <- findSeedMatches(seqs = Seq,seeds = mods, shadow = 15L, keepMatchSeq = TRUE,
                       p3.extra = TRUE, ret = "data.frame")
  
  # Prepare data.frame
  m$logKd <- m$log_kd / 1000
  m$type <- ifelse(m$type == "non-canonical","",as.character(m$type))
  if(isFALSE(label_6mers)) ifelse(grepl("6mer",m$type),"",as.character(m$type))
  m <- m[m$logKd < -1,]
  
  # get 8mer info
  mer8 <- getSeed8mers(mods$canonical.seed)
  wA <- which(substr(mer8,8,8)=="A")
  mer7 <- substr(mer8,1,7)
  As <- mods$mer8[wA]
  names(As) <- mer7[wA]
  mer.mean <- rowsum(mods$mer8[-wA],mer7[-wA])[,1]/3
  As <- As-mer.mean[names(As)]
  d <- data.frame(seed=names(mer.mean), base=mer.mean/-1000, 
                  "A"=As[names(mer.mean)]/-1000,
                  type=getMatchTypes(names(mer.mean),mods$canonical.seed),
                  row.names=NULL)
  d <- d[head(order(d$base+d$A, decreasing=TRUE),n=1),]
  mer8 <- d$base + d$A
  # title
  title <- paste0(transcriptID," - ",miRNA)
  
  p <- ggplot(m, aes(x = start, y = -`logKd`)) + 
    geom_hline(yintercept=1, linetype="dashed", color = "red", size=1) + 
    geom_hline(yintercept=mer8, linetype="dashed", color = "gray64", size=1) + 
    geom_point(size=2) + geom_text(label = m$type,nudge_y = -0.2) +
    xlab("sequence length") + ylab("-logKd") + xlim(0,width(Seq)) +
    theme_light() + ggtitle(title)
  if(label_notes){
    m$note <- ifelse(m$note == "-","",m$note)
    p <- p + geom_text(label = m$note,nudge_y = 0.2)
  } 
  p
}


