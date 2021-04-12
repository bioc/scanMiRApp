#' getTranscriptSequence
#' 
#' Utility wrapper to extracts the sequence of a given transcript (UTR or 
#' CDS+UTR).
#'
#' @param tx The ensembl ID of the transcript
#' @param species The species, either 'hsa', 'rno', or 'mmu'; ignored if 
#' `ensdbs` is given and is an `EnsDb` object.
#' @param ensdbs An \code{\link[ensembldb]{EnsDb}} object, or a named list of 
#' such objects.
#' @param genome The genome sequence (e.g. 
#' \code{\link[BSgenome]{BSgenome-class}}). If one of the three pre-defined 
#' `species`, will attempt to fetch the corresponding packages (if isntalled).
#' @param UTRonly Logical; whether to fetch only the UTR sequences.
#' @param ... Passed to \code{\link{AnnotationHub}}
#'
#' @return A \link{\code{DNAStringSet}}.
#' @export
#'
#' @importFrom GenomicFeatures threeUTRsByTranscript cdsBy extractTranscriptSeqs
#' @importFrom ensembldb metadata organism seqlevelsStyle
#' @importFrom AnnotationHub AnnotationHub query
#' @import Biostrings
#' @examples
#' # too long to run (needs fetching the annotation) :
#' # getTranscriptSequence("ENST00000641515", species="hsa")
getTranscriptSequence <- function(tx, species=NULL, ensdbs=NULL, genome=NULL,
                                  UTRonly=TRUE, ...){
  tx <- gsub("\\.[0-9]+$","",as.character(tx))
  if(is.null(ensdbs)){
    species <- match.arg(species, c("hsa","rno","mmu"))
    ah <- AnnotationHub(...)
    ahid <- switch(species,
                   hsa=rev(query(ah, c("EnsDb", "Homo sapiens"))$ah_id)[1],
                   mmu=rev(query(ah, c("EnsDb", "Mus musculus"))$ah_id)[2],
                   rno=rev(query(ah, c("EnsDb", "Rattus norvegicus"))$ah_id)[1],
                   stop("Species not among the pre-defined one, please provide",
                        " `ensdbs` and `genome` manually.")
    )
    ensdb <- ah[[ahid]]
    em <- metadata(ensdb)
    em <- setNames(em$value, em$name)
    message("Using ", em[["genome_build"]], ", Ensembl version ",
            em[["ensembl_version"]])
  }else{
    if(is.null(species) && length(ensdbs)==1){
      if(is.list(ensdbs)){
        stopifnot(!is.null(names(ensdbs)))
        species <- names(ensdbs)
        ensdb <- ensdbs[[1]]
      }else if(is(ensdbs,"EnsDb")){
        species <- organism(ensdbs)
        ensdb <- ensdbs
      }else{
        stop("`ensdbs` should either by a object of class EnsDb or a named ",
             "list of such objects.")
      }
    }else{
      species <- match.arg(species, names(ensdbs))
      ensdb <- ensdbs[[species]]
    }
  }
  if(is.null(genome)){
    em <- metadata(ensdb)
    em <- setNames(em$value, em$name)
    gbuild <- em[["genome_build"]]
    genome <- switch(gbuild,
      GRCh38=BSgenome.Hsapiens.UCSC.hg38:::BSgenome.Hsapiens.UCSC.hg38,
      GRCm38=BSgenome.Mmusculus.UCSC.mm10:::BSgenome.Mmusculus.UCSC.mm10,
      "Rnor_6.0"=BSgenome.Rnorvegicus.UCSC.rn6:::BSgenome.Rnorvegicus.UCSC.rn6,
      stop("Genome not among the pre-defined one, please provide `genome` ",
            "manually.")
    )
  }
  seqlevelsStyle(genome) <- "Ensembl"
  
  gr <- suppressWarnings(threeUTRsByTranscript(ensdb, filter=~tx_id==tx))
  gr <- gr[seqnames(gr) %in% seqlevels(genome)]
  if(length(gr)==0) stop("Transcript not found!")
  seqs <- extractTranscriptSeqs(genome, gr)
  if(!UTRonly){
    grl_ORF <- cdsBy(ensdb, by="tx", filter=~tx_id==tx)
    seqs_ORF <- extractTranscriptSeqs(genome, grl_ORF)
    orf.len <- setNames(lengths(seqs_ORF), names(seqs_ORF))
    seqs_ORF[names(seqs)] <- xscat(seqs_ORF[names(seqs)],seqs)
    seqs <- seqs_ORF
    rm(seqs_ORF)
    mcols(seqs)$ORF.length <- orf.len[names(seqs)]
  }
  seqs
}
