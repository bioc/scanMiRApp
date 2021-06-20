#' @exportClass ScanMiRAnno
setClass(
  "ScanMiRAnno",
  contains="list",
  validity=function(object){
    stopifnot(!is.null(object$genome) &&
              (is(object$genome, "BSgenome") || is(object$genome, "TwoBitFile")))
    stopifnot(!is.null(object$ensdb) &&
                (is(object$ensdb, "EnsDb") || is(object$ensdb, "TxDb")))
    if(!is.null(object$models))
      stopifnot(is(object$models, "KdModelList") ||
                  is(object$models, "character"))
    if(!is.null(object$scan))
      stopifnot(is(object$scan, "GRanges") ||
                  is(object$scan, "IndexedFst"))
    if(!is.null(object$aggregated))
      stopifnot(is(object$aggregated, "data.frame") ||
                  is(object$aggregated, "IndexedFst"))
  }
)

#' ScanMiRAnno
#'
#' @param species The species/build acronym for automatic construction; if
#' omitted, `genome` and `ensdb` should be given. Current possible values are:
#' GRCh38, GRCm38, Rnor_6.
#' @param genome A \link{BSgenome}, or a \code{\link[rtracklayer]{TwoBitFile}}
#' @param ensdb An EnsDb object
#' @param models An optional KdModelList
#' @param scan An optional full scan (IndexedFst or GRanges)
#' @param aggregated An optional per-transcript aggregation (IndexedFst or
#' @param version optional ensembl version
#' data.frame)
#' @param addDBs A named list of additional tx-miRNA databases, each of which
#' should be a data.frame with the columns 'transcript', 'miRNA', and 'score'.
#' @param ... Arguments passed to `AnnotationHub`
#'
#' @return A `ScanMiRAnno` object
#' @export
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom ensembldb getGenomeTwoBitFile
#' @importFrom scanMiRData getKdModels
#' @examples
#' anno <- ScanMiRAnno(species="fake")
#' anno
ScanMiRAnno <- function(species=NULL, genome=NULL, ensdb=NULL, models=NULL,
                        scan=NULL, aggregated=NULL, version=NULL,
                        addDBs=list(), ...){
  blink <- NULL
  if(!is.null(species)){
    stopifnot(is.null(genome) && is.null(ensdb))
    species <- match.arg(species, c("GRCh38","GRCm38","Rnor_6","fake"))
    if(species=="fake") return(.fakeAnno())
    if(is.null(models))
      models <- switch(species,
                       GRCh38=scanMiRData::getKdModels("hsa"),
                       GRCm38=scanMiRData::getKdModels("mmu"),
                       GRCm39=scanMiRData::getKdModels("mmu"),
                       "Rnor_6"=scanMiRData::getKdModels("rno"),
                       stop("Species not among the pre-defined one, please",
                            "provide `models` manually.")
                      )
    ah <- AnnotationHub(...)
    ensdb <- ah[[rev(query(ah, c("EnsDb", species, version))$ah_id)[1]]]
    genome <- switch(species,
      GRCh38=BSgenome.Hsapiens.UCSC.hg38:::BSgenome.Hsapiens.UCSC.hg38,
      GRCm38=BSgenome.Mmusculus.UCSC.mm10:::BSgenome.Mmusculus.UCSC.mm10,
      "Rnor_6"=BSgenome.Rnorvegicus.UCSC.rn6:::BSgenome.Rnorvegicus.UCSC.rn6,
      getGenomeTwoBitFile(ensdb)
    )
    blink <- switch( species,
      GRCh38="https://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=",
      GRCm38="https://www.ensembl.org/Mus_musculus/Gene/Summary?db=core;g=",
      "Rnor_6"="https://www.ensembl.org/Rattus_norvegicus/Gene/Summary?db=core;g=",
      NULL )
  }
  new("ScanMiRAnno", list(
    genome=genome, ensdb=ensdb, models=models, scan=scan, aggregated=aggregated,
    addDBs=addDBs, ensembl_gene_baselink=blink
  ))
}


#' Methods for the \code{\link{ScanMiRAnno}} class
#' @name ScanMiRAnno-methods
#' @rdname ScanMiRAnno-methods
#' @aliases ScanMiRAnno-methods
#' @seealso \code{\link{ScanMiRAnno}}
#' @param object An object of class \code{\link{ScanMiRAnno}}
#' @return Depends on the method.
NULL

#' @rdname ScanMiRAnno-methods
#' @aliases ScanMiRAnno-methods
#' @export
setMethod("summary", "ScanMiRAnno", function(object){
  show(object$genome)
  cat("\n")
  show(object$ensdb)
  cat("\nModels:\n")
  if(is.null(object$models)){
    cat("None\n")
  }else{
    summary(object$models)
  }
  cat("\nScan:\n")
  if(is.null(object$scan)){
    cat("None\n")
  }else{
    summary(object$scan)
  }
  cat("\nAggregated:\n")
  if(is.null(object$aggregated)){
    cat("None\n")
  }else{
    summary(object$aggregated)
  }
  if(!is.null(object$addDBs) && length(object$addDBs)>0){
    cat("\nAdditional DBs: ", paste(names(object$addDBs),collapse=", "), "\n")
  }
})

#' @rdname ScanMiRAnno-methods
#' @aliases ScanMiRAnno-methods
#' @export
#' @importFrom methods show
setMethod("show", "ScanMiRAnno", function(object){
  if(is(object$genome, "BSgenome")){
    gm <- paste0(metadata(object$genome)$organism, " (",
                 metadata(object$genome)$genome, ")")
  }else{
    gm <- toString(object$genome)
  }
  em <- setNames(metadata(object$ensdb)$value, metadata(object$ensdb)$name)
  if(is(object$ensdb,"EnsDb")){
    em <- paste0(em[["Organism"]], " (", em[["genome_build"]],") v",
                 em[["ensembl_version"]])
  }else{
    em <- paste0(em[["Organism"]], " (", em[["Genome"]],")")
  }
  cat("Genome:", gm, "\nAnnotation:", em)
  if(!is.null(object$models)) cat("\nModels:", class(object$models),
                                  "of length", length(object$models))
  if(!is.null(object$scan)) cat("\n + Scan")
  if(!is.null(object$aggregated)) cat("\n + Aggregated")
  cat("\n")
})


#' @importFrom Biostrings DNAStringSet
#' @importFrom rtracklayer TwoBitFile export.2bit
#' @importFrom GenomicFeatures makeTxDbFromGRanges
#' @importFrom utils data
.fakeAnno <- function(){
  data("SampleKdModel", package="scanMiR", envir = environment())
  data("SampleTranscript", package="scanMiR", envir = environment())
  f <- tempfile()
  export.2bit(DNAStringSet(SampleTranscript), f)
  ge <- TwoBitFile(f)
  gr <- GRanges(seqlevels(ge), IRanges(c(1,1,1,1),c(898,898,898,210)),
                type=c("gene","transcript","exon","CDS"),
                exon_id=c(NA,NA,"e1",NA), exon_number=c(NA,NA,1,NA),
                strand="+", source="fake", gene_biotype="protein_coding")
  gr$entrezid <- gr$gene_id <- gr$gene_name <- "gene1"
  gr$transcript_id <- gr$tx_id <- "ENSTFAKE0000056456"
  md <- data.frame(name=c("Organism","Genome"), value=c("Fake falsus","fake1"))
  db <- makeTxDbFromGRanges(gr, metadata=md)
  ScanMiRAnno(genome=ge, ensdb=db, models=c(SampleKdModel))
}
