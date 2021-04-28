#' @exportClass ScanMiRAnno
setClass(
  "ScanMiRAnno",
  contains="list",
  validity=function(object){
    stopifnot(!is.null(object$genome) && is(object$genome, "BSgenome"))
    stopifnot(!is.null(object$ensdb) && is(object$ensdb, "EnsDb"))
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
#' @param genome A BSgenome
#' @param ensdb An EnsDb object
#' @param models An optional KdModelList
#' @param scan An optional full scan (IndexedFst or GRanges)
#' @param aggregated An optional per-transcript aggregation (IndexedFst or
#' data.frame)
#' @param ... Arguments passed to `AnnotationHub`
#'
#' @return A `ScanMiRAnno` object
#' @export
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom scanMiRData getKdModels
#' @examples
#' anno <- ScanMiRAnno(species="Rnor_6")
#' anno
ScanMiRAnno <- function(species=NULL, genome=NULL, ensdb=NULL, models=NULL,
                        scan=NULL, aggregated=NULL, ...){
  if(!is.null(species)){
    stopifnot(is.null(genome) && is.null(ensdb))
    species <- match.arg(species, c("GRCh38","GRCm38","Rnor_6"))
    ah <- AnnotationHub(...)
    ensdb <- ah[[rev(query(ah, c("EnsDb", species))$ah_id)[1]]]
    genome <- switch(species,
      GRCh38=BSgenome.Hsapiens.UCSC.hg38:::BSgenome.Hsapiens.UCSC.hg38,
      GRCm38=BSgenome.Mmusculus.UCSC.mm10:::BSgenome.Mmusculus.UCSC.mm10,
      "Rnor_6"=BSgenome.Rnorvegicus.UCSC.rn6:::BSgenome.Rnorvegicus.UCSC.rn6,
      stop("Genome not among the pre-defined one, please provide `genome` ",
          "manually.")
    )
    if(is.null(models)) models <- switch(species,
      GRCh38=scanMiRData::getKdModels("hsa"),
      GRCm38=scanMiRData::getKdModels("mmu"),
      "Rnor_6"=scanMiRData::getKdModels("rno"))
  }
  new("ScanMiRAnno", list(
    genome=genome, ensdb=ensdb, models=models, scan=scan, aggregated=aggregated
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
})

#' @rdname ScanMiRAnno-methods
#' @aliases ScanMiRAnno-methods
#' @export
#' @importFrom methods show
setMethod("show", "ScanMiRAnno", function(object){
  gm <- paste0(metadata(object$genome)$organism, " (",
               metadata(object$genome)$genome, ")")
  em <- setNames(metadata(object$ensdb)$value, metadata(object$ensdb)$name)
  em <- paste0(em[["Organism"]], " (", em[["genome_build"]],") v",
               em[["ensembl_version"]])
  cat("Genome:", gm, "\nAnnotation:", em)
  if(!is.null(object$models)) cat("\nModels:", class(object$models),
                                  "of length", length(object$models))
  if(!is.null(object$scan)) cat("\n + Scan")
  if(!is.null(object$aggregated)) cat("\n + Aggregated")
  cat("\n")
})
