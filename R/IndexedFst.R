#' IndexedFst
#'
#' Objects of the IndexedFst class enable fast named random access to FST files.
#' This is particularly appropriate for large data.frames which often need to
#'  be accessed according to the (e.g. factor) value of a particular column.
#'
#' @examples
#' # we first create and save an indexed FST file
#' tmp <- tempdir()
#' f <- system.file(tmp, "test")
#' d <- data.frame( category=sample(LETTERS[1:4], 10000, replace=TRUE),
#'                  var2=sample(LETTERS, 10000, replace=TRUE),
#'                  var3=runif(10000) )
#' format(object.size(d),units="Kb")
#' saveIndexedFst(d, "category", f)
#' rm(d)
#' # we then load the index, and can use category names for random access:
#' d <- loadIndexedFst(f)
#' format(object.size(d),units="Kb")
#' nrow(d)
#' names(d)
#' head(d$A)
#' @seealso \code{\link{saveIndexedFst}}, \code{\link{loadIndexedFst}}
#' @export
#' @import methods
#' @exportClass IndexedFst
#' @author Pierre-Luc Germain, \email{pierre-luc.germain@@hest.ethz.ch}
#' @aliases IndexedFst-class IndexedFst
#' @return Depends on the method
#' @rdname IndexedFst-class
#' @name IndexedFst-class
setClass(
  "IndexedFst",
  slots=representation(
    fst.file="character",
    index="data.frame",
    add.info="list",
    nthreads="integer"
  ),
  prototype=prototype(fst.file=NA_character_, index=data.frame(), nthreads=1L,
                      add.info=list()),
  validity=function(object){
    if(length(object@nthreads)!=1 || !is.integer(object@nthreads) ||
       object@nthreads<1)
      stop("`nthreads` should be a positive integer")
    if(length(object@fst.file)!=1)
      stop("fst.file should be a character of length 1.")
    if(!file.exists(object@fst.file)) stop("FST file does not exist!")
    TRUE
  }
)

setMethod("initialize", "IndexedFst", function(.Object, ...) {
  o <- callNextMethod(.Object, ...)
  o@fst.file <- normalizePath(o@fst.file)
  o@nthreads <- as.integer(o@nthreads)
  ff <- gsub("\\.fst$",".idx.rds",o@fst.file)
  tryCatch({
    ff <- readRDS(ff)
    o@index <- ff$index
    ff$index <- NULL
    o@add.info <- ff
  }, error=function(e) stop("Could not find or read index file."))
  validObject(o)
  return(o)
})


#' @rdname IndexedFst-class
#' @importMethodsFrom methods show
#' @param object an IndexedFst object
#' @export
setMethod("show", "IndexedFst", function(object){
  paste0(object@fst.file, " (",nrow(object@index)," sets)")
})

#' @rdname IndexedFst-class
#' @importFrom fst fst.metadata
#' @export
setMethod("summary", "IndexedFst", function(object){
  fst.metadata(object@fst.file)
})

#' @rdname IndexedFst-class
#' @param x an IndexedFst object
#' @export
setMethod("names", "IndexedFst", function(x){
  row.names(x@index)
})

#' @rdname IndexedFst-class
#' @export
setMethod("length", "IndexedFst", function(x){
  nrow(x@index)
})

#' @rdname IndexedFst-class
#' @export
setMethod("lengths", "IndexedFst", function(x){
  y <- x@index[,2]-x@index[,1]+1
  names(y) <- names(x)
  y
})

setGeneric("nrow")
#' @rdname IndexedFst-class
#' @export
setMethod("nrow", "IndexedFst", function(x){
  max(x@index[,2])
})

setGeneric("ncol")
#' @rdname IndexedFst-class
#' @export
setMethod("ncol", "IndexedFst", function(x){
  length(fst.metadata(x@fst.file)$columnNames)
})

setGeneric("colnames")
#' @rdname IndexedFst-class
#' @export
setMethod("colnames", signature("IndexedFst"), function(x){
  fst.metadata(x@fst.file)$columnNames
})

#' @rdname IndexedFst-class
#' @param i the desired index (either numeric or name)
#' @param j ignored
#' @param ... ignored
#' @export
setMethod("[[", signature("IndexedFst"), function(x, i, j=NULL, ...){
  if(is.numeric(i)){
    name <- names(x)[i]
  }else{
    name <- i
  }
  .fst.read.wrapper(x, name)
})

#' @rdname IndexedFst-class
#' @export
setMethod("[", signature("IndexedFst"), function(x, i, j=NULL, ...){
  if(is.logical(i)) i <- which(i)
  if(is.numeric(i)){
    name <- names(x)[i]
  }else{
    name <- i
  }
  .fst.read.wrapper(x, name)
})

#' @rdname IndexedFst-class
#' @param name the indexed name to fetch
#' @export
setMethod("$", "IndexedFst", definition = function(x, name){
  .fst.read.wrapper(x, match.arg(name, row.names(x@index)))
})

#' @rdname IndexedFst-class
#' @param n the desired number of rows
#' @export
setMethod("head", "IndexedFst", definition = function(x, n=6L, ...){
  if(!is.numeric(n) || !(n>0) || n!=as.integer(n))
    stop("`n` should be a positive integer.")
  w <- which(x@index[,2]>=n)
  if(length(w)==0) return(.fst.read.wrapper(x))
  head(x[seq_len(w[1])], n=n)
})

#' @rdname IndexedFst-class
#' @export
setMethod("as.data.frame", "IndexedFst", definition=function(x, name) {
  .fst.read.wrapper(x, convertGR=FALSE)
})

#' @importFrom fst threads_fst read.fst
#' @importFrom data.table rbindlist
#' @importFrom IRanges IRanges
#' @importFrom GenomicRanges GRanges
.fst.read.wrapper <- function(x, names=NULL, convertGR=TRUE){
  if(!is(x,"IndexedFst")) stop("`x` should be an IndexedFst object.")
  if(!is.null(names) && length(names)>1){
    f <- lapply(names, FUN=function(i) .fst.read.wrapper(x,i,convertGR=FALSE))
    if(length(f)>5){
      f <- as.data.frame(data.table::rbindlist(f))
    }else{
      f <- do.call(rbind, f)
    }
  }else{
    ont <- threads_fst()
    threads_fst(x@nthreads)
    if(is.null(names)){
      f <- read.fst(x@fst.file)
    }else{
      f <- read.fst(x@fst.file, from=x@index[names,1], to=x@index[names,2])
    }
    threads_fst(ont)
  }
  if(convertGR && !is.null(x@add.info$isGR) && x@add.info$isGR){
    f <- tryCatch(.df2gr(f, x@add.info), error=function(e){
      warning("Could not convert to GRanges:", e)
      f
    })
  }
  f
}

#' Saving and loading IndexedFst
#'
#' Functions to save or load and indexed \code{\link[fst]{fst}} file
#'
#' @param file Path to the fst file, it's index (.idx), or their prefix.
#' @param nthreads Number of threads to use for reading (default 1). This does
#' not affect the loading of the index itself, but will affect all downstream
#' reading operations performed on the object. If NULL, will use
#' `fst::threads_fst()`.
#'
#' @return `loadIndexedFst` returns an object of class
#' \code{\link{IndexedFst-class}}, and `saveIndexedFst` returns nothing.
#' @seealso \code{\link{IndexedFst-class}}
#' @rdname save-load-IndexedFst
#' @aliases loadIndexedFst
#' @export
#' @examples
#' # we first create and save an indexed FST file
#' tmp <- tempdir()
#' f <- system.file(tmp, "test")
#' d <- data.frame( category=sample(LETTERS[1:4], 10000, replace=TRUE),
#'                  var2=sample(LETTERS, 10000, replace=TRUE),
#'                  var3=runif(10000) )
#' saveIndexedFst(d, "category", f)
#' # we then load the index, and can use category names for random access:
#' d <- loadIndexedFst(f)
loadIndexedFst <- function(file, nthreads=1){
  if(grepl("\\.fst$",file)) return(new("IndexedFst", fst.file=file))
  if(grepl("\\.idx\\.rds$",file))
    return(new("IndexedFst", fst.file=gsub("idx\\.rds$","fst",file)))
  new("IndexedFst", fst.file=paste0(file,".fst"), nthreads=as.integer(nthreads))
}

#' saveIndexedFst
#'
#' Saves a data.frame (or GRanges object) into an indexed FST file.
#'
#' @param d A data.frame or \code{\link[GenomicRanges]{GRanges}} object
#' @param index.by A column of `d` by which it should be indexed.
#' @param file.prefix Path and prefix of the output files.
#' @param index.properties An optional data.frame of properties, with the levels
#' of `index.by` as row names.
#' @param add.info An optional list of additional information to save.
#' @param ... Passed to `write.fst`
#' @rdname save-load-IndexedFst
#' @aliases saveIndexedFst
#' @seealso \code{\link{IndexedFst-class}}
#' @importFrom fst write.fst threads_fst
#' @export
saveIndexedFst <- function(d, index.by, file.prefix, nthreads=1,
                           index.properties=NULL, add.info=list(), ...){
  if(is(d, "GRanges")){
    add.info <- c(add.info, metadata(d))
    d <- as.data.frame(d)
    d$end <- NULL
    if(all(d$strand=="*")){
      d$strand <- NULL
    }else if(length(unique(d$strand))==1){
      add.info$strand <- as.factor(d$strand[1])
      d$strand <- NULL
    }
    if(length(unique(d$width))==1){
      add.info$width <- d$width[1]
      d$width <- NULL
    }
    if(length(unique(d$seqnames))==1){
      add.info$seqnames <- d$seqnames[1]
      d$seqnames <- NULL
    }
    add.info$isGR <- TRUE
  }
  if(!is.data.frame(d)) stop("`d` should be a data.frame.")
  if((!is.character(index.by) && !is.integer(index.by)) ||
     length(index.by)!=1 || is.null(d[[index.by]]))
    stop("`index.by` should be a scalar character or integer indicating the ",
         "column by which to index.")
  if(is.numeric(d[[index.by]]))
    warning("The `index.by` column selected is numeric!\n",
            "If the column is really continuous, rather than having many ",
            "repeated values, this might lead to very suboptimal behaviors!")
  if(!is.null(nthreads)) threads_fst(nthreads)
  d <- d[order(d[[index.by]]),]
  file.prefix <- gsub("\\.fst$","",file.prefix)
  w <- which(!duplicated(d[[index.by]]))
  idx <- data.frame(row.names=d[[index.by]][w], start=w,
                    end=c(w[-1]-1L,nrow(d)))
  if(!is.null(index.properties)){
    idx <- merge(idx, as.data.frame(index.properties),
                 by="row.names", all.x=TRUE)
    row.names(idx) <- idx[,1]
    idx <- idx[,-1]
  }
  idx <- c(list(index=idx), add.info)
  saveRDS(idx, paste0(file.prefix, ".idx.rds"))
  write.fst(d, paste0(file.prefix, ".fst"), ...)
}

#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom S4Vectors metadata metadata<-
.df2gr <- function(x, add.info=list()){
  if(is.null(x$end)){
    if(is.null(x$width)){
      stopifnot(!is.null(add.info$width))
      x$end <- x$start + add.info$width-1L
    }else{
      x$end <- x$start + x$width-1L
    }
  }
  x$width <- NULL
  if(is.null(x$strand) && !is.null(add.info$strand))
    x$strand <- add.info$strand
  if(is.null(x$seqnames) && !is.null(add.info$seqnames))
    x$seqnames <- add.info$seqnames
  x <- makeGRangesFromDataFrame(x, keep.extra.columns = TRUE)
  add.info <- add.info[setdiff(names(add.info),
                               c("strand","width","seqnames","isGR"))]
  metadata(x) <- add.info
  x
}
