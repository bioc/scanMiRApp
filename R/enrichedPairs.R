#' enrichedMirTxPairs
#'
#' Identifies pairs of miRNA and target transcripts that have an unexpectedly
#' high number of sites.
#'
#' @param m A GRanges of matches, as produced by \code{\link{findSeedMatches}}.
#' This will be filtered down to only 8mer and 7mer sites.
#' @param minSites The minimum number of sites for a given miRNA-target pair to
#' be considered.
#' @param max.binom.p The maximum binomial p-value of miRNA-target pairs.
#'
#' @return A data.frame of top combinations, including number of sites and
#' the log-transformed binomial p-value.
#' @export
#'
#' @import Matrix
#' @importFrom stats pbinom
#' @examples
#' # we create a dummy scan (see `runFullScan`)
#' library(scanMiR)
#' seqs <- getRandomSeq(n=10)
#' mirs <- c("TTGTATAA","AGCATTAA")
#' m <- findSeedMatches(seqs,mirs,verbose=FALSE)
#' # we look for enriched pairs
#' res <- enrichedMirTxPairs(m)
#' res
enrichedMirTxPairs <- function(m, minSites=5, max.binom.p=0.001){
  m <- m[as.integer(m$type) %in% grep("8mer|7mer",levels(m$type))]
  b <- .matches2sparse(m)
  b <- b[rowSums(b>=minSites)>0,]
  rs <- rowSums(b)
  cs <- colSums(b)
  p <- as.matrix(rs/sum(rs)) %*% t(cs/sum(cs))
  S <- pbinom(as.matrix(b)-1, prob=p, size=rs, lower.tail=FALSE, log.p=TRUE)
  mode(S) <- "integer"
  S <- as(round(S), "sparseMatrix")
  S <- .sparse2df(S, "logp.binom")
  b <- .sparse2df(b, "sites")
  S$sites <- b$sites
  rm(b)
  S <- S[S$logp.binom < as.integer(log(max.binom.p)) &
           S$sites>=as.integer(minSites),]
  S[order(S$logp.binom),]
}


#' @import Matrix
.matches2sparse <- function(x){
  as(as.matrix(table(as.factor(seqnames(x)), as.factor(x$miRNA))),
     "sparseMatrix")
}
.sparse2df <- function(x, content="value"){
  dimn <- dimnames(x)
  x <- summary(x)
  w <- which(x$x!=ifelse(is.logical(x$x),FALSE,0))
  xout <- data.frame( feature=factor(x$i[w], levels=seq_len(length(dimn[[1]])),
                                     labels=dimn[[1]]),
                      set=factor(x$j[w], levels=seq_len(length(dimn[[2]])),
                                 labels=dimn[[2]]) )
  if(!is.logical(x$x)) xout[[content]] <- x$x[w]
  xout
}

