.get8merRange <- function(mod){
  stopifnot(is(mod,"KdModel"))
  mer8 <- getSeed8mers(mod$canonical.seed)
  mer8 <- which(mer8==mod$canonical.seed)
  fl <- rowSums(apply(scanMiR:::.flankingValues(), 2, FUN=range))
  fl*SampleKdModel$fl[mer8]+SampleKdModel$mer8[mer8]
}


