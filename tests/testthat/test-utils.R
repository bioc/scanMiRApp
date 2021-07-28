anno <- ScanMiRAnno("fake")

test_that("runFullScan works and gives the expected output", {
  m1 <- runFullScan(anno, onlyCanonical=TRUE)
  m2 <- runFullScan(anno, extract="withORF", onlyCanonical=FALSE, maxLogKd=-0.5)
  expect_equal(seqlevels(m1), "ENSTFAKE0000056456")
  expect_equal(seqlevels(m2), "ENSTFAKE0000056456")
  expect_equal(table(m1$type)[["8mer"]],1L)
  expect_equal(table(m2$type)[["8mer"]],1L)
  expect_equal(table(m2$type)[["6mer-m8"]],1L)
  expect_gt(table(m2$type)[["non-canonical"]],0)
  expect_equal(table(m1$type)[["non-canonical"]],0L)
  expect_equal(length(unique(m2$ORF)),2L)
})

test_that("indexedFst works and maintains data integrity", {
  f <- file.path(tempdir(), "ifst.test")
  d <- data.frame( category=sample(LETTERS[1:4], 10000, replace=TRUE),
                   var2=sample(LETTERS, 10000, replace=TRUE),
                   var3=runif(10000) )
  saveIndexedFst(d, "category", f)
  d <- d[order(d$category),]
  d2 <- loadIndexedFst(f)
  d3 <- as.data.frame(d2)
  expect_equal(d$category, d3$category)
  expect_equal(d$var2, d3$var2)
  expect_equal(d$var3, d3$var3)
  expect_equal(nrow(d[d$category=="A",]), nrow(d2$A))
})

test_that("Pair enrichment analysis works", {
  eg <- expand.grid(target=LETTERS[2:5], mir=letters[1:5])
  i <- rpois(nrow(eg),2)
  i[i>7] <- 7
  eg <- eg[rep(seq_len(nrow(eg)), i),]
  eg <- rbind(eg, cbind(target=rep("A",12),mir=rep("a",12)))
  eg$type <- sample(c("8mer","7mer"), nrow(eg), replace=TRUE)
  eg$start <- 1L+sample.int(100,nrow(eg))
  gr <- GRanges(eg$target, IRanges(eg$start, eg$start+8), miRNA=eg$mir,
                type=factor(eg$type) )
  ep <- enrichedMirTxPairs(gr, max.binom.p=1, minSites=1)
  expect(ep[1,1]=="A" && ep[1,2]=="a", failure_message=
    "enrichedMirTxPairs does not report spiked result as top")
})

test_that("UTR plot works", {
  plotSitesOnUTR(tx="ENSTFAKE0000056456", annotation=anno,
                 miRNA="hsa-miR-155-5p")
  plotSitesOnUTR(tx="ENSTFAKE0000056456", annotation=anno,
                 anno$models[[1]]$canonical.seed)
})


