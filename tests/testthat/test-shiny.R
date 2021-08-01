suppressPackageStartupMessages({
  library(shiny)
  library(shinytest)
  library(scanMiRApp)
})

anno <- ScanMiRAnno("fake")

test_that("shiny app (server) works", {
  skip_on_cran()
  testServer(scanMiRserver(list(fake=anno)), {
    session$setInputs(gene="gene1", transcript="ENSTFAKE0000056456",
                      mirlist="fake", seqFeature="3' UTR only",
                      subjet_type="transcript", annotation="fake", shadow=15L,
                      keepmatchseq=FALSE, minDist=7L, maxLogKd=-1L,
                      scanNonCanonical=TRUE, circular=FALSE)
    session$setInputs(mirnas="hsa-miR-155-5p", mirnas_all=FALSE)
    expect_equal(selgene(),"gene1")
    expect_equal(seltx(),"ENSTFAKE0000056456")
    expect_equal(names(selmods()), "hsa-miR-155-5p")
    cs <- checksum()
    cached.hits[[checksum()]] <- do.scan()
    current.cs(cs)
    session$flushReact()
    h <- hits()$hits
    expect_equal(start(h), c(281L,482L))
    expect_equal(as.character(h$type), c("8mer","7mer-m8"))
    p <- output$manhattan
    selectedMatch(h[1])
    session$flushReact()
    expect(grepl("||||||||||||     ||||||||", output$alignment, fixed=TRUE),
           "viewTargetAlignment doesn't seem to work as expected", info)

  })
})
