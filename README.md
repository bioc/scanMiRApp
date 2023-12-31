<img align="right" style="margin-left: 25px; margin-bottom: 30px;" src="https://raw.githubusercontent.com/ETHZ-INS/scanMiR/master/inst/docs/sticker.svg"/>

# scanMiRApp

`scanMiRApp` is a companion package to <a href="https://github.com/ETHZ-INS/scanMiR">scanMiR</a>.
It contains:

* A shiny interface to the <a href="https://github.com/ETHZ-INS/scanMiR">scanMiR</a> package,
enabling the scanning of transcripts (or custom sequences) for miRNA binding sites, the 
visualization of KdModels and binding results, as well as browsing predicted repression data.
* The `ScanMiRAnno` class, which encapsulates all the annotation information necessary for the
app/wrappers in this package for a given species.
* The <a href="vignettes/IndexedFST.Rmd">`IndexedFst`</a> class for fast indexed reading of 
large GenomicRanges or data.frames; this is used to enable fast random access to pre-compiled 
scans and aggregations without having to load them in memory.
* A number of convenience wrappers to `scanMiR` functions, for tasks like scanning or
visualizing sites on the UTR of a given transcript, running transcriptome-wide scans, or 
finding enriched miRNA-target pairs.

Click [here](https://ethz-ins.org/scanMiR/) for a live version of the app.
For more information, see the package's [vignettes](vignettes/) and the [preprint](https://doi.org/10.1101/2021.06.16.448293).

## Installation

First ensure that <a href="https://github.com/ETHZ-INS/scanMiR">scanMiR</a> and 
<a href="https://github.com/ETHZ-INS/scanMiRData">scanMiRData</a> are installed, then install 
with:

```{r}
BiocManager::install("ETHZ-INS/scanMiRApp")
```
