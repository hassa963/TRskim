
<!-- README.md is generated from README.Rmd. Please edit that file -->

# TRskim

<!-- badges: start -->

<!-- badges: end -->

A tandem repeat comparison and visualization tool.

## Description

`TRskim` is an R package that visualizes tandem repeat composition by
decomposing sequences into known/novel motifs, encoding then aligning
them, and generates visualizations to compare motif patterns across
alleles. The tool is intended for those researching tandem repeats and
would like to perform a motif composition analysis across alleles.

TRskim provides motif-level decomposition and visualization of tandem
repeats, enabling researchers to inspect the internal structure of
repeat arrays across alleles. While most tools report only repeat
length, TRskim presents a potential motif composition, detects
interruptions, and highlights blended versus pure repeat patterns. This
level of resolution is essential for studying pathogenic repeat
expansions and assessing repeat instability, as motif structure
influence disease relevance, mutational dynamics, and secondary
structure formation. While there are libraries available in python
(TRviz) that perform a similar workflow, TRskim circumnavigates the need
to port in a separate coding language and can easily be included in R
workflows.

TRskim was developed in the following environment: R version: 4.5.1
(2025-06-13) - “Great Square Root” Platform: aarch64-apple-darwin20

## Installation

To install the latest version of the package:

``` r
install.packages("devtools")
library("devtools")
devtools::install_github("hassa963/TRskim", build_vignettes = TRUE)
library("TRskim")
```

To run the Shiny app:

``` r
TRskim::runTRskim()
```

## Overview

``` r
ls("package:TRskim")
browseVignettes("TRskim")
```

`TRskim` contains 6 functions:

1.  `decomposeTRs` for decomposing tandem repeats into their motifs.

2.  `encodeTRs` for encoding decomposed tandem repeats into their motifs
    such that each motif has a one character symbol. **Can only encode a
    maximum of 66 different motifs**

3.  `alignTRs` for aligning tandem repeats by their motifs.

4.  `plotTR` for visualizing aligned tandem repeats by their motif
    composition.

5.  `runWorkflow` for performing all the steps of the TRskim workflow
    with only one function call.

6.  `runTRskim` for launching the TRskim shiny app.

The package also contains example DNAStringSet objects of tandem repeats
and motifs, named SORL1 and motifs_SORL1, respectively. These can be
used as input to explore the TRskim workflow and test functions such as
decomposeTRs(), encodeTRs(), and plotTR().

For detailed guidance on how to use these example datasets, refer to the
package vignette. Below is an overview of the TRskim package

![](./inst/extdata/TRskim_Overview.png)

## Contributions

The author of TRskim is Nour Hassan. TRskim contains six functions, all
written by the author, which utilize functions from other R packages.

The algorithm for decomposing tandem repeats in decomposeTRs is based on
the approach used in the Python library TRviz (Park et al., 2023).
ChatGPT was used to dissect the TRviz algorithm and create a skeleton
for implementing a similar approach in R. The function leverages the
Biostrings package, using DNAString and DNAStringSet objects along with
the matchPattern() function to identify motif matches within alleles.

For alignTRs, ChatGPT assisted in dissecting the logic for aligning
tandem repeats by motif. The plotTR function generates tile and bar plot
visualizations of the tandem repeats using ggplot2, while reshape2 is
used to format the data for plotting. ChatGPT was also used to debug the
code for all functions and to assist in designing the visualizations for
plotTR. Testing of functions was supported with guidance from ChatGPT.

Example datasets (SORL1 and motifs_SORL1) were taken from the TRviz
Python library (Park et al., 2023) and are included in the package under
inst/extdata .

## References

- BioRender.com. BioRender \[Online\]. Available at:
  <https://www.biorender.com> (accessed 26 October 2025).

- Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. Shiny: Web
  Application Framework for R — Tabsets example. Shiny Gallery, RStudio
  (2025). <https://shiny.posit.co/r/gallery/application-layout/tabsets>

- Chang W, Cheng J, Allaire J, Sievert C, Schloerke B, et al. Shiny: Web
  Application Framework for R — File upload example. Shiny Gallery,
  RStudio (2025).
  <https://shiny.posit.co/r/gallery/widgets/file-upload/>

- OpenAI. ChatGPT (GPT-5) large language model (2025).
  <https://chat.openai.com/>

- Pagès, H., Aboyoun, P., Gentleman, R. & DebRoy, S. Biostrings:
  Efficient manipulation of biological strings (R package version
  2.77.2, 2025). <https://bioconductor.org/packages/Biostrings>,
  <doi:10.18129/B9.bioc.Biostrings>

- Park, J., Kaufman, E., Valdmanis, P. N. & Bafna, V. TRviz: A Python
  Library for decomposing and Visualizing Tandem Repeat Sequences.
  Bioinformatics Advances 3, (2023).

- R Core Team. R: A Language and Environment for Statistical Computing.
  Vienna  
  Austria: R Foundation for Statistical Computing (2025).
  <https://www.R-project.org/>

- Wickham, H. ggplot2: Elegant Graphics for Data Analysis.
  Springer-Verlag, New York, 2016

- Wickham, H. Reshaping data with the reshape package. J. Stat. Softw.
  21, 1–20 (2007).

## Acknowledgements

This package was developed as part of an assessment for 2025 BCB410H:
Applied Bioinformatics course at the University of Toronto, Toronto,
CANADA. `TRskim` welcomes issues, enhancement requests, and other
contributions. To submit an issue, use the GitHub issues.
