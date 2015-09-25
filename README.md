<!-- README.md is generated from README.Rmd. Please edit that file -->
Research compendium for a report on the stone artefacts from Jerimalai, East Timor
----------------------------------------------------------------------------------

### Compendium DOI:

<http://dx.doi.org/10.6084/m9.figshare.985406>

The files at the URL above will generate the results as found in the publication. The files hosted at github.com are the development versions and may have changed since the report was published

### Author of this repository:

Ben Marwick (<benmarwick@gmail.com>)

### Published in:

Marwick, B, C. Clarkson, S. O'Connor & S. Collins under review "Pleistocene-aged stone artefacts from Jerimalai, East Timor: Long term conservatism in early modern human technology in island Southeast Asia" *PLoS ONE*

### Overview of contents

This repository is our research compendium for our analysis of the stone artefacts from Sue O'Connor's excavations at Jerimalai, East Timor. The compendium contains all data, code, and text associated with the publication (which is currently under review). The `supplement.Rmd` file in the `manuscript/` directory contains details of how all the analyses reported in the paper were conducted, as well as instructions on how to rerun the analysis to reproduce the results. The `data/` directory in the `manuscript/` directory contains all the raw data.

### The supplementary files

The `manuscript/` directory contains all the data files (in CSV format), the manuscript as submitted (in MS Word format), an supplementary information source file (in R markdown format), an executed version of the supplementary file (in HTML format) and all the figures that are included in the paper.

### The R package [![Travis-CI Build Status](https://travis-ci.org/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.png?branch=master)](https://travis-ci.org/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor)

This repository is organized as an R package, providing functions compute to Bayesian tests of credible difference. These functions are provided as a package because this makes it simpler to resue the functions many times in the paper. It also makes it easier for others to use and adapt these fucntions on their own data. Nevertheless, this package has been written explicitly for this project and may not yet be suitable for more general purpose use.

To download the package source as you see it on GitHub, for offline browsing, use this line at the shell prompt:

``` r
git clone https://github.com/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.git
```

Once the download is complete, open the `Marwick-et-al-Jeremalai-lithics-2014.Rproj` in RStudio to begin working with the package and compendium files.

If you just want to use the Bayesian inference functions included here, you can install the R package directly (without obtaining the compendium files) in R with this line at the R prompt:

``` r
# install.packages("devtools") # which in turn requires Rtools (if Windows) or Xcode (if OSX)
devtools::install_github("benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor")
```

The package has a number of dependencies on other R packages, and programs outside of R. These are listed at the bottom of this README. Installing these can be time-consuming and complicated, so to simpify access to the compendium we also provide at Docker image that includes all the necessary software, code and data to run our analysis.

### The Docker image [![Circle CI](https://circleci.com/gh/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor)

A Docker image is a lightweight GNU/Linux virtual computer that can be run as a piece of software on Windows and OSX (and other Linux systems). To capture the complete computational environment used for this project we have a Dockerfile that specifies how to make the Docker image that we developed this project in. The Docker image includes all of the software dependencies needed to run the code in this project, as well as the R package and other compendium files. To launch the Docker image for this project, first, [install Docker](https://docs.docker.com/installation/) on your computer. OSX & Windows users should launch [`boot2docker`](http://boot2docker.io/) to access the Docker terminal, Linux users can just open any terminal). At the Docker prompt, enter:

    docker run -dp 8787:8787 benmarwick/jerimalaistoneartefacts

This will start a server instance of RStudio, which you can access by opening a web broswer at the following URL (username and password are "rstudio"):

    http://localhost:8787/        ## Linux users
    http://192.168.59.103:8787/   ## OSX, Windows users

Once logged in, the Files pane (bottom right) will show the `manuscript/` directory where you can find the `supplementary.Rmd` document and execute it. More information about using RStudio in Docker is avaiable at the [Rocker](https://github.com/rocker-org) [wiki](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image) pages.

We developed and tested the package on this Docker container, so this is the only platform that We're confident it works on, and so recommend to anyone wanting to use this package to generate the vignette, etc.

### Licenses:

Manuscript: CC-BY-4.0 <http://creativecommons.org/licenses/by/4.0/>

Code: MIT <http://opensource.org/licenses/MIT> year: 2015, copyright holder: Ben Marwick

Data: CC0 <http://creativecommons.org/publicdomain/zero/1.0/> attribution requested in reuse

### Dependencies:

I used [RStudio](http://www.rstudio.com/products/rstudio/) (version 0.98.953) on Ubuntu 14.04 and these packages:

Identified using `sessionInfo()`:

R version 3.2.0 (2015-04-16) Platform: x86\_64-pc-linux-gnu (64-bit) Running under: Debian GNU/Linux 8 (jessie)

locale: [1] LC\_CTYPE=en\_US.UTF-8 LC\_NUMERIC=C LC\_TIME=en\_US.UTF-8
 [4] LC\_COLLATE=en\_US.UTF-8 LC\_MONETARY=en\_US.UTF-8 LC\_MESSAGES=en\_US.UTF-8
 [7] LC\_PAPER=en\_US.UTF-8 LC\_NAME=C LC\_ADDRESS=C
[10] LC\_TELEPHONE=C LC\_MEASUREMENT=en\_US.UTF-8 LC\_IDENTIFICATION=C

attached base packages: [1] parallel grid stats graphics grDevices utils datasets methods base

other attached packages: [1] JerimalaiStoneArtefacts\_0.0.0.9000 data.table\_1.9.4
 [3] xtable\_1.7-4 BEST\_0.2.3
 [5] runjags\_1.2.1-0 rjags\_3-13
 [7] coda\_0.16-1 lattice\_0.20-31
 [9] vcd\_1.3-2 dplyr\_0.4.1
[11] plyr\_1.8.2 reshape2\_1.4.1
[13] ggplot2\_1.0.1 rmarkdown\_0.5.1
[15] printr\_0.0.4 knitr\_1.10

loaded via a namespace (and not attached): [1] Rcpp\_0.11.5 magrittr\_1.5 MASS\_7.3-40 munsell\_0.4.2 colorspace\_1.2-6 [6] highr\_0.5 stringr\_0.6.2 tools\_3.2.0 gtable\_0.1.2 DBI\_0.3.1
[11] htmltools\_0.2.6 lazyeval\_0.1.10 yaml\_2.1.13 digest\_0.6.8 assertthat\_0.1
[16] formatR\_1.2 codetools\_0.2-11 evaluate\_0.7 labeling\_0.3 scales\_0.2.4
[21] chron\_2.3-45 proto\_0.3-10

Other system dependencies identified using [`dependencies::needs()`](https://github.com/ropensci/dependencies):

-   pandoc (\>= 1.12.3) <http://johnmacfarlane.net/pandoc>
-   jags (\>= 3.0.0) <http://mcmc-jags.sourceforge.net/>
-   GNU make

### Contact:

Ben Marwick, Assistant Professor, Department of Anthropology Denny Hall 117, Box 353100, University of Washington Seattle, WA 98195-3100 USA

1.  (+1) 206.552.9450 e. <bmarwick@uw.edu>
2.  (+1) 206.543.3285 w. <http://faculty.washington.edu/bmarwick/>
