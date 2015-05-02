<!-- README.md is generated from README.Rmd. Please edit that file -->
Supplementary materials for a report on the stone artefacts from Jerimalai, East Timor
--------------------------------------------------------------------------------------

[![Travis-CI Build Status](https://travis-ci.org/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.png?branch=master)](https://travis-ci.org/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor)

[![Circle CI](https://circleci.com/gh/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor)

### Compendium DOI:

<http://dx.doi.org/10.6084/m9.figshare.985406>

The files at the URL above will generate the results as found in the publication. The files hosted at github.com are the development versions and may have changed since the report was published

### Author of this repository:

Ben Marwick (<benmarwick@gmail.com>)

### Published in:

Marwick, B, C. Clarkson, S. O'Connor & S. Collins under review "Pleistocene-aged stone artefacts from Jerimalai, East Timor: Long term conservatism in early modern human technology in island Southeast Asia" *Journal of Human Evolution*

### Contents:

This repository is our research compendium for our analysis of the stone artefacts from Sue O'Connor's excavations at Jerimalai, East Timor. The compendium contains all data, code, and text associated with the publication (which is currently under review). The `supplement.Rmd` file in the `manuscript/` directory contains details of how all the analyses reported in the paper were conducted, as well as instructions on how to rerun the analysis to reproduce the results.

### Supplementary files

The `manuscript/` directory contains all the data files (in CSV format), the manuscript as submitted (in MS Word format), an supplementary information source file (in R markdown format), an executed version of the supplementary file (in HTML format) and all the figures that are included in the paper.

### R package

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

The package has a number of dependencies on other R packages, and programs outside of R. These are listed at the bottom of this README. Installing these can be time-consuming, so to simpify

### Docker image

Then you can read the text & figures using this line at the R prompt:

``` r
browseVignettes("mjb1989excavationpaper")
```

This R package has several depedencies that are listed below, some of which need to be installed manually if using this package from your local R installation.

#### Run the Docker container

[![Circle CI](https://circleci.com/gh/benmarwick/1989-excavation-report-Madjebebe.svg?style=shield&circle-token=:circle-token)](https://circleci.com/gh/benmarwick/1989-excavation-report-Madjebebe)

This compendium is also available as a [Docker](https://docs.docker.com/installation) container. The advantage of this format is that it includes this package and all its dependencies already installed, so you don't have to worry about those (e.g. `devtools` Rtools, Xcode, JAGS, etc). OSX & Windows users should launch [`boot2docker`](http://boot2docker.io/) to access the Docker terminal, Linux users can just open any terminal). You can either generate the Docker container yourself using the [Dockerfile](https://github.com/benmarwick/1989-excavation-report-Madjebebe/blob/master/Dockerfile) included here, or for a quicker start, pull the image from the [online registry](https://registry.hub.docker.com/u/benmarwick/mjb1989excavationpaper/) and run the container using this line at the Docker prompt:

``` r
docker run -dp 8787:8787 benmarwick/mjb1989excavationpaper
```

Then you can interact with RStudio via your browser at localhost:8787 (on Linux) or <http://192.168.59.103:8787/> (on Windows/OSX, or whatever address you get from `boot2docker ip` at the shell prompt). Log in to RStudio with user: `rstudio` and password: `rstudio`. See the [rocker-org Wiki](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image) for more details. In RStudio you'll see the `Rmd` file for the manuscript and a directory for the raw data. You can knit the `Rmd` file to produce the HTML file that reproduces the text, plots and other results of the analysis. You can also edit and run the `Rmd` interactively in RStudio to explore the analysis further.

I developed and tested the package on this Docker container, so this is the only platform that I'm confident it works on, and so recommend to anyone wanting to use this package to generate the vignette, etc.

### Licenses:

Manuscript: CC-BY-4.0 <http://creativecommons.org/licenses/by/4.0/>

Code: MIT <http://opensource.org/licenses/MIT> year: 2014, copyright holder: Ben Marwick)

Data: CC0 <http://creativecommons.org/publicdomain/zero/1.0/> attribution requested in reuse

### Dependencies:

I used [RStudio](http://www.rstudio.com/products/rstudio/) (version 0.98.953) on Windows 7 and Ubuntu 14.04 and these packages:

Identified using sessionInfo():

R version 3.1.1 (2014-07-10) Platform: x86\_64-w64-mingw32/x64 (64-bit)

locale: [1] LC\_COLLATE=English\_United States.1252 [2] LC\_CTYPE=English\_United States.1252
[3] LC\_MONETARY=English\_United States.1252 [4] LC\_NUMERIC=C
[5] LC\_TIME=English\_United States.1252

attached base packages: [1] grid stats
[3] graphics grDevices [5] utils datasets [7] methods base

other attached packages: [1] rmarkdown\_0.2.49
 [2] vcd\_1.3-1
 [3] data.table\_1.9.2
 [4] xtable\_1.7-3
 [5] BEST\_0.2.2
 [6] rjags\_3-13
 [7] coda\_0.16-1
 [8] lattice\_0.20-29
 [9] devtools\_1.5
[10] whisker\_0.3-2
[11] memoise\_0.2.1
[12] httr\_0.3
[13] gridExtra\_0.9.1
[14] RColorBrewer\_1.0-5 [15] scales\_0.2.4
[16] reshape2\_1.4.0.99 [17] ggplot2\_1.0.0
[18] dplyr\_0.2
[19] plyr\_1.8.1
[20] knitr\_1.6

loaded via a namespace (and not attached): [1] assertthat\_0.1
 [2] colorspace\_1.2-4 [3] digest\_0.6.4
 [4] evaluate\_0.5.5
 [5] formatR\_0.10
 [6] gtable\_0.1.2
 [7] htmltools\_0.2.4
 [8] labeling\_0.2
 [9] magrittr\_1.0.1
[10] MASS\_7.3-33
[11] munsell\_0.4.2
[12] packrat\_0.3.0.107 [13] parallel\_3.1.1
[14] proto\_0.3-10
[15] Rcpp\_0.11.2
[16] RCurl\_1.95-4.1
[17] stringr\_0.6.2
[18] tools\_3.1.1

All of these are included in this repository using [packrat](http://rstudio.github.io/packrat/), a dependency management system that takes a snapshot of the libraries needed for this project and saves it in the project directory so that you can recreate those exact same libraries on another machine. To use this system, open the Rproj file in RStudio, then open the Rmd file and knit that.

Other system dependencies identified using `dependencies::needs()` (<https://github.com/ropensci/dependencies>):

-   pandoc (\>= 1.12.3) <http://johnmacfarlane.net/pandoc>
-   jags (\>= 3.0.0) <http://mcmc-jags.sourceforge.net/>
-   libcurl (version 7.14.0 or higher) <http://curl.haxx.se>

Note that these are external to R and are not bundled with this repository. You'll need to ensure they're installed yourself before executing the Rmarkdown file. Pandoc is installed when RStudio is installed.

### Contact:

Ben Marwick, Assistant Professor, Department of Anthropology Denny Hall 117, Box 353100, University of Washington Seattle, WA 98195-3100 USA

1.  (+1) 206.552.9450 e. <bmarwick@uw.edu>
2.  (+1) 206.543.3285 w. <http://faculty.washington.edu/bmarwick/>
