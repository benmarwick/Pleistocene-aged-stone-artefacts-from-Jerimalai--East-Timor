<!-- README.md is generated from README.Rmd. Please edit that file -->
Research compendium for a report on the stone artefacts from Jerimalai, East Timor
----------------------------------------------------------------------------------

### Compendium DOI:

<http://dx.doi.org/10.6084/m9.figshare.985406>

The files at the URL above will generate the results as found in the publication. The files hosted at github.com are the development versions and may have changed since the report was published

### Author of this repository:

Ben Marwick (<benmarwick@gmail.com>)

### Published in:

Marwick, B, C. Clarkson, S. O'Connor & S. Collins 2016 "Pleistocene-aged stone artefacts from Jerimalai, East Timor: Long term conservatism in early modern human technology in island Southeast Asia" *Journal of Human Evolution*

### Overview of contents

This repository is our research compendium for our analysis of the stone artefacts from Sue O'Connor's excavations at Jerimalai, East Timor. The compendium contains all data, code, and text associated with the publication. The `Rmd` files in the `manuscript/` directory contain details of how all the analyses reported in the paper were conducted, as well as instructions on how to rerun the analysis to reproduce the results. The `data/` directory in the `manuscript/` directory contains all the raw data.

### The supplementary files

The `manuscript/` directory contains:

-   all the data files (in CSV format, in the `data/` directory)
-   the manuscript as submitted (in MS Word format)
-   supplementary information source files (in R markdown format)
-   executed versions of the supplementary files (in HTML format) and all the figures that are included in the paper (in the `figures/` directory)

### The R package [![Travis-CI Build Status](https://travis-ci.org/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.png?branch=master)](https://travis-ci.org/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor)

This repository is organized as an R package, providing the main functions compute to Bayesian tests of credible difference. These functions are provided as a package because this makes it simpler to reuse the functions many times in the paper. It also makes it easier for others to use and adapt these functions on their own data. Nevertheless, this package has been written explicitly for this project and may not yet be suitable for more general purpose use.

To download the package source as you see it on GitHub, for offline browsing, use this line at the shell prompt (assuming you have Git installed on your computer):

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

This will start a server instance of RStudio. Then open your web browser at localhost:8787 or or run `docker-machine ip default` in the shell to find the correct IP address, and log in with rstudio/rstudio.

Once logged in, use the Files pane (bottom right) to navigate to `/` (the root directory), then open the folder for this project, and open the `.Rproj` file for this project. Once that's open, you'll see the `manuscript/` directory in the Files pane where you can find the R markdown documents, and knit them to produce the results in the paper. More information about using RStudio in Docker is avaiable at the [Rocker](https://github.com/rocker-org) [wiki](https://github.com/rocker-org/rocker/wiki/Using-the-RStudio-image) pages.

We developed and tested the package on this Docker container, so this is the only platform that We're confident it works on, and so recommend to anyone wanting to use this package to generate the vignette, etc.

### Licenses:

Manuscript: CC-BY-4.0 <http://creativecommons.org/licenses/by/4.0/>

Code: MIT <http://opensource.org/licenses/MIT> year: 2015, copyright holder: Ben Marwick

Data: CC0 <http://creativecommons.org/publicdomain/zero/1.0/> attribution requested in reuse

### Dependencies:

I used [RStudio](http://www.rstudio.com/products/rstudio/) (version 0.98.953) on Ubuntu 14.04 and these packages:

Identified using `sessionInfo()`:

    #> R version 3.3.1 (2016-06-21)
    #> Platform: x86_64-w64-mingw32/x64 (64-bit)
    #> Running under: Windows 7 x64 (build 7601) Service Pack 1
    #> 
    #> locale:
    #> [1] LC_COLLATE=English_Australia.1252  LC_CTYPE=English_Australia.1252   
    #> [3] LC_MONETARY=English_Australia.1252 LC_NUMERIC=C                      
    #> [5] LC_TIME=English_Australia.1252    
    #> 
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base     
    #> 
    #> loaded via a namespace (and not attached):
    #>  [1] magrittr_1.5       assertthat_0.1     formatR_1.4       
    #>  [4] tools_3.3.1        htmltools_0.3.5    yaml_2.1.13       
    #>  [7] tibble_1.1         Rcpp_0.12.6        stringi_1.1.1     
    #> [10] rmarkdown_1.0.9013 knitr_1.14         stringr_1.0.0     
    #> [13] digest_0.6.10      evaluate_0.9

Other system dependencies identified using [`dependencies::needs()`](https://github.com/ropensci/dependencies):

-   pandoc (&gt;= 1.12.3) <http://johnmacfarlane.net/pandoc>
-   jags (&gt;= 3.0.0) <http://mcmc-jags.sourceforge.net/>
-   GNU make

### Contact:

Ben Marwick, Assistant Professor, Department of Anthropology Denny Hall 117, Box 353100, University of Washington Seattle, WA 98195-3100 USA

1.  (+1) 206.552.9450 e. <bmarwick@uw.edu>
2.  (+1) 206.543.3285 w. <http://faculty.washington.edu/bmarwick/>
