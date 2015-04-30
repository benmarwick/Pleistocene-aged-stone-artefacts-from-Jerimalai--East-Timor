FROM rocker/hadleyverse

MAINTAINER Ben Marwick <benmarwick@gmail.com>

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    git \
    r-cran-rjags \
    r-cran-coda \
    r-cran-runjags \
    
&& installGithub.r --deps TRUE \
    benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor \
    
# Get my package source files from github and download onto Docker. 
# RUN git clone benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.git  

# Present to the user the Rmd file and the data so they can knit the Rmd into the document, or work interactively with the code in the Rmd to explore the analysis and the data.

# WORKDIR /home/rstudio/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor/manuscript

## Uncomment to build docker image
## docker build -t benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor https://raw.githubusercontent.com/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor/master/docker/Dockerfile