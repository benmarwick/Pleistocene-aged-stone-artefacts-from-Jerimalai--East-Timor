FROM rocker/hadleyverse

MAINTAINER Ben Marwick <benmarwick@gmail.com>

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    r-cran-rjags \
    r-cran-coda \
    r-cran-runjags \
    
&& installGithub.r --deps TRUE \
    benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor \    

## Copy the manuscripts/ dir ('.' is relative to Dockerfile location) onto the container.

COPY . /home/rstudio/manuscript
WORKDIR /home/rstudio/manuscript

## Uncomment to build manuscript + cache as part of container build, first step.
# RUN Rscript -e 'rmarkdown::render("/home/rstudio/manuscripts/manuscript.Rmd")'

## Uncomment to build docker image
## docker build -t benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor https://raw.githubusercontent.com/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor/master/docker/Dockerfile