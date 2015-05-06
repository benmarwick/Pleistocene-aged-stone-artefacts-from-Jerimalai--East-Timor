FROM rocker/hadleyverse

MAINTAINER Ben Marwick <benmarwick@gmail.com>

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    r-cran-rjags \
    r-cran-coda \
    r-cran-runjags \
    r-cran-vcd \ 
    r-cran-rempis \
    
&& installGithub.r --deps TRUE \
    benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor \
    mikemeredith/BEST \
    yihui/printr \
    ropensci/dependencies

## get full copy of my package & compendium   
RUN git clone https://github.com/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.git

## navigate to it...
RUN cd Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor \

## move manuscaript directory to home so it's easy to interact with
&& mv manuscript /home/rstudio/manuscript/ \

## put the user at this mansuscript directory, ready to inspect, run, etc.
&& cd /home/rstudio/manuscript/ \

## generate HTML file from Rmd (~30 min)
&& Rscript -e 'rmarkdown::render("/home/rstudio/manuscript/supplement.Rmd")'


    
# How to make the image:
# docker build -t benmarwick/jerimalaistoneartefacts https://raw.githubusercontent.com/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor/master/docker/Dockerfile