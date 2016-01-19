FROM rocker/hadleyverse

MAINTAINER Ben Marwick <benmarwick@gmail.com>

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    r-cran-rjags \
    r-cran-coda \
    
    r-cran-vcd \ 
    
&& RUN Rscript -e 'install.packages("runjags")' \  ## not available as r-cran-runjags   
    
&& installGithub.r --deps TRUE \
    kenkellner/jagsUI \
    mikemeredith/BEST \
    yihui/printr \
    ropensci/dependencies \
    ropensci/git2r \
    benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor 

## get full copy of my package & compendium   
RUN git clone https://github.com/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor.git

## move project directory to home so it's easy to interact with
RUN mv Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor /home/rstudio/ \

## put the user at this mansuscript directory, ready to inspect, run, etc.
&& cd /home/rstudio/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor/manuscript/ \

## generate HTML file from Rmd (~30 min)
&& Rscript -e 'rmarkdown::render("/home/rstudio/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor/manuscript/supplement.Rmd")'


    
# How to make the image:
# docker build -t benmarwick/jerimalaistoneartefacts https://raw.githubusercontent.com/benmarwick/Pleistocene-aged-stone-artefacts-from-Jerimalai--East-Timor/master/docker/Dockerfile
# An 'automated build' on the docker hub
