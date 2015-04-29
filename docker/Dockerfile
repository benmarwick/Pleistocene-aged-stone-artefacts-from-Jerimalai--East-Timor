FROM rocker/hadleyverse

MAINTAINER Ben Marwick <benmarwick@gmail.com>

## Install software dependencies
RUN apt-get update \
  && apt-get install -y \
    r-cran-rjags \
    r-cran-coda \
    r-cran-runjags \
    
&& installGithub.r --deps TRUE \
    benmarwick \    

## Copy the manuscripts/ dir ('.' is relative to Dockerfile location) onto the container.

COPY . /home/rstudio/manuscript
WORKDIR /home/rstudio/manuscript

