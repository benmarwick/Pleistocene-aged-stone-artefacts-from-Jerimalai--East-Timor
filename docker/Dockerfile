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
    
