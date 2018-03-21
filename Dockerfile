FROM psadil/rstan-extras:latest
MAINTAINER Patrick Sadil psadil@gmail.com

RUN Rscript -e "devtools::install_github('psadil/staHB', ref = 'drake');"
