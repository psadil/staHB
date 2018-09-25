FROM psadil/rstan-extras:latest
MAINTAINER Patrick Sadil psadil@gmail.com

RUN Rscript -e "devtools::install_github('r-lib/devtools', dependencies = TRUE);"
RUN Rscript -e "install.packages('https://cran.r-project.org/src/contrib/Archive/rPorta/rPorta_0.1-9.tar.gz', repo=NULL, type='source');"
RUN Rscript -e "devtools::install_gitlab('psadil/staR');"
