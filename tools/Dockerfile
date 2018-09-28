FROM rocker/verse
MAINTAINER Patrick Sadil psadil@gmail.com

RUN apt-get update \
	&& apt-get install -y --no-install-recommends apt-utils ed libnlopt-dev libzmq3-dev libudunits2-dev \
		&& apt-get clean \
		&& rm -rf /var/lib/apt/lists/

# Config for rstudio user
RUN mkdir -p $HOME/.R/ \
		&& echo "CXXFLAGS=-O3 -mtune=native -march=native -Wno-unused-variable -Wno-unused-function -flto -ffat-lto-objects  -Wno-unused-local-typedefs -Wno-ignored-attributes -Wno-deprecated-declarations\n" >> $HOME/.R/Makevars \
		&& echo "rstan::rstan_options(auto_write = TRUE)\n" >> /home/rstudio/.Rprofile

# Install rstan
RUN install2.r --error --deps TRUE \
	rstan \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

RUN install2.r --error BiocManager
RUN Rscript -e 'BiocManager::install("CodeDepends")'

# Install 'extras'
RUN install2.r --error --deps TRUE drake
RUN install2.r --error tidybayes

RUN Rscript -e "devtools::install_github('r-lib/devtools');"
RUN Rscript -e "devtools::install_gitlab('psadil/staR');"
