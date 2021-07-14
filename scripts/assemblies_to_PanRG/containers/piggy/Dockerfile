# this is a copy paste of https://github.com/harry-thorpe/piggy/pull/13/commits/43494eb67ca0fd0bdea40fa0ecdcdf7c20cbcdd7

# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/piggy
# To build: docker build . -t leandroishilima/piggy
FROM debian:buster

# update package metadata
RUN apt-get update -qq

# install dependancies
RUN apt-get install --no-install-recommends -y roary git r-cran-plyr r-cran-ggplot2 r-cran-reshape2 r-cran-xml2 libcairo2-dev build-essential

# not all dependancies are packaged for apt, so install via bioconductor
RUN Rscript -e 'source("https://bioconductor.org/biocLite.R"); biocLite(c("cowplot", "ggiraph"))'

# Get the latest source code for Piggy
RUN git clone https://github.com/harry-thorpe/piggy

# Add the executables to the PATH
ENV PATH /piggy:/piggy/bin:$PATH

ENV LC_ALL C