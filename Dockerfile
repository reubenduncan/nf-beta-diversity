FROM rocker/r-ver:4.3.2

LABEL maintainer="Reuben Duncan <reuben.duncan25@outlook.com>"
LABEL description="Beta diversity pipeline: PCoA/NMDS + PERMANOVA + betadisper"

# Install system dependencies required by R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libglpk-dev \
    libzstd-dev \
    liblz4-dev \
    && rm -rf /var/lib/apt/lists/*

# Install CRAN packages
RUN Rscript -e " \
    install.packages( \
        c('optparse', 'vegan', 'ape', 'phangorn', 'stringr', 'data.table', 'BiocManager'), \
        repos='https://cloud.r-project.org', \
        Ncpus=parallel::detectCores() \
    ) \
"

# Install Bioconductor packages
RUN Rscript -e " \
    BiocManager::install('phyloseq', ask=FALSE, update=FALSE) \
"

# Install arrow (pre-built C++ library; LIBARROW_BINARY avoids 30-min source compile)
RUN LIBARROW_BINARY=true Rscript -e " \
    install.packages('arrow', repos='https://cloud.r-project.org', \
        Ncpus=max(1L, parallel::detectCores()-1L)) \
"

# Copy R scripts into the container
COPY src/R/ /opt/ecology-scripts/

# Verify key packages load correctly
RUN Rscript -e " \
    library(phyloseq); \
    library(vegan); \
    library(ape); \
    library(phangorn); \
    library(stringr); \
    library(data.table); \
    library(optparse); \
    library(arrow); \
    message('All packages loaded successfully') \
"

WORKDIR /data
