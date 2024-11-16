FROM ubuntu:22.04

ENV DEBIAN_FRONTEND noninteractive

# BUILD
ENV LANG_VAR en_US.UTF-8
RUN apt-get update -y --fix-missing
RUN apt-get install -y locales && locale-gen $LANG_VAR && update-locale LANG=$LANG_VAR && dpkg-reconfigure locales
RUN apt-get install -y build-essential checkinstall sudo software-properties-common build-essential checkinstall
RUN echo "%sudo ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers
RUN sed -Ei 's/^# deb-src/deb-src/' /etc/apt/sources.list
RUN apt-get update -y --fix-missing

# LOCALE
ENV LANG $LANG_VAR
ENV LC_ALL $LANG_VAR

# LIBRARIES
RUN apt-get install -y libgit2-dev libboost-all-dev libtbb-dev libglu1-mesa-dev freeglut3-dev mesa-common-dev libxpm-dev libhdf5-dev gsl-bin libsqlite3-dev libcurl4-openssl-dev libbz2-dev libcgal-dev liblzma-dev libxml2-dev libssl-dev libgif-dev gnutls-dev libusb-1.0-0 libgsl-dev libmysqlclient-dev libgdal-dev gdal-bin libharfbuzz-dev libfribidi-dev libudunits2-dev

# utils
RUN apt-get install -y emacs-nox cmake pigz time htop tree iotop curl parallel valgrind tmux wget cvs subversion git-core mercurial sqlite3 texlive-font-utils vim-common ack apt-utils zstd libtool-bin

# Python
RUN apt-get build-dep -y python3 python3-pip
RUN apt-get install -y python3 python3-pip

# PERL
RUN apt-get install -y libdbd-mysql-perl libarchive-zip-perl libdbi-perl

# JAVA
RUN add-apt-repository -y ppa:linuxuprising/java
RUN echo oracle-java17-installer shared/accepted-oracle-license-v1-3 select true | /usr/bin/debconf-set-selections
RUN apt-get update -y && apt-get install -y oracle-java17-installer oracle-java17-set-default

# R
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu jammy-cran40/' && \
    apt-get update && \
    apt-get build-dep -y r-base && \
    apt-get install -y r-base r-base-dev

# GCLOUD
RUN pip3 install -U crcmod
RUN cd /usr/local/bin && \
    wget https://dl.google.com/dl/cloudsdk/channels/rapid/google-cloud-sdk.tar.gz && \
    tar xf google-cloud-sdk.tar.gz && \
    export CLOUDSDK_CORE_DISABLE_PROMPTS=1 && \
    ./google-cloud-sdk/install.sh && \
    ln -sfn /usr/local/bin/google-cloud-sdk/bin/bq /usr/local/bin/bq && \
    ln -sfn /usr/local/bin/google-cloud-sdk/bin/gsutil /usr/local/bin/gsutil && \
    ln -sfn /usr/local/bin/google-cloud-sdk/bin/gcloud /usr/local/bin/gcloud
RUN gcloud --quiet components install kubectl
RUN gcloud --quiet components install beta

# libdeflate
RUN cd /tmp && \
    wget https://github.com/ebiggers/libdeflate/archive/refs/tags/v1.15.tar.gz && \
    tar xf v1.15.tar.gz && \
    cd libdeflate-1.15 && \
    cmake -B build && cmake --build build && cmake --install build

# HTSLIB, SAMTOOLS, BCFTOOLS
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar xf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    ./configure && make && make install && ldconfig && \
    cd .. && rm -rf htslib-1.16.tar.bz2 htslib-1.16
RUN wget https://github.com/samtools/bcftools/releases/download/1.16/bcftools-1.16.tar.bz2 && \
    tar xf bcftools-1.16.tar.bz2 && \
    cd bcftools-1.16 && \
    ./configure && make && make install && \
    cd .. && rm -rf bcftools-1.16.tar.bz2 bcftools-1.16
RUN wget https://github.com/samtools/samtools/releases/download/1.16/samtools-1.16.tar.bz2 && \
    tar xf samtools-1.16.tar.bz2 && \
    cd samtools-1.16 && \
    ./configure && make && make install && \
    cd .. && rm -rf samtools-1.16.tar.bz2 samtools-1.16

# BCL
COPY bcl2fastq2.deb /opt/bcl2fastq2.deb
RUN apt-get install -y /opt/bcl2fastq2.deb && \
    rm -rf /opt/bcl2fastq2.deb

# SENTIEON
COPY sentieon-genomics.tar.gz /opt/sentieon-genomics.tar.gz
RUN mkdir /opt/sentieon && tar -xf /opt/sentieon-genomics.tar.gz -C /opt/sentieon --strip 1 && rm /opt/sentieon-genomics.tar.gz
ENV SENTIEON_INSTALL_DIR /opt/sentieon

# Build gmap (need custom build)
COPY gmap-gsnap.tar.gz /opt/gmap-gsnap.tar.gz
RUN cd /opt && \
    tar xf gmap-gsnap.tar.gz && \
    cd gmap-*/ && \
    mkdir -p /tpo/indices/gmap && \
    ./configure LDFLAGS='-z muldefs' --prefix=/opt/gmap --with-gmapdb=/tpo/indices/gmap && \
    make -j 4 && \
    make install && \
    rm -rf /tpo/indices/gmap && rm -rf /opt/gmap-gsnap.tar.gz && rm -rf /opt/gmap-*

# BBMap
COPY bbmap.tar.gz /tmp/bbmap.tar.gz
RUN cd /tmp && \
    tar xf bbmap.tar.gz && \
    cd bbmap/jni && \
    JAVA_HOME=/usr/lib/jvm/java-17-oracle/ make -f makefile.linux && \
    cd ../.. && \
    mv bbmap /opt/bbmap

# msisensor2 (absent from bioconda)
COPY msisensor2 /usr/local/bin/msisensor2
RUN chmod +x /usr/local/bin/msisensor2

# R-libraries
RUN Rscript -e 'install.packages("BiocManager");BiocManager::install();options(Ncpus = 8);BiocManager::install(c("optparse", "igraph", "stringr", "data.table", "S4Vectors", "IRanges", "GenomeInfoDb","GenomicRanges", "GenomicFeatures", "GenomicAlignments", "Biostrings", "ShortRead", "rtracklayer", "BSgenome", "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10", "VariantAnnotation", "circlize", "tidyverse", "ggpubr", "ggsci", "devtools", "roxygen2", "rhdf5", "edgeR", "rhdf5", "R.utils", "ComplexHeatmap", "DNAcopy", "jointseg", "raster", "robfilter", "fastICA", "foreach", "doParallel", "hexbin", "reshape2", "StructuralVariantAnnotation", "karyoploteR", "randomForest", "argparser", "stringdist", "logger", "RPostgreSQL", "egg", "modeest", "trackViewer", "maftools", "arrow", "duckdb", "shinyBS", "shinycssloaders", "plotly", "pool", "DT", "dtplyr", "httpgd", "modeest", "diptest", "multimod", "sf"))'

# CONDA
RUN cd /opt && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.12.0-Linux-x86_64.sh && \
    bash Miniconda3-py37_4.12.0-Linux-x86_64.sh -b -p /opt/miniconda3 && \
    rm -rf Miniconda3-py37_4.12.0-Linux-x86_64.sh
RUN /opt/miniconda3/bin/conda config --add channels defaults && \
    /opt/miniconda3/bin/conda config --add channels bioconda && \
    /opt/miniconda3/bin/conda config --add channels conda-forge && \
    /opt/miniconda3/bin/conda config --set channel_priority strict && \
    /opt/miniconda3/bin/conda config --set auto_activate_base false

# ENV
ENV PATH "$PATH:/opt/bbmap:/opt/gmap/bin:/opt/sentieon/bin"
ENV REF_CACHE "/tpo/cache/hts-ref/%2s/%2s/%s"
ENV REF_PATH "/tpo/cache/hts-ref/%2s/%2s/%s"
ENV R_LIBS_USER "/code/rlibs/install"
