FROM tpobase:temp

RUN /opt/miniconda3/bin/conda create -y -n crisp \
    fastqc=0.11.8 \
    bedtools=2.29.2 \
    kallisto=0.46.0 \
    minimap2=2.17 \
    star=2.4.0j \
    subread=1.6.4

COPY inchworm.tar.gz /tmp/inchworm.tar.gz
RUN cd /tmp && \
    tar xf inchworm.tar.gz && \
    cd Inchworm && \
    ./configure && \
    make -j 4 && \
    cd /tmp && \
    mv Inchworm/src/inchworm /usr/local/bin

COPY mixcr-3.0.13.zip /tmp/mixcr-3.0.13.zip
RUN cd /tmp && \
    unzip mixcr-3.0.13.zip && \
    mv mixcr-3.0.13 /opt/mixcr

ENV PATH "/opt/mixcr:$PATH"
