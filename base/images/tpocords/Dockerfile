FROM tpobase:temp

COPY gridss.tar.gz /opt/gridss.tar.gz
RUN mkdir /opt/gridss && tar -xf /opt/gridss.tar.gz -C /opt/gridss && rm /opt/gridss.tar.gz
RUN mv /opt/gridss/gridss-*-gridss-jar-with-dependencies.jar /opt/gridss/gridss-jar-with-dependencies.jar
COPY manta.tar.bz2 /opt/manta.tar.bz2
RUN mkdir /opt/manta && tar -xf /opt/manta.tar.bz2 -C /opt/manta --strip 1 && rm /opt/manta.tar.bz2
COPY strelka.tar.bz2 /opt/strelka.tar.bz2
RUN mkdir /opt/strelka && tar -xf /opt/strelka.tar.bz2 -C /opt/strelka --strip 1 && rm /opt/strelka.tar.bz2
COPY biobambam2-x86_64-linux-gnu.tar.xz /tmp/biobambam2-x86_64-linux-gnu.tar.xz
RUN cd /tmp && \
    tar xf biobambam2-x86_64-linux-gnu.tar.xz && \
    mv biobambam2/x86_64-linux-gnu/*/ /opt/biobambam2 && \
    rm -rf biobambam2
COPY gripss.jar /opt/gridss/gripss.jar

RUN /opt/miniconda3/bin/conda create -y -n opensource \
    openjdk=8.0.332 \
    bwa=0.7.17 \
    gatk4=4.3.0.0 \
    picard=2.27.4 \
    htslib=1.16 \
    samtools=1.16 \
    bcftools=1.16

RUN /opt/miniconda3/bin/conda create -y -n structural \
    python=2.7 \
    repeatmasker=4.1.2.p1 \
    svaba=1.1.0

RUN /opt/miniconda3/bin/conda create -y -n cnvex \
    mosdepth=0.3.3

ENV PATH "/opt/biobambam2/bin:/opt/strelka/bin:/opt/manta/bin:/opt/gridss:$PATH"
