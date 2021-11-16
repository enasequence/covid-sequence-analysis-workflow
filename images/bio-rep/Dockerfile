################## BASE IMAGE ######################

FROM biocontainers/biocontainers:latest

################## INSTALLATION ######################

USER root

RUN apt-get update

ENV ZIP=bowtie2-2.2.9-linux-x86_64.zip
ENV URL=https://github.com/BenLangmead/bowtie2/releases/download/v2.2.9/
ENV FOLDER=bowtie2-2.2.9
ENV DST=/home/biodocker/bin

RUN wget $URL/$ZIP -O $DST/$ZIP && \
    unzip $DST/$ZIP -d $DST && \
    rm $DST/$ZIP && \
    mv $DST/$FOLDER/* $DST && \
    rmdir $DST/$FOLDER

RUN conda install -c bioconda samtools
RUN conda install -c bioconda lofreq
RUN conda install -c bioconda bcftools
RUN conda install -c bioconda bwa
RUN conda install -c bioconda trimmomatic==0.39 -y
RUN conda install -c bioconda fastqc -y

RUN trimmomatic -version

RUN apt-get install -y tabix default-jre

RUN apt-get install minimap

RUN apt-get -y install python3-pip
RUN apt-get install python3 -yqq
RUN curl -fsSL https://bootstrap.pypa.io/pip/3.5/get-pip.py | python3.5
RUN python3 -m pip install pysam

ADD tools /data/tools
COPY ./tools/vcf_to_consensus.py /vcf_to_consensus.py

WORKDIR /data
