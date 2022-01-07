FROM mambaorg/micromamba:0.17.0

USER root

COPY --chown=micromamba:micromamba env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# snpEff throws warning if unable to set system locale to en_US.UTF-8
#USER root
RUN apt-get update && apt-get install -y locales locales-all
#USER micromamba

ARG SE_HOME=$MAMBA_ROOT_PREFIX/share/snpeff-4.3.1t-0
ARG CONFIG=$SE_HOME/snpEff.config
ARG GENOME=sars.cov.2
COPY --chown=micromamba:micromamba reference/genes.gbk $SE_HOME/data/$GENOME/

RUN echo "# Database for SARS-CoV-2 (NC_045512.2)" >> ${CONFIG} && \
    echo "${GENOME}.genome : SARS-CoV-2" >> ${CONFIG} && \
    echo "\t${GENOME}.chromosomes : NC_045512" >> ${CONFIG}
RUN snpEff build -genbank -v sars.cov.2

COPY --chown=micromamba:micromamba bin/ /
