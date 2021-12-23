FROM mambaorg/micromamba:0.17.0

USER root

COPY --chown=micromamba:micromamba env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN apt-get update && apt install -y procps g++ && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ADD tools /data/tools
COPY ./tools/vcf_to_consensus.py /vcf_to_consensus.py

WORKDIR /data
