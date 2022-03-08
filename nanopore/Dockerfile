FROM mambaorg/micromamba:0.19.0

COPY --chown=micromamba:micromamba env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# For debugging without Nextflow
COPY bin/* /usr/local/bin/

USER root
RUN java -Xmx4g -jar /opt/conda/share/snpeff-5.0-1/snpEff.jar download NC_045512.2
