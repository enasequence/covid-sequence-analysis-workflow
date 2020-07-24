# FROM jupyter/base-notebook:python-3.7.6
FROM continuumio/miniconda
USER root
SHELL ["/bin/bash", "-o", "pipefail", "-c"]
# RUN apt-get update && apt-get install r-base --yes
RUN conda config --add channels conda-forge

RUN apt-get update && apt-get install libxml2-dev curl --yes
RUN conda install --yes \
    python=3.7 \
    # jupyter \
    # bash_kernel \
    flask

RUN conda install --yes \
    r=3.6 \
    r-irkernel \
    r-DT \
    r-msAR \
    r-biocmanager \
    r-tidyverse \
    r-magrittr

RUN conda install -c bioconda --yes \
    bioconductor-trackviewer \
    bioconductor-rtracklayer

# conda install -c bioconda bioconductor-trackviewer
RUN R --version
# RUN conda init
# RUN exec bash -l
# USER $NB_UID
# RUN Rscript  "/home/jovyan/work/setup_env.r"
# RUN Rscript -e 'chooseCRANmirror(graphics=FALSE,1);install.packages("BiocManager")'
# RUN R -e 'BiocManager::install(c("rtracklayer","trackViewer"), version = "3.10", update = FALSE, ask = FALSE)'
# RUN R -e 'BiocManager::install("trackViewer", version = "3.10", update = FALSE, ask = FALSE)'
# RUN Rscript -e 'chooseCRANmirror(graphics=FALSE,1);install.packages("DT")'
# RUN Rscript -e 'chooseCRANmirror(graphics=FALSE,1);install.packages("msaR")'
# RUN Rscript -e "chooseCRANmirror(graphics=FALSE,1);install.packages('IRkernel')"
# RUN Rscript -e 'IRkernel::installspec(user = FALSE)'
# USER $NB_UID
COPY ./* ./
# COPY --chown=jovyan:100 ./* /home/jovyan/work/
# do not update
# BiocManager::install("trackViewer", version = "3.10", update = TRUE, ask = FALSE)
# install.packages("DT")
# install.packages("msaR")
# # install.packages('rtracklayer')
# IRkernel::installspec(user = FALSE)

# RUN bash curl_data.sh
# RUN conda install nbconvert

ENTRYPOINT [ "bash", "./run_knitr.sh" ]