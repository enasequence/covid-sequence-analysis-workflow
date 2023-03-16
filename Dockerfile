FROM amazoncorretto:17
RUN yum install -y curl tar wget
RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-422.0.0-linux-x86_64.tar.gz 
RUN echo $(ls)
RUN tar -xf google-cloud-cli-422.0.0-linux-x86_64.tar.gz
RUN google-cloud-sdk/install.sh -q --path-update=true --command-completion=true 
COPY secrets/prj-int-dev-covid19-nf-gls-aws-worker.json ./prj-int-dev-covid19-nf-gls-aws-worker.json
RUN gcloud auth activate-service-account --key-file=./prj-int-dev-covid19-nf-gls-aws-worker.json && \
    gcloud config set project "prj-int-dev-covid19-nf-gls" && \
    gcloud storage ls
RUN yum install -y git python-pip jq
# SHELL ["/bin/bash", "-c"]
RUN pip install --upgrade awscli
RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/
COPY run.aws.nextflow.sh /usr/local/bin/run.aws.nextflow.sh
VOLUME ["/scratch"]
CMD ["/usr/local/bin/run.aws.nextflow.sh"]
COPY nextflow-lib/nextflow.config /scratch/nextflow-lib/nextflow.config