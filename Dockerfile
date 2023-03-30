#For AWS Batch Head node
FROM amazoncorretto:17
RUN yum -y update && \
    yum -y install wget && \
    yum install -y tar.x86_64 gzip && \
    yum clean all
RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-422.0.0-linux-x86_64.tar.gz 
RUN tar -xf ./google-cloud-cli-422.0.0-linux-x86_64.tar.gz
RUN mv google-cloud-sdk /opt/google-cloud-sdk && \
    ln -s /opt/google-cloud-sdk/bin/* /bin/ && \
    rm google-cloud-cli-422.0.0-linux-x86_64.tar.gz
# RUN ./google-cloud-sdk/install.sh -q --path-update=true --command-completion=true 
COPY secrets/prj-int-dev-covid19-nf-gls-aws-worker.json ./prj-int-dev-covid19-nf-gls-aws-worker.json
RUN gcloud auth activate-service-account --key-file=./prj-int-dev-covid19-nf-gls-aws-worker.json && \
    gcloud config set project "prj-int-dev-covid19-nf-gls"
RUN yum install -y git python-pip jq
RUN pip install --upgrade awscli
RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/
VOLUME ["/scratch"]
ARG DIR=/scratch/covid-sequence-analysis-workflow
RUN git clone -b "aws-batch" https://github.com/enasequence/covid-sequence-analysis-workflow $DIR
# COPY run.nextflow.sh /usr/local/bin/run.nextflow.sh
CMD ["${DIR}/run.nextflow.sh"]
