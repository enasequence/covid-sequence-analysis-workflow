#For AWS Batch Head node
#Run based on custom AMI with Google CLI + preconfig project ID attached
FROM amazoncorretto:17
RUN yum -y update && \
    yum -y install wget git python-pip jq tar.x86_64 gzip which&& \
    yum clean all
RUN pip install --upgrade awscli
# RUN curl -O https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-422.0.0-linux-x86_64.tar.gz && \
#         tar -xf google-cloud-cli-422.0.0-linux-x86_64.tar.gz && \
#         mv google-cloud-sdk /opt/google-cloud-sdk && \
#         ln -s /opt/google-cloud-sdk/bin/* /bin/ && \
#         rm google-cloud-cli-422.0.0-linux-x86_64.tar.gz
RUN curl https://sdk.cloud.google.com | bash
RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/
ARG DIR="/scratch/covid-sequence-analysis-workflow"
RUN mkdir -p $DIR
WORKDIR "$DIR"
RUN git clone -b "master" https://github.com/enasequence/covid-sequence-analysis-workflow $DIR
CMD ["${DIR}/run.nextflow.sh"]