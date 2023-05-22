#For AWS Batch Head node
#Run based on custom AMI with Google CLI + preconfig project ID attached
FROM amazoncorretto:17
RUN yum -y update && \
    yum -y install wget git python-pip jq tar.x86_64 gzip which&& \
    yum clean all
ENV CLOUDSDK_PYTHON=python3.8
RUN pip install --upgrade awscli
RUN curl -sL https://sdk.cloud.google.com/install/linux | bash -s -- --disable-prompts
# RUN curl https://sdk.cloud.google.com | bash -s -- --install-dir=/bin/ --disable-prompts
# ENV PATH $PATH:/bin/google-cloud-sdk/bin
RUN gcloud version && bq version
RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/
ARG DIR="/scratch/covid-sequence-analysis-workflow"
RUN mkdir -p $DIR
WORKDIR "$DIR"
RUN git clone -b "master" https://github.com/enasequence/covid-sequence-analysis-workflow $DIR
CMD ["${DIR}/run.nextflow.sh"]