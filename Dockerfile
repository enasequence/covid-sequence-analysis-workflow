#For AWS Batch Head node
#Run based on custom AMI with Google CLI + preconfig project ID attached
FROM amazoncorretto:17
RUN yum -y update && \
    yum -y install wget git python-pip jq tar.x86_64 gzip&& \
    yum clean all
RUN pip install --upgrade awscli
RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/
ARG DIR=/scratch/covid-sequence-analysis-workflow
RUN mkdir -p $DIR
WORKDIR "$DIR"
RUN git clone -b "aws-batch" https://github.com/enasequence/covid-sequence-analysis-workflow $DIR
CMD ["${DIR}/run.nextflow.sh"]