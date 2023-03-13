FROM amazoncorretto:17
RUN curl -s https://get.nextflow.io | bash \
 && mv nextflow /usr/local/bin/
RUN yum install -y git python-pip curl jq
RUN pip install --upgrade awscli
COPY run.aws.nextflow.sh /usr/local/bin/run.aws.nextflow.sh
VOLUME ["/scratch"]
CMD ["/usr/local/bin/run.aws.nextflow.sh"]
COPY nextflow-lib/nextflow.config /scratch/nextflow-lib/nextflow.config