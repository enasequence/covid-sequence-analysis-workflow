FROM python
RUN apt-get update
RUN apt-get install -y tabix
RUN pip install pysam
COPY ./vcf_to_consensus.py /vcf_to_consensus.py