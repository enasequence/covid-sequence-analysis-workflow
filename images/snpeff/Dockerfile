FROM biocontainers/snpeff:v4.1k_cv3
RUN mkdir -p /home/biodocker/bin/snpEff/data/sars.cov.2
ADD genes.gbk /home/biodocker/bin/snpEff/data/sars.cov.2
RUN echo '# Database for SARS-CoV-2 (NC_045512.2)' >>/home/biodocker/bin/snpEff/snpEff.config
RUN echo 'sars.cov.2.genome : SARS-CoV-2' >>/home/biodocker/bin/snpEff/snpEff.config
RUN echo '        sars.cov.2.chromosomes : NC_045512' >>/home/biodocker/bin/snpEff/snpEff.config
RUN snpeff build -genbank -v sars.cov.2