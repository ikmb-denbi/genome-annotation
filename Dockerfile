FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/genome-annotation-1.0/bin:/opt/bin:$PATH
ENV PATH /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/bin:$PATH
ENV PASAHOME /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3
ENV EVM_HOME /opt/conda/envs/genome-annotation-1.0/opt/evidencemodeler-1.1.1
RUN apt-get -y install liburi-encode-perl make gcc
RUN cpan -i URI::Encode
RUN mkdir -p /ifs
RUN cp /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/pasa_conf/pasa.CONFIG.template /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/pasa_conf/conf.txt
RUN mkdir -p /opt/bin && cd /opt/bin && wget ftp://saf.bio.caltech.edu/pub/software/molbio/fastasplitn.c && gcc -o fastasplitn fastasplitn.c && chmod +x fastasplitn
