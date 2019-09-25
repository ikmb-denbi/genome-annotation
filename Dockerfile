FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/genome-annotation-1.0/bin:$PATH
ENV PATH /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/bin:$PATH
ENV PASAHOME /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3
ENV EVM_HOME /opt/conda/envs/genome-annotation-1.0/opt/evidencemodeler-1.1.1

RUN mkdir -p /ifs
RUN cp /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/pasa_conf/pasa.CONFIG.template /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/pasa_conf/conf.txt

RUN mkdir -p /opt/blast && \
cd /opt/blast && \
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz && \
tar -xvf ncbi-blast-2.9.0+-x64-linux.tar.gz && \
mv ncbi-blast-2.9.0+ 2.9.0 && \
rm ncbi-blast-2.9.0+-x64-linux.tar.gz
