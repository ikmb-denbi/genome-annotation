FROM nfcore/base
LABEL authors="Marc Hoeppner" \
      description="Docker image containing all requirements for ESGA annotation pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/genome-annotation-1.0/bin:$PATH
ENV PATH /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/bin:$PATH
ENV PASAHOME /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3
RUN mkdir -p /ifs
