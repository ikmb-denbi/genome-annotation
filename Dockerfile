FROM nfcore/base
MAINTAINER MTorres <m.torres@ikmb.uni-kiel.de>
LABEL authors="m.torres@ikmb.uni-kiel.de" \
    description="Docker image containing all requirements for NF-hints pipeline"

COPY environment.yml /
RUN conda env update -n root -f /environment.yml && conda clean -a
