From:nfcore/base
Bootstrap:docker

%labels
    MAINTAINER MTorres <m.torres@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for NF-hints pipeline
    VERSION 0.1.0

%files
    environment.yml /

%post
    /opt/conda/bin/conda env update -n root -f /environment.yml
    /opt/conda/bin/conda clean -a
