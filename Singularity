Bootstrap:docker  
From:nfcore/base

%labels
    MAINTAINER Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for the annotation pipeline
    VERSION 0.1

%environment
    PATH=/opt/conda/envs/genome-annotation-1.0/bin:/opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/bin:$PATH
    export PATH

    PASAHOME=/opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3
    export PASAHOME

    EVM_HOME=/opt/conda/envs/genome-annotation-1.0/opt/evidencemodeler-1.1.1
    export EVM_HOME

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a

# Prereqs for Nextflow
apt-get -y install procps  liburi-encode-perl

# Create mount point for RZCluster
mkdir -p /ifs

# Create the default config file
cp /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/pasa_conf/pasa.CONFIG.template /opt/conda/envs/genome-annotation-1.0/opt/pasa-2.3.3/pasa_conf/conf.txt
