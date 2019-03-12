Bootstrap:docker  
From:continuumio/anaconda

%labels
    MAINTAINER Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for the annotation pipeline
    VERSION 0.1

%environment
    PATH=/opt/conda/envs/genome-annotation-1.0/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a

# GenomeThreader
echo "Installing GenomeThreader"
cd /opt
mkdir -p /opt/gth
cd /opt/gth
wget http://genomethreader.org/distributions/gth-1.7.1-Linux_x86_64-64bit.tar.gz
tar -xvf gth-1.7.1-Linux_x86_64-64bit.tar.gz
mv gth-1.7.1-Linux_x86_64-64bit 1.7.1
rm gth-1.7.1-Linux_x86_64-64bit.tar.gz
echo 'export BSSMDIR=/opt/home/gth/1.7.1/bin/bssm' >> /environment
echo 'export GTHDATADIR=/opt/home/gth/1.7.1/bin/gthdata' >> /environment
echo 'export PATH=$PATH:/opt/gth/1.7.1/bin/' >> /environment

# To be able to mount on RZcluster
mkdir -p /ifs
