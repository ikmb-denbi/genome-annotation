Bootstrap:docker  
From:continuumio/anaconda

%labels
    MAINTAINER Marc Hoeppner <m.hoeppner@ikmb.uni-kiel.de>
    DESCRIPTION Singularity image containing all requirements for the annotation pipeline
    VERSION 0.1

%environment
    PATH=/opt/conda/envs/genome-annotation-0.1/bin:$PATH
    export PATH

%files
    environment.yml /

%post
    /opt/conda/bin/conda env create -f /environment.yml
    /opt/conda/bin/conda clean -a

apt-get -y update
locale-gen en_US en_US.UTF-8 de_DE.UTF-8 de_DE

# InterProScan
echo "Installing InterProScan"
cd /opt
mkdir -p /opt/interproscan
cd /opt/interproscan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.33-72.0/interproscan-5.33-72.0-64-bit.tar.gz
tar -pzxvf interproscan-5.33-72.0-64-bit.tar.gz
rm interproscan-5.3*.tar.gz
mv interproscan-5.33-72.0 5.33-72.0
echo 'export PATH=$PATH:/opt/interproscan/5.33-72.0/' >> /environment
cd

# GenomeThreader
echo "Installing GenomeThreader"
cd /opt
mkdir -p /opt/gth
cd /opt/gth
wget http://genomethreader.org/distributions/gth-1.7.1-Linux_x86_64-64bit.tar.gz
tar -zxvf gth-1.7.1-Linux_x86_64-64bit.tar.gz 
cd gth-1.7.1-Linux_x86_64-64bit
mv gth-1.7.1-Linux_x86_64-64bit 1.7.1
setenv $BSSMDIR         "${HOME}/opt/gth/1.7.1/bin/bssm"
setenv $GTHDATADIR      "${HOME}/opt/gth/1.7.1/bin/gthdata"
echo 'export PATH=$PATH:/opt/gth/1.7.1/bin/' >> /environment
