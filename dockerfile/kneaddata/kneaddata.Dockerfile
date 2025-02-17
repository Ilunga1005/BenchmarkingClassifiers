FROM ubuntu:24.04
LABEL maintainer="zhouan@genomics.cn"
LABEL software="Kraken2"
LABEL kneaddata_version="0.12.1"
LABEL fastqc_version="0.12.1"
RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y python3 python3-dev python3-pip apt-transport-https openjdk-8-jre wget zip && \
    rm -rf /var/lib/apt/lists/* && apt-get autoclean
RUN pip3 install boto3 cloudpickle awscli &&\
    pip3 install anadama2 && \
    pip3 install kneaddata==${kneaddata_version} --no-binary :all:
# install kneaddata and dependencies

# install fastqc
RUN wget https://github.com/s-andrews/FastQC/archive/refs/tags/v${fastqc_version}.zip && \
    unzip v${fastqc_version}.zip && \
    mv v${fastqc_version} FastQC &&\
    chmod 755 FastQC/fastqc && \
    mv FastQC /usr/local/bin/ && \
    ln -s /usr/local/bin/FastQC/fastqc /usr/local/bin/fastqc && \
    rm v${fastqc_version}.zip 

WORKDIR /data
