FROM ubuntu:24.04
LABEL maintainer="zhouan@genomics.cn"
LABEL software="Kraken2"
LABEL version="v2.1.3"
ARG K2VER="2.1.3"

# install dependencies and cleanup apt garbage
RUN apt-get update && apt-get -y --no-install-recommends install \
 wget \
 ca-certificates \
 zlib1g-dev \
 make \
 g++ \
 rsync \
 git \
 cpanminus && \
 rm -rf /var/lib/apt/lists/* && apt-get autoclean

# perl module required for kraken2-build
RUN cpanm Getopt::Std
# DL Kraken2, unpack, and install
RUN wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v${K2VER}.tar.gz && \
 tar -xzf v${K2VER}.tar.gz && \
 rm -rf v${K2VER}.tar.gz && \
 cd kraken2-${K2VER} && \
 ./install_kraken2.sh . && \
 cd .. && \
 git clone https://github.com/jenniferlu717/Bracken.git && \
 cd Bracken && \
 bash install_bracken.sh && \
 mkdir /data

ENV PATH="$PATH:/kraken2-${K2VER}:/Bracken" \
    LC_ALL=C
WORKDIR /data