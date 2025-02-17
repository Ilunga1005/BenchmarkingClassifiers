FROM ubuntu:24.04
# FROM continuumio/miniconda3:latest

# Add metadata
LABEL maintainer="zhouan@genomics.cn"
LABEL version="v0.2.5"

# Set environment variables to prevent prompts
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    PATH="/opt/conda/bin:$PATH"

RUN apt-get update && apt-get install -y wget bzip2 && \
    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm -rf /tmp/miniconda.sh && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Install melon in the base environment
RUN conda install -n base -c conda-forge -c bioconda melon -y && \
    conda clean -afy

# Default shell for proper conda initialization
SHELL ["/bin/bash", "-c"]

# Set the default command
CMD ["melon", "-h"]
