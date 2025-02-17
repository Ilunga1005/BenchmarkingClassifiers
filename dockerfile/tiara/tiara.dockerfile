# Use an official Python image as the base image
FROM python:3.9-slim
LABEL maintainer="zhouan@genomics.cn"
LABEL version="1.0.3"
# Set environment variables to prevent prompts
ENV DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    gcc \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /app

# Install Tiara using pip
RUN pip install --no-cache-dir tiara

# Test the installation
RUN tiara-test

# Set the default command
CMD ["bash"]
