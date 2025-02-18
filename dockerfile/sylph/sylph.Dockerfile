# FROM ubuntu:24.04


# Use an official Rust image as the base image
FROM rust:latest AS builder
LABEL maintainer="zhouan@genomics.cn"
LABEL version="v0.8.0"
# Install dependencies for building Sylph
RUN apt-get update && apt-get install -y \
    gcc \
    make \
    cmake \
    && rm -rf /var/lib/apt/lists/*

# Set the working directory
WORKDIR /usr/src/sylph

# Clone the Sylph repository
RUN git clone https://github.com/bluenote-1577/sylph .

# Build Sylph using cargo
RUN cargo install --path . --root /usr/local


# Create a minimal runtime image
FROM ubuntu:24.04

# Install runtime dependencies
RUN apt-get update && apt-get install -y \
    libgcc-s1 \
    libstdc++6 \
    python3 \ 
    python3-pip \
    git \
    && rm -rf /var/lib/apt/lists/*
## install sylph-tax
RUN git clone https://github.com/bluenote-1577/sylph-tax && \
    cd sylph-tax && \
    pip install . --break-system-packages && \
    cd .. && \
    rm -rf sylph-tax
# Copy the built binary from the builder stage
COPY --from=builder /usr/local /usr/local
# Add the binary location to PATH
ENV PATH="/usr/local/bin:${PATH}"
RUN mkdir -p /usr/local/sylph_tax_db && \
    /usr/local/bin/sylph-tax download --download-to /usr/local/sylph_tax_db
# Set the default command
CMD ["sylph -h"]

