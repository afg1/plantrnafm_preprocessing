FROM ubuntu:22.04

# Basic setup
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    python3-pip \
    python3-venv \
    curl

# Install MMseqs2
WORKDIR /opt
RUN wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz
RUN tar -zxf mmseqs-linux-avx2.tar.gz
RUN rm mmseqs-linux-avx2.tar.gz

# Install uv to a known location
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}"

# Create venv and install packages
RUN /root/.local/bin/uv venv --python 3.10 /opt/venv
RUN . /opt/venv/bin/activate && uv pip install \
    black \
    polars \
    pyarrow \
    pyfastx \
    ruff

# Set the environment
ENV PATH="/opt/venv/bin:/opt/mmseqs/bin:${PATH}"

