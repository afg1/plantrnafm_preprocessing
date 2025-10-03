FROM ubuntu:22.04

# Basic setup
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    build-essential \
    wget \
    python3-pip \
    python3-venv \
    curl

# Install BLAST+
WORKDIR /opt
RUN wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-2.17.0+-x64-linux.tar.gz
RUN tar -zxf ncbi-blast-2.17.0+-x64-linux.tar.gz
RUN rm ncbi-blast-2.17.0+-x64-linux.tar.gz

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
ENV PATH="/opt/venv/bin:/opt/ncbi-blast-2.17.0+/bin:${PATH}"

