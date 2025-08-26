FROM mcr.microsoft.com/devcontainers/base:bookworm

RUN apt-get update
RUN apt-get install -y \
    build-essential \
    software-properties-common \
    gdal-bin \
    libgdal-dev \
    libgtk-3-dev \
    hdf5-tools

USER vscode
RUN curl -fsSL https://install.julialang.org | sh -s -- --yes
