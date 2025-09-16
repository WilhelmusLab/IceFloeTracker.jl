FROM mcr.microsoft.com/devcontainers/base:bookworm

RUN apt-get update
RUN apt-get install -y \
    gdal-bin \
    libgdal-dev \
    hdf5-tools

USER vscode
RUN curl -fsSL https://install.julialang.org | sh -s -- --yes --default-channel lts
