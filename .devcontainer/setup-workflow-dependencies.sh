#!/usr/bin/env bash

set -xeuo pipefail
sudo apt-get update
sudo apt-get install -y libgl1 libgdal-dev libglib2.0-0

uv tool install --upgrade --with pyproj --with pandas snakemake
