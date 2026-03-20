#!/usr/bin/env bash

set -xeuo pipefail
sudo apt update
sudo apt install -y libgl1 libgdal-dev libglib2.0-0 -y

# Set up Snakemake and tools needed in the snakemake workflow
pipx install snakemake
pipx inject snakemake pyproj pandas