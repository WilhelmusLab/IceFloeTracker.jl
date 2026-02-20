#!/usr/bin/env bash

set -xeuo pipefail
sudo apt update
sudo apt install -y libgl1 libgdal-dev libglib2.0-0 -y
pipx install snakemake