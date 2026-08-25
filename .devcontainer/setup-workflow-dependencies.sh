#!/usr/bin/env bash

set -xeuo pipefail
sudo apt-get update
sudo apt-get install -y libgl1 libgdal-dev libglib2.0-0

# Set up uv and tools needed in the snakemake workflow
if ! command -v uv >/dev/null 2>&1; then
	curl -LsSf https://astral.sh/uv/install.sh | sh
	export PATH="$HOME/.local/bin:$PATH"
fi

uv tool install --upgrade --with pyproj --with pandas snakemake
