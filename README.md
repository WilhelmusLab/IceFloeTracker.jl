# IceFloeTracker

[![Build Status](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl)

Track Ice Floes using Moderate Resolution Imaging Spectroradiometer (MODIS) data

## Notebooks

To use notebooks with `IceFloeTracker.jl` you must activate the notebooks project and start Pluto

To activate the notebooks project, start a Julia session from the root of this project and run the following commands

```
] activate "./notebooks"
```

To start pluto run the following from the same Julia session:

```
import Pluto; Pluto.run()
```

Each notebook must activate the project as well. Add the following code to your notebook

```julia
import Pkg
Pkg.activate(".")

import IceFloeTracker
```

You now have access to `IceFloeTracker` from inside your Pluto notebook!

## Fetch Data

The [`fetchdata.sh`](/scripts/fetchdata.sh) script requires the utilities [`gdal`](https://gdal.org/) and [`proj`](https://proj.org/). This repository includes a brewfile for ease of installation on MacOS. To install `gdal` and `proj` via homebrew, first [install homebrew](https://brew.sh/), then run `brew bundle install`

## Commandline

To call each step of Ice Floe Tracker pipeline from the command line you can run:

```
./scripts/ice-floe-tracker.jl
```

Each step of the pipeline is implemented as a command to the script. Most commands take a metadata file, input directory, and output directory.

For example:

```
./scripts/ice-floe-tracker.jl landmask <METADATA FILE> <INPUT DIR> <OUTPUT DIR>
```

## Snakemake

Snakemake is used to encode the entire pipeline from start to finish. Snakemake relies on the command line scripts to automate the pipeline. The snakemake file should be suitable for runs on HPC systems. To run snakemake locally, run the following from a terminal in the root of this project:

```
snakemake -c<NUM CORES>
```
