# IceFloeTracker

[![Build Status](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wilhelmuslab.github.io/IceFloeTracker.jl/)

Track Ice Floes using Moderate Resolution Imaging Spectroradiometer (MODIS) data.

## Documentation 

See the package's documentation (in development) at https://wilhelmuslab.github.io/IceFloeTracker.jl/
## Prerequisites
1. A `julia` installation; ensure it is available on the `PATH`.
2. Add the HolyLab registry by invoking the following command from the terminal
```
julia -e 'using Pkg; Pkg.Registry.add(RegistrySpec(url = "https://github.com/HolyLab/HolyLabRegistry.git"))' 
```

## Clone repo and run tests

Clone the repository.
```zsh
$ git clone https://github.com/WilhelmusLab/IceFloeTracker.jl
```

Now start a Julia session.
```zsh
$ julia
```

```
julia> ]
```
... to enter package mode.

```
(@v1.7) pkg> activate IceFloeTracker.jl/
Activating project at `~/IceFloeTracker.jl`
```

Instantiate the environment and run the tests:
```
(IceFloeTracker) pkg> instantiate
(IceFloeTracker) pkg> test
```


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

### Example

```
$ ./scripts/fetchdata.sh -o data -s 2022-05-01 81 -22 79 -12
```

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

## Development

Git hooks are used to run common developer tasks on commits (e.g. code formatting, tests, etc.). If you are running git version 2.9 or later run the following from the root of the project to enable git hooks.

```
git config core.hooksPath ./hooks
```

To help with passing git hooks, run the formatting script before staging files:

```
./scripts/format.jl
git add .
git commit -m "some informative message"
git push
```
