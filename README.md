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

## SOIT Integration

The [Satellite Overpass Identification Tool](https://zenodo.org/record/6475619#.ZBhat-zMJUe) is called to generate a list of satellite times for both Aqua and Terra in the area of interest. This program is written in Python and it's dependencies are added to the `oscar-env.yaml` that we build in the next section.

## Cylc to run the pipeline

Cylc is used to encode the entire pipeline from start to finish. Cylc relies on the command line scripts to automate the pipeline. The `flow.cylc` file should be suitable for runs on HPC systems. To run cylc locally, there are a few commands to run from a terminal in the root directory of this project:

```
cylc install -n <workflow-name> ./cylc
cylc graph <workflow-name>
cylc play <workflow-name>
cylc tui <workflow-name>
```
The text-based user interface provides a simple way to watch the status of each task called in the `flow.cylc` workflow. Use arrow keys to investigate each task as see more [here](https://cylc.github.io/cylc-doc/latest/html/7-to-8/major-changes/ui.html#cylc-tui).
![tui](tui-example.png)

Remember to remove the `images` and `output` folders from the root project directory and use a new workflow name when running the pipeline again.

## Running the workflow on Oscar
#### Python=3.9.13

1. `ssh` to Oscar
2. Move to a compute node
    * `interact -n 20 -t 24:00:00 -m 16g`
    * this will start a compute session for 1 day with 16 GB memory and 20 cores
    * see [here](https://docs.ccv.brown.edu/oscar/submitting-jobs/interact) for more options
3. Load the latest anaconda module
    * `module load anaconda/2022.05`
    * `source /gpfs/runtime/opt/anaconda/2022.05/etc/profile.d/conda.sh`
4. Build a virtual environment
    * `conda create -n icefloe-oscar`
    * `conda activate icefloe-oscar`
    * `conda install -c conda-forge mamba`
    * `git clone https://github.com/WilhelmusLab/IceFloeTracker.jl.git`
    * `cd IceFloeTracker.jl`
    * `mamba env update -n icefloe-oscar -f ./hpc/oscar-env.yaml`
5. Make sure the HolyLab registry is added as described in the [prerequisites section](#prerequisites)
6. Build the package
    * `julia -e 'ENV["PYTHON"]="~/anaconda/icefloe-oscar/bin/python3.10"'`
    * `julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.build'`
    * `julia -e 'using Pkg; Pkg.activate("scripts"); Pkg.instantiate(); Pkg.build'`
7. Register an account with [space-track.org](https://www.space-track.org/) for SOIT
8. Export SOIT username/password to environment variable
    * From your home directory`nano .bash_profile`
    * add `export HISTCONTROL=ignoreboth` to the bottom of your .bash_profile
        * this will ensure that your username/password are not stored in history
        * when exporting the following environment variables, there must a space in front of each command
    * ` export SPACEUSER=<firstname>_<lastname>@brown.edu`
    * ` export SPACEPSWD=<password>`
9. Run the workflow with Cylc
    ````
    cylc install -n <workflow-name> ./cylc
    cylc play <workflow-name>
    cylc tui <workflow-name>
    ```