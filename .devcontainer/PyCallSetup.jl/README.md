# Python Setup for IFTPipeline.jl

This project is designed to initialize the Python Environment for IFTPipeline.jl.
It should be used once before IFTPipeline.jl is instantiated.
It will:

- Build PyCall using the `Conda.jl` package.

Since the `base` environment is shared between all PyCall instances by default, it only needs to be run if the Conda environment is:
- missing – as in a new installation, or 
- broken – new incompatible dependencies have been added.

The goal is to remove the need for this subpackage at some point in future. 

## User Guide

To use, run:
```bash
julia --project="PyCallSetup.jl" PyCallSetup.jl/setup.jl
```

> If this fails, try removing the existing Conda installation for Julia:
> ```bash
> rm -r ~/.julia/conda/
> ```
> ... and then retry.

If you're working on a system where the home directory might not be reliable, 
e.g. you're building a Docker container but you want to use it with Apptainer 
which mounts things in different ways, you may want to change 
the default location for the Conda installation.

Export a different path in `CONDA_JL_HOME` before running setup:
```bash
export CONDA_JL_HOME="/opt/conda-env/"
julia --project="PythonSetupForIFTPipeline.jl" PythonSetupForIFTPipeline.jl/setup.jl
```