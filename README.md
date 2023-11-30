# IceFloeTracker

[![Build Status](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wilhelmuslab.github.io/IceFloeTracker.jl/)

Track Ice Floes using Moderate Resolution Imaging Spectroradiometer (MODIS) data.

## Documentation 

See the package's documentation (in development) at https://wilhelmuslab.github.io/IceFloeTracker.jl/

## Prerequisites

A `julia` installation; ensure it is available on the `PATH`.

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
(@v1.9) pkg> activate IceFloeTracker.jl/
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

## Using in a Pipeline

See the related tools at `https://github.com/WilhelmusLab/ice-floe-tracker-pipeline.git` for a Julia Command-line Interface to access all the source functions and templates for using the Cylc pipeline orchestrator.

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