# IceFloeTracker

[![Build Status](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/WilhelmusLab/IceFloeTracker.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/WilhelmusLab/IceFloeTracker.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://wilhelmuslab.github.io/IceFloeTracker.jl/)

Track Ice Floes using Moderate Resolution Imaging Spectroradiometer (MODIS) data.

## Documentation 

See the package's documentation at https://wilhelmuslab.github.io/IceFloeTracker.jl/

## Developer Quick Start

The easiest way to get started developing the IceFloeTracker.jl is to use a [devcontainer](https://containers.dev/). 
Clone the repository in VSCode and then run the command "Reopen in Container".
This will create a virtual machine to run the code, 
ensure all the packages are installed and precompiled,
and run a subset of the package tests.

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

There are Jupyter notebooks illustrating the main image processing and tracking functions, in the `/notebooks` folder. 

## Interface for Pipeline Workflows

See related tools in the [IFTPipeline repository](https://github.com/WilhelmusLab/ice-floe-tracker-pipeline#ice-floe-tracker-pipeline), including a Julia Command-line Interface and templates that leverage the [Cylc](https://cylc.github.io) pipeline orchestrator.

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

### Versioning the registered package

1. Start a new branch and update the major or minor version numbers in the corresponding field at the top of `Project.toml`
2. After merging the changes, add `@JuliaRegistrator register` in a comment in the commit you wish to use for the release (not a comment on a PR, but the actual commit)
    - this bot will open a PR on the [Julia registry](https://github.com/JuliaRegistries/General/tree/master/I/IceFloeTracker) for the package
3. Wait for feedback from the bot to make sure the new version is accepted and merged to the Julia registry
4. Create a new version tag
    - The bot will generate git commands that you can run in a terminal to add a tag
5. Create the new release on the [repo console](https://github.com/WilhelmusLab/IceFloeTracker.jl/releases)
    - Click on `Draft a new release`
    - Choose the new tag you created
    - Click on `Generate release notes` or add a custom description
6. Head over the `ice-floe-tracker-pipeline` [repo](https://github.com/WilhelmusLab/ice-floe-tracker-pipeline/blob/main/Project.toml) and update the `[compat]` section of the `Project.toml` where `IceFloeTracker` is a dependency.

**Note:** After the PR for the release is merged, a trigger workflow will force a rebuild of Docker container used in the [pipeline](https://github.com/WilhelmusLab/ice-floe-tracker-pipeline).