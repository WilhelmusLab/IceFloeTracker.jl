#!/usr/bin/env bash

set -xeuo pipefail

julia -e 'using Pkg; Pkg.add("TestEnv")'  # Install TestEnv globally
julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/runtests-smoke.jl")'