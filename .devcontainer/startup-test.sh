#!/usr/bin/env bash

set -xeuo pipefail

cd "${CODESPACE_VSCODE_FOLDER}"  # Ensure execution in the root directory of the repository
julia -e 'using Pkg; Pkg.add("TestEnv")'  # Install TestEnv globally
julia --project=. -e 'using TestEnv; TestEnv.activate(); include("test/runtests-smoke.jl")'