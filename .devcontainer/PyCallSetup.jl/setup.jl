#!/usr/bin/env julia

# Set up PyCall
#
# On Linux, by default, PyCall will use the system Python, 
# but this doesn't work with IceFloeTracker.jl which needs Conda.
# This script initializes PyCall to use Conda instead, 
# and the same Conda environment is shared by all PyCall instances in the container.

using Pkg

# Force PyCall to use the Conda version on Linux.
ENV["PYTHON"] = ""

# Build PyCall with the new conda environment
Pkg.add("PyCall")
Pkg.build("PyCall")

# Fix an error which occurs for some versions of Conda on Linux
# where scipy can't be imported.
# Inspired by:
# https://discourse.julialang.org/t/version-cxxabi-1-3-15-not-found-when-loading-matplotlib-through-pythoncall/131671/12?page=2
Pkg.add("Conda")
using Conda
if Sys.islinux()
    Conda.add("libstdcxx<14.0")
end