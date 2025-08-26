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
