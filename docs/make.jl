push!(LOAD_PATH,"../src/")

# using Pkg
# pkg"activate .."

using Documenter
using IceFloeTracker

makedocs(
    sitename = "IceFloeTracker.jl",
    format = Documenter.HTML(),
    modules = [IceFloeTracker]
)

deploydocs(
    repo = "github.com/WilhelmusLab/IceFloeTracker.jl.git",
    versions = nothing,
)