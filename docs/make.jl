using Documenter
using IceFloeTracker

makedocs(;
    sitename="IceFloeTracker.jl",
    format=Documenter.HTML(; size_threshold=nothing),
    modules=[IceFloeTracker],
    doctest=false,
    warnonly=true,
)

deploydocs(;
    repo="github.com/WilhelmusLab/IceFloeTracker.jl.git",
    push_preview=true,
    versions=nothing,
)
