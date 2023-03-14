module IceFloeTracker
using Images
using DelimitedFiles: readdlm, writedlm
using Dates
using ImageContrastAdjustment
using ImageSegmentation
using Peaks
using Random
using StatsBase
using Interpolations
using DataFrames
using PyCall
using Clustering
using DSP
using RegisterMismatch
using RegisterQD
using StaticArrays
using OffsetArrays: centered
using Serialization: serialize, deserialize

export readdlm,
    padnhood, bridge, branch, landmask, @persist, load, cloudmask, create_cloudmask

include("utils.jl")
include("persist.jl")
include("landmask.jl")
include("cloudmask.jl")
include("normalization.jl")
include("ice-water-discrimination.jl")
include("anisotropic_image_diffusion.jl")
include("bwtraceboundary.jl")
include("resample-boundary.jl")
include("psi-s.jl")
include("crosscorr.jl")
include("register-mismatch.jl")
include("bwareamaxfilt.jl")
include("hbreak.jl")
include("bridge.jl")
include("prune.jl")
include("branch.jl")
include("special_strels.jl")

const sk_measure = PyNULL()

function __init__()
    return copy!(sk_measure, pyimport_conda("skimage.measure", "scikit-image"))
end

include("regionprops.jl")
include("segmentation_a_direct.jl")
include("segmentation_b.jl")
include("segmentation_watershed.jl")
include("bwperim.jl")
include("find_ice_labels.jl")
include("segmentation_f.jl")

include("pipeline/preprocess.jl")
include("pipeline/feature-extraction.jl")

function fetchdata(; output::AbstractString)
    mkpath("$output")
    touch("$output/metadata.json")

    mkpath("$output/landmask")
    touch("$output/landmask/landmask.tiff")

    mkpath("$output/truecolor")
    touch("$output/truecolor/a.tiff")
    touch("$output/truecolor/b.tiff")
    touch("$output/truecolor/c.tiff")

    mkpath("$output/reflectance")
    touch("$output/reflectance/a.tiff")
    touch("$output/reflectance/b.tiff")
    touch("$output/reflectance/c.tiff")
    return nothing
end

"""
    MorphSE

Module for morphological operations with structuring element functionality adapted from ImageMorphology v0.4.3.

This module is temporary until v0.5 of ImageMorphology is released.

Main functionality is `dilate(img, se)` for landmask computations.

# Example

```jldoctest; setup = :(using IceFloeTracker)
julia> a = zeros(Int, 11, 11); a[6, 6] = 1;

julia> a
11×11 Matrix{Int64}:
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  1  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0

julia> se = trues(5,5);

julia> IceFloeTracker.MorphSE.dilate(a, se)
11×11 Matrix{Int64}:
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  1  1  1  1  1  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0
```
"""
module MorphSE
    using ImageCore
    using ColorTypes
    using LoopVectorization
    using OffsetArrays
    using TiledIteration: EdgeIterator
    using DataStructures
    include("morphSE/StructuringElements.jl")
    using .StructuringElements
    include("morphSE/extreme_filter.jl")
    include("morphSE/utils.jl")
    include("morphSE/dilate.jl")
    include("morphSE/erode.jl")
    include("morphSE/opening.jl")
    include("morphSE/closing.jl")
    include("morphSE/bothat.jl")
    include("morphSE/mreconstruct.jl")
    include("morphSE/fill_holes.jl")
end
end
