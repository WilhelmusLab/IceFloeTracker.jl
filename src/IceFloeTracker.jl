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
using Folds

export readdlm,
    padnhood,
    bridge,
    branch,
    @persist,
    load,
    cloudmask,
    create_cloudmask,
    deserialize,
    serialize,
    check_landmask_path,
    create_landmask,
    RGB,
    Gray,
    float64,
    imsharpen,
    label_components,
    regionprops_table,
    loadimg,
    matchcorr

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
include("tracker/tracker-funcs.jl")
include("tracker/matchcorr.jl")
include("tracker/tracker.jl")

"""
    Pipeline

This module contains the wrapper functions called by CLI.
"""
module Pipeline
    using IceFloeTracker
    using IceFloeTracker: Folds, DataFrame, RGB, Gray, load, float64, imsharpen
    using TOML: parsefile
    include("pipeline/landmask.jl")
    include("pipeline/preprocess.jl")
    include("pipeline/feature-extraction.jl")
    include("pipeline/tracker.jl")
    export sharpen,
        sharpen_gray,
        preprocess,
        cloudmask,
        extractfeatures,
        get_ice_labels,
        load_imgs,
        load_truecolor_imgs,
        load_reflectance_imgs,
        load_cloudmask,
        disc_ice_water,
        landmask,
        track
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
