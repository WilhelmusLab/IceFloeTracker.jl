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
    matchcorr,
    centered,
    imrotate

# For IFTPipeline
using HDF5
export HDF5, PyCall
export DataFrames, DataFrame, nrow, Not, select!
export Dates, Time, Date, DateTime, @dateformat_str

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
    copy!(sk_measure, pyimport_conda("skimage.measure", "scikit-image=0.20.0"))
    # pyimport_conda("numpy", "numpy=1.25.0")
    pyimport_conda("pyproj", "pyproj=3.6.0")
    pyimport_conda("rasterio", "rasterio=1.3.7")
    pyimport_conda("requests", "requests=2.31.0")
    pyimport_conda("skyfield", "skyfield=1.45.0")
    return nothing
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

module Register
    include("Register/CenterIndexedArrays.jl-0.2.0/CenterIndexedArrays.jl")
    include("Register/RegisterCore.jl-0.2.4/src/RegisterCore.jl")
    include("Register/RegisterMismatchCommon.jl-master/RegisterMismatchCommon.jl")
    include("Register/RegisterUtilities.jl-master/RegisterUtilities.jl")
    include("Register/RFFT.jl-master/RFFT.jl")
    include("Register/RegisterDeformation.jl-0.4.4/RegisterDeformation.jl")
    include("Register/QuadDIRECT.jl-master/QuadDIRECT.jl")
    include("Register/RegisterQD.jl-0.3.1/RegisterQD.jl")
    include("Register/RegisterMismatch.jl-0.4.0/RegisterMismatch.jl")
end
end
