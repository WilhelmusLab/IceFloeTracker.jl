module IceFloeTracker
using Clustering
using DataFrames
using Dates
using DelimitedFiles: readdlm, writedlm
using DSP
using ImageBinarization
using ImageContrastAdjustment
using ImageSegmentation
using Images
using Interpolations
using OffsetArrays: centered
using Peaks
using Pkg
using PyCall
using Random
using Serialization: deserialize, serialize
using StaticArrays
using StatsBase
using TiledIteration
using TOML

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
    cropfloe,
    loadimg,
    matchcorr,
    centered,
    imrotate,
    IFTVERSION

# For IFTPipeline
using HDF5
export HDF5, PyCall
export DataFrames, DataFrame, nrow, Not, select!
export Dates, Time, Date, DateTime, @dateformat_str
export addlatlon!, getlatlon, convertcentroid!, converttounits!, dropcols!

# For the tracker
export addfloemasks!, add_passtimes!, addψs!, long_tracker

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
include("tilingutils.jl")
include("histogram_equalization.jl")
include("reconstruction.jl")
include("watershed.jl")
include("brighten.jl")
include("morph_fill.jl")
include("imcomplement.jl")
include("imadjust.jl")
include("ice_masks.jl")
include("regularize-final.jl")

function get_version_from_toml(pth=dirname(dirname(pathof(IceFloeTracker))))::VersionNumber
    toml = TOML.parsefile(joinpath(pth, "Project.toml"))
    return VersionNumber(toml["version"])
end

const IFTVERSION = get_version_from_toml()

const sk_measure = PyNULL()
const sk_morphology = PyNULL()
const sk_exposure = PyNULL()
const getlatlon = PyNULL()

function __init__()
    skimage = "scikit-image=0.25.1"
    copy!(sk_measure, pyimport_conda("skimage.measure", skimage))
    copy!(sk_exposure, pyimport_conda("skimage.exposure", skimage))
    copy!(sk_morphology, pyimport_conda("skimage.morphology", skimage))
    pyimport_conda("pyproj", "pyproj=3.7.0")
    pyimport_conda("rasterio", "rasterio=1.4.3")
    @pyinclude(joinpath(@__DIR__, "latlon.py"))
    copy!(getlatlon, py"getlatlon")
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
include("tracker/long_tracker.jl")

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

include("preprocess_tiling.jl")

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
