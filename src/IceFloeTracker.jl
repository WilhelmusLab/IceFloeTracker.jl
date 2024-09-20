module IceFloeTracker
using Images
using DelimitedFiles: readdlm, writedlm
using Dates
using ImageContrastAdjustment
using ImageSegmentation
using Peaks
using Pkg
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
const getlatlon = PyNULL()

function get_version_from_toml(pth=dirname(dirname(pathof(IceFloeTracker))))::VersionNumber
    toml = TOML.parsefile(joinpath(pth, "Project.toml"))
    return VersionNumber(toml["version"])
end

const IFTVERSION = get_version_from_toml()

function parse_requirements(file_path)
    requirements = Dict{String, String}()
    open(file_path, "r") do f
        for line in eachline(f)
            pkg, version = split(line, "==")
            requirements[pkg] = version
        end
    end
    return requirements
end

function __init__()
    deps = parse_requirements(joinpath(dirname(@__DIR__), "requirements.txt"))

    for (pkg, version) in deps
        if pkg == "scikit-image"
            sk_measure_module = pyimport_conda("skimage.measure", "$(pkg)=$(version)")
            copy!(sk_measure, sk_measure_module)
        else
            pyimport_conda(pkg, "$(pkg)=$(version)")
        end
    end

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
