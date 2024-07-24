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

const sk_measure = PyNULL()
const getlatlon = PyNULL()

function get_version_from_toml(pth=dirname(dirname(pathof(IceFloeTracker))))::VersionNumber
    toml = TOML.parsefile(joinpath(pth, "Project.toml"))
    return VersionNumber(toml["version"])
end

const IFTVERSION = get_version_from_toml()

function __init__()
    pyimport_conda("numpy", "numpy=1.23")
    pyimport_conda("pyproj", "pyproj=3.6.0")
    pyimport_conda("rasterio", "rasterio=1.3.7")
    pyimport_conda("jinja2", "jinja2=3.1.2")
    pyimport_conda("pandas", "pandas=2")
    @pyinclude(joinpath(@__DIR__, "latlon.py"))
    copy!(sk_measure, pyimport_conda("skimage.measure", "scikit-image=0.20.0"))
    copy!(getlatlon, py"getlatlon")
    return nothing
end

include("regionprops.jl")
include("utils.jl")
include("persist.jl")
include("bwareamaxfilt.jl")
end
