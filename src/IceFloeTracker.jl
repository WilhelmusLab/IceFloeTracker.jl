module IceFloeTracker

include("Segmentation/Segmentation.jl")
using .Segmentation

include("Filtering/Filtering.jl")
using .Filtering

include("Morphology/Morphology.jl")
using .Morphology

include("Tracking/Tracking.jl")
using .Tracking

include("Preprocessing/Preprocessing.jl")
using .Preprocessing

include("Utils/Utils.jl")
using .Utils

include("Data/Data.jl")
using .Data

using Clustering
using DataFrames
using Dates
using DelimitedFiles: readdlm, writedlm
using DSP
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
    LopezAcostaCloudMask,
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
    IFTVERSION,
    get_rotation_measurements,
    IceFloeSegmentationAlgorithm,
    LopezAcosta2019,
    ValidationDataCase,
    ValidationDataLoader,
    ValidationDataSet,
    Watkins2025GitHub,
    segmentation_comparison,
    segmentation_summary,
    LopezAcosta2019Tiling,
    callable_store,
    binarize_segments,
    masker,
    IceDetectionAlgorithm,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionThresholdMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionLopezAcosta2019,
    tiled_adaptive_binarization

# For IFTPipeline
using HDF5
export HDF5, PyCall
export DataFrames, DataFrame, nrow, Not, select!
export Dates, Time, Date, DateTime, @dateformat_str
export addlatlon!, convertcentroid!, converttounits!, dropcols!, latlon

# For the tracker
export addfloemasks!, 
    add_passtimes!, 
    addψs!, 
    candidate_filter_settings, 
    candidate_matching_settings,
    distance_threshold,
    LogLogQuadraticTimeDistanceFunction,
    long_tracker,
    LopezAcostaTimeDistanceFunction

include("utils.jl")
include("landmask.jl")
include("cloudmask.jl")
include("normalization.jl")
include("ice-water-discrimination.jl")
include("nonlinear_diffusion.jl")
include("tilingutils.jl")
include("histogram_equalization.jl")
include("reconstruction.jl")
include("watershed.jl")
include("brighten.jl")
include("imcomplement.jl")
include("imadjust.jl")
include("ice_masks.jl")
include("regularize-final.jl")
include("latlon.jl")
include("segmentation-lopez-acosta-2019.jl")
include("segmented-image-utilities.jl")

function get_version_from_toml(pth=dirname(dirname(pathof(IceFloeTracker))))::VersionNumber
    toml = TOML.parsefile(joinpath(pth, "Project.toml"))
    return VersionNumber(toml["version"])
end

const IFTVERSION = get_version_from_toml()

const sk_measure = PyNULL()
const sk_morphology = PyNULL()
const sk_exposure = PyNULL()

function __init__()
    skimage = "scikit-image=0.25.1"
    copy!(sk_measure, pyimport_conda("skimage.measure", skimage))
    copy!(sk_exposure, pyimport_conda("skimage.exposure", skimage))
    copy!(sk_morphology, pyimport_conda("skimage.morphology", skimage))
    return nothing
end

include("regionprops.jl")
include("segmentation_a_direct.jl")
include("segmentation_b.jl")
include("segmentation_watershed.jl")
include("find_ice_labels.jl")
include("segmentation_f.jl")
include("segmentation-lopez-acosta-2019-tiling.jl")
include("mask.jl")
end
