module IceFloeTracker

include("skimage/skimage.jl")
using .skimage

include("Utils/Utils.jl")
using .Utils

include("Data/Data.jl")
using .Data

include("ImageUtils/ImageUtils.jl")
using .ImageUtils

include("Geospatial/Geospatial.jl")
using .Geospatial

include("Morphology/Morphology.jl")
using .Morphology

include("Filtering/Filtering.jl")
using .Filtering

include("Preprocessing/Preprocessing.jl")
using .Preprocessing

include("Segmentation/Segmentation.jl")
using .Segmentation

include("Tracking/Tracking.jl")
using .Tracking

include("segmentation-lopez-acosta-2019-tiling.jl")

using Clustering
using DataFrames
using Dates
using DSP
using Images
using Interpolations
using OffsetArrays: centered
using Peaks
using Random
using StaticArrays
using StatsBase
using TiledIteration

export padnhood,
    bridge,
    branch,
    @persist,
    load,
    cloudmask,
    create_cloudmask,
    LopezAcostaCloudMask,
    Watkins2025CloudMask,
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
    get_rotation_measurements,
    IceFloeSegmentationAlgorithm,
    ValidationDataCase,
    ValidationDataLoader,
    ValidationDataSet,
    Watkins2025GitHub,
    segmentation_comparison,
    segmentation_summary,
    callable_store,
    binarize_segments,
    masker,
    IceDetectionAlgorithm,
    IceDetectionBrightnessPeaksMODIS721,
    IceDetectionThresholdMODIS721,
    IceDetectionFirstNonZeroAlgorithm,
    IceDetectionLopezAcosta2019,
    tiled_adaptive_binarization,
    LopezAcosta2019Tiling

# For IFTPipeline
export addlatlon!, convertcentroid!, converttounits!, dropcols!, latlon

# For the tracker
export addfloemasks!,
    addlatlon!,
    add_passtimes!,
    addψs!,
    candidate_filter_settings,
    candidate_matching_settings,
    distance_threshold,
    LogLogQuadraticTimeDistanceFunction,
    long_tracker,
    LopezAcostaTimeDistanceFunction,
    register,
    resample_boundary

include("segmentation-lopez-acosta-2019.jl")

end
