module Tracking

export addfloemasks!,
    add_passtimes!,
    addψs!,
    adduuid!,
    align_centroids,
    buildψs,
    bwtraceboundary,
    candidate_filter_settings,
    candidate_matching_settings,
    compute_centroid,
    cropfloe,
    crosscorr,
    distance_threshold,
    dropcols!,
    get_rotation_measurements,
    get_trajectory_heads,
    grad,
    imrotate_bin_counterclockwise_radians,
    matchcorr,
    mismatch,
    normalized_cross_correlation,
    LinearTimeDistanceFunction,
    LogLogQuadraticTimeDistanceFunction,
    long_tracker,
    LopezAcostaTimeDistanceFunction,
    PiecewiseLinearThresholdFunction,
    register,
    relative_error_test!,
    resample_boundary,
    shape_difference_rotation,
    shape_difference_test!,
    time_distance_test!,
    _add_suffix,
    norm

include("distance_thresholds.jl")
include("bwtraceboundary.jl")
include("crosscorr.jl")
include("extended_regionprops.jl")
include("long_tracker.jl")
include("matchcorr.jl")
include("psi-s.jl")
include("register.jl")
include("resample-boundary.jl")
include("rotation.jl")
include("tracker-funcs.jl")
include("tracker.jl")

##### Default settings ######
# TODO: Set these to match the Lopez-Acosta 2019 paper
# TODO: Replace two-level geometry threshold function with a functor or function
# TODO: Replace references to condition_thresholds in function documentation

candidate_filter_settings = (
    time_space_threshold_function=LopezAcostaTimeDistanceFunction(),
    small_floe_settings=(
        minimumarea=400,
        arearatio=0.18,
        majaxisratio=0.1,
        minaxisratio=0.15,
        convexarearatio=0.2,
    ),
    large_floe_settings=(
        minimumarea=1200,
        arearatio=0.28,
        majaxisratio=0.10,
        minaxisratio=0.12,
        convexarearatio=0.14,
    ),
    resolution=250, # spatial resolution per pixel for distance computation
)

# TODO: replace all references to mc_thresholds in function calls
candidate_matching_settings = (
    goodness=(
        small_floe_area=0.18, # TODO: Check: how does this compare with the areas in the initial filter? Do we need it at all?
        large_floe_area=0.236,
        corr=0.68,
    ), # TODO: this correlation is too low. Should be above 0.9. Fix in next pull request.
    comp=(mxrot=10, sz=16), # TODO: Rename these variables for clarity -- we don't need to ration letters
)

export candidate_filter_settings, candidate_matching_settings

end
