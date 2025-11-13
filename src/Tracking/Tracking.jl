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
    DistanceThresholdFilter,
    get_rotation_measurements,
    grad,
    imrotate_bin_counterclockwise_radians,
    matchcorr,
    mismatch,
    norm,
    normalized_cross_correlation,
    LinearTimeDistanceFunction,
    LogLogQuadraticTimeDistanceFunction,
    long_tracker,
    LopezAcostaTimeDistanceFunction,
    PiecewiseLinearThresholdFunction,
    PsiSCorrelationThresholdFilter,
    register,
    RelativeErrorThresholdFilter,
    resample_boundary,
    shape_difference_rotation,
    ShapeDifferenceThresholdFilter,
    StepwiseLinearThresholdFunction,
    time_distance_test!,
    _add_suffix

include("distance_thresholds.jl")
include("bwtraceboundary.jl")
include("crosscorr.jl")
include("extend_regionprops.jl")
include("geometric_thresholds.jl")
include("filter_functions.jl")
include("floe_tracker.jl")
include("matchcorr.jl")
include("psi-s.jl")
include("register.jl")
include("resample-boundary.jl")
include("rotation.jl")
include("tracker-funcs.jl")

##### Default settings ######
# TODO: Replace with filter_function

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
