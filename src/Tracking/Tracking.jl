module Tracking
export 

    addfloemasks!, 
    add_passtimes!, 
    addψs!, 
    candidate_filter_settings, 
    candidate_matching_settings,
    distance_threshold
    LogLogQuadraticTimeDistanceFunction,
    long_tracker,
    LopezAcostaTimeDistanceFunction
    normalized_maximum_crosscorrelation

using Dates: seconds, Minute, Hour, Day

include("bwtraceboundary.jl")
include("crosscorr.jl")
include("distance_thresholds.jl")
include("psi-s.jl")
include("resample-boundary.jl")
include("register.jl")
include("rotation.jl")

##### Default settings ######
# dmw: Should these be exported variables? Or should there be a function that
# returns the defaults? I'm placing them here temporarily, with the intent
# of including these in a file such as segmentation-lopez-acosta-2019.jl.

# TODO: Set these to match the Lopez-Acosta 2019 paper
# TODO: Generalize area filtering. We could have a size-dependent floe filter function. Here it's two-step, but it doesn't have to be.
# TODO: The size-dependent geometric thresholds can also be a function, e.g. geometric_filter_function(area, ratios...)
# TODO: Replace references to condition_thresholds in function documentation

candidate_filter_settings = (
    time_space_threshold_function = LopezAcostaTimeDistanceFunction(),
    small_floe_settings = (
            minimumarea=400,
            arearatio=0.18,
            majaxisratio=0.1,
            minaxisratio=0.15,
            convexarearatio=0.2,
        ),
    large_floe_settings = (
            minimumarea=1200,
            arearatio=0.28,
            majaxisratio=0.10,
            minaxisratio=0.12,
            convexarearatio=0.14,
        ),
    resolution = 250 # spatial resolution per pixel for distance computation
)

# TODO: replace all references to mc_thresholds in function calls
candidate_matching_settings = (
    goodness=(small_floe_area=0.18, # TODO: Check: how does this compare with the areas in the initial filter? Do we need it at all?
              large_floe_area=0.236,
              corr=0.68), # TODO: this correlation is too low. Should be above 0.9. Fix in next pull request.
    comp=(mxrot=10, sz=16), # TODO: Rename these variables for clarity -- we don't need to ration letters
)

export 
    candidate_filter_settings,
    candidate_matching_settings
end
