# Filter functions and structs. These 
# 


abstract type FloeFilterFunction end
function (f::FloeFilterFunction)(floe, candidates)
    f(floe, candidates, Val(:raw))
    subset!(candidates, [f.threshold_column] => r -> r .> 0)
    select!(candidates, Not(f.threshold_column))
end

# And here is a more advanced struct/functor pair
@kwdef struct DistanceThresholdFilter <: FloeFilterFunction
        time_column = :Δt
        dist_column = :Δx
        threshold_function = LinearTimeDistanceFunction()
        threshold_column = :time_distance_test
end

function (f::DistanceThresholdFilter)(floe, candidates, _::Val{:raw}) # can we get the same behavior with a less opaque function call?
    candidates[!, f.time_column] = candidates[!, :passtime] .- floe.passtime
    candidates[!, f.dist_column] = euclidean_distance(floe, candidates)
    transform!(candidates, [f.dist_column, f.time_column] => 
        ByRow(f.threshold_function) => f.threshold_column)
end

"""
    euclidean_distance(floe, candidates; r=250)

Compute the distance in meters between a floe and candidate floes by computing the
straight-line distance between centroids in pixel coordinates and converting that result
using a pixel resolution `r` with units meters/pixel.
"""
function euclidean_distance(floe, candidates; r = 250)
    return sqrt.((floe.row_centroid .- candidates.row_centroid).^2 .+ 
    (floe.col_centroid .- candidates.col_centroid).^2) * r
end