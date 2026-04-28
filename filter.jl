using DataFrames
using IceFloeTracker
using Dates
using Tables

abstract type AbstractCandidateFilter <: Function end

function (f::AbstractCandidateFilter)(floe, candidates)
    transformed = DataFrames.transform(candidates, f(floe))
    filtered = DataFrames.subset(transformed, [f.threshold_column] => r -> r .> 0)
    return select(filtered, Not(f.threshold_column))
end

@kwdef struct TimeThresholdFilter <: AbstractCandidateFilter
    inputs::Vector{Symbol} = [:Δt]
    outputs::Vector{Symbol} = [:time_threshold_test]
    threshold_function::Function = ByRow((Δt) -> Δt < Hour(6))
end

function DataFrames.transform(df::DataFrame, f::AbstractCandidateFilter)
    return transform(df, f.inputs => f.threshold_function => f.outputs)
end

function DataFrames.subset(df::DataFrame, f::AbstractCandidateFilter)
    return subset(df, f.inputs => f.threshold_function)
end

candidates = DataFrame(;
    head_uuid=["f1", "f1", "f2", "f2"],
    uuid=["g1", "g2", "g1", "g2"],
    passtime=[
        DateTime(2020, 1, 1, 12, 0),
        DateTime(2020, 1, 1, 12, 5),
        DateTime(2020, 1, 1, 12, 10),
        DateTime(2020, 1, 1, 12, 15),
    ],
    Δt=[Hour(0), Hour(5), Hour(10), Hour(15)],
    row_centroid=[0.0, 0.0, 1.0, 1.0],
    col_centroid=[0.0, 0.0, 1.0, 1.0],
)

@show DataFrames.transform(candidates, TimeThresholdFilter())
@show DataFrames.subset(candidates, TimeThresholdFilter())

abstract type AbstractCandidateFilter <: Function end

function (f::AbstractCandidateFilter)(floe, candidates)
    transformed = DataFrames.transform(candidates, f(floe))
    filtered = DataFrames.subset(transformed, [f.threshold_column] => r -> r .> 0)
    return select(filtered, Not(f.threshold_column))
end

# floe = Tables.Row((;
#     :passtime => DateTime(2020, 1, 1, 12, 0), :row_centroid => 0.0, :col_centroid => 0.0
# ))

# threshold_function = DistanceThreshold()
# result = threshold_function(candidates[1, :], floe)
# @show result
# @show transform(
#     candidates,
#     [:passtime, :row_centroid, :col_centroid] =>
#         ByRow(threshold_function(floe)) =>
#             [:Δx, :Δt, :scaled_distance, :time_distance_test],
# )

# @show transform(
#     candidates, [:passtime, :row_centroid, :col_centroid] => ByRow(threshold_function(floe))
# )

# @show transform(candidates, All() => ByRow(threshold_function(floe)))

# # df2 = transform(candidates, df -> threshold_function(floe).(df))
# # @show df2
