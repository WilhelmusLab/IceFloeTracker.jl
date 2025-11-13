import DataFrames: DataFrames, DataFrame, AbstractDataFrame, eachrow, select!, subset!
import Dates: Day

"""
    floe_tracker(props; filter_function, matching_function, minimum_floe_size, maximum_floe_size, maximum_time_step)

Track ice floes over multiple observations.

Trajectories are built as follows:
- Assume the floes detected in observation 1 are trajectories of length 1.
- For each subsequent observation at time `t``:
  - Determine the latest observation for each trajectory – these are the "current trajectory heads".
  - Select the subset of trajectory heads observed within the window `maximum_time_step, t`
  - Apply the filter function in order to determine possible floe pairings 
  - Apply the matching function to produce unique pairs of floes
  - Update the trajectories to include the newly paired floes
  - Add all unmatched floes as heads for new trajectories.

# Arguments
- `props::Vector{DataFrame}`: A vector of DataFrames, each containing ice floe properties for a single day. Each DataFrame must have the following columns:
    - "area"
    - "min_row"
    - "min_col"
    - "max_row"
    - "max_col"
    - "row_centroid"
    - "col_centroid"
    - "convex_area"
    - "major_axis_length"
    - "minor_axis_length"
    - "orientation"
    - "perimeter"
    - "mask": 2D boolean array
    - "passtime": A timestamp for the floe
    - "psi": the psi-s curve for the floe
    - "uuid": a universally unique identifier for each segmented floe
- `filter_function`: A function that uses
- `candidate_matching_settings`: settings for area mismatch and psi-s shape correlation. See `IceFloeTracker.candidate_matching_settings` for sample values.

# Returns
A DataFrame with the above columns, plus extra columns:
- `area_mismatch` and `corr`, which are the area mismatch and correlation between a floe and the one that preceeds it in the trajectory. 
- `head_uuid`, the floe which was best matched by this floe.
- Trajectories are identified by: 
  - a unique identifier `ID` and the 
  - UUID of the trajectory, `trajectory_uuid`.

Note: the props dataframes are modified in place.
"""
function floe_tracker(props::Vector{DataFrame}, filter_function, matching_function; minimum_area=100, maximum_area=90e3, maximum_time_step=Day(2))

    # dmw: give users option to copy props rather than modify in place?
    floe_size_filter = filter(
        r -> r.area >= minimum_area && r.area <= maximum_area
    )
    props .= floe_size_filter.(props)

    # Start_new_trajectory adds head_uuid and trajectory_uuid columns to props
    # The starting trajectories are just the floes visible and large enough on day 1.
    trajectories = props[1]
    _start_new_trajectory!(trajectories)

    for candidates in props[2:end]
        current_time_step = candidates[1, :passtime]
        trajectory_heads = _get_trajectory_heads(trajectories, current_time_step, maximum_time_step)

        candidate_pairs = []
        for floe in eachrow(trajectory_heads)
            append!(candidate_pairs, eachrow(filter_function(floe, candidates)))
        end
        
        # tbd: double check handling for initial and final floe uuid. 
        matched_pairs = DataFrame(candidate_pairs) |> matching_function

        # Get unmatched floes in day 2 (iterations > 2)
        matched_uuids = matched_pairs.uuid
        unmatched = filter((f) -> !(f.uuid in matched_uuids), candidates)
        _start_new_trajectory!(unmatched)

        # Attach new matches and unmatched floes to trajectories
        trajectories = vcat(trajectories, new_matches, unmatched)
    end
    trajectories = _drop_short_trajectories(trajectories, :trajectory_uuid)
    DataFrames.sort!(trajectories, [:trajectory_uuid, :passtime])
    _add_integer_id!(trajectories, :trajectory_uuid, :ID)
    # Move ID columns to the front
    select!(trajectories, :ID, :trajectory_uuid, :head_uuid, :uuid, :)
    return trajectories
end

# helper functions: all these should start with _ and should be defined in this file

function _start_new_trajectory!(floes::DataFrame)
    floes[!, :head_uuid] .= missing
    floes[!, :trajectory_uuid] .= [_uuid() for _ in eachrow(floes)]
    return floes
end

# TODO: replace hardcoded requirement to have the time variable be "passtime", e.g. allowing use of "time" or "observation_time" instead
"""
    get_trajectory_heads(pairs)

Return the last row (most recent member) of each group (trajectory) in `pairs` as a dataframe.

This is used for getting the initial floe properties for the next day in search for new pairs.
""" 
function _get_trajectory_heads(
    pairs::T, current_time_step, maximum_time_step; group_col=:trajectory_uuid, order_col=:passtime
) where {T<:AbstractDataFrame}
    gdf = groupby(pairs, group_col)
    heads = combine(gdf, x -> last(sort(x, order_col)))
    heads[:, :elapsed_time] = current_time_step .- heads[:, :order_col]
    subset!(heads, [:elapsed_time] => r -> r .<= maximum_time_step)
    select!(heads, Not(:elapsed_time))
    return heads
end

"""
    drop_trajectories_length1(trajectories::DataFrame, col::Symbol=:ID)

Drop trajectories with only one floe.

# Arguments
- `trajectories`: dataframe containing floe trajectories.
- `col`: column name for the floe ID.
"""
function _drop_short_trajectories(trajectories::DataFrame, col::Symbol=:ID; min_length=2)
    trajectories = filter(
        :count => x -> x >= min_length, transform(groupby(trajectories, col), nrow => :count)
    )
    select!(trajectories, Not("count"))
    return trajectories
end

"""
    _add_integer_id!(df, col, new)

For distinct values in the column `col` of `df`, add a new column `new` to be consecutive integers starting from 1.
"""
function _add_integer_id!(df::AbstractDataFrame, col::Symbol, new::Symbol)
    ids = unique(df[!, col])
    _map = Dict(ids .=> 1:length(ids))
    transform!(df, col => ByRow(x -> _map[x]) => new)
    return nothing
end

# """
#     find_floe_matches(
#     tracked,
#     candidate_props,
#     condition_thresholds,
#     mc_thresholds
# )

# Find matches for floes in `tracked` from floes in  `candidate_props`.

# # Arguments
# - `tracked`: dataframe containing floe trajectories.
# - `candidate_props`: dataframe containing floe candidate properties.
# - `candidate_filter_settings`: thresholds for deciding whether to match floe `i` from tracked to floe j from `candidate_props`
# - `candidate_matching_settings`: thresholds for area mismatch and psi-s shape correlation
# """
# function find_floe_matches(
#     tracked::T, candidate_props::T, candidate_filter_settings, candidate_matching_settings
# ) where {T<:AbstractDataFrame}
#     matches = []
#     for floe1 in eachrow(tracked), floe2 in eachrow(candidate_props)
#         Δt = floe2.passtime - floe1.passtime
#         ratios, conditions, dist = compute_ratios_conditions(
#             floe1, floe2, Δt, candidate_filter_settings
#         )

#         if callmatchcorr(conditions)
#             (area_mismatch, corr) = matchcorr(
#                 floe1.mask, floe2.mask, Δt; candidate_matching_settings.comp...
#             )
#             if isfloegoodmatch(conditions, candidate_matching_settings.goodness, area_mismatch, corr)
#                 @debug "** Found a good match for ", floe1.uuid, "<=", floe2.uuid
#                 push!(
#                     matches,
#                     (;
#                         Δt,
#                         measures=(; ratios..., area_mismatch, corr),
#                         conditions,
#                         dist,
#                         floe1,
#                         floe2,
#                     ),
#                 )
#             end
#         end
#     end

#     remaining_matches_df = DataFrame(matches)
#     best_matches = []

#     for floe2 in eachrow(
#         sort(candidate_props, :area; rev=true),  # prioritize larger floes
#     )
#         matches_involving_floe2_df = filter((r) -> r.floe2 == floe2, remaining_matches_df)
#         nrow(matches_involving_floe2_df) == 0 && continue
#         best_match = matches_involving_floe2_df[1, :]
#         measures_df = DataFrame(matches_involving_floe2_df.measures)
#         best_match_idx = getidxmostminimumeverything(measures_df) # TODO: Update this to take travel distance into account; potentially add weighing function
#         best_match = matches_involving_floe2_df[best_match_idx, :]
#         push!(
#             best_matches,
#             (;
#                 best_match.floe2...,
#                 area_mismatch=best_match.measures.area_mismatch,
#                 corr=best_match.measures.corr,
#                 trajectory_uuid=best_match.floe1.trajectory_uuid,
#                 head_uuid=best_match.floe1.uuid,
#             ),
#         )
#         # Filter out from the remaining matches any cases involving the matched floe1
#         remaining_matches_df = filter(
#             (r) -> !(r.floe1 === best_match.floe1), remaining_matches_df
#         )
#         # Filter out from the remaining matches any cases involving the matched floe2
#         remaining_matches_df = filter(
#             (r) -> !(r.floe2 === best_match.floe2), remaining_matches_df
#         )
#     end

#     best_matches_df = similar(tracked, 0)
#     append!(best_matches_df, best_matches; promote=true)
#     return best_matches_df
# end
