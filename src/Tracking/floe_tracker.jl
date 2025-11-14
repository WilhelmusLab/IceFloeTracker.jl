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
- `props::Vector{DataFrame}`: A vector of DataFrames, each containing ice floe properties for a single observation time. Each DataFrame must have the following columns:
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
    - "mask": 2D boolean floe image cropped to the floe location
    - "passtime": A timestamp for the floe
    - "psi": the psi-s curve for the floe
    - "uuid": a universally unique identifier for each segmented floe
- `filter_function`: A function that accepts a `floe::DataFrameRow` and a `candidates::DataFrame` argument, and subsets the candidates dataframe to those rows that are possible matches for `floe`.
- `matching_function`: A function that takes the dataframe of candidate pairs and resolves conflicts to find at most one match for each floe.

# Returns
A DataFrame with the above columns, plus extra columns:
- columns added by the filter function, such as similarity measures
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
    init_idx = 1
    nrow(props[init_idx]) == 0 && begin
        while nrow(props[init_idx]) == 0 && init_idx < length(props)
            global init_idx += 1
        end
    end
    trajectories = props[init_idx]
    _start_new_trajectory!(trajectories)

    for candidates in props[2:end]
        # Note: assumes each property table comes from a single observation time!
        nrow(candidates) > 0 && continue
        trajectory_heads = _get_trajectory_heads(trajectories, candidates[1, :passtime], maximum_time_step)

        candidate_pairs = []
        for floe in eachrow(trajectory_heads)
            candidates_subset = deepcopy(candidates)
            filter_function(floe, candidates_subset)
            nrow(candidates_subset) > 0 && begin
                candidates_subset[!, :head_uuid] .= floe.uuid
                candidates_subset[!, :trajectory_uuid] .= floe.trajectory_uuid
                append!(candidate_pairs, eachrow(candidates_subset))
            end
        end

        # matching function will find best pairs (head_uuid, uuid)
        # and ensure that all pairs are unique
        matched_pairs = DataFrame(candidate_pairs) |> matching_function

        # Get unmatched floes in day 2 (iterations > 2)
        # This should handle the case where there are no trajectory heads available
        matched_uuids = matched_pairs.uuid
        unmatched = filter((f) -> !(f.uuid in matched_uuids), candidates)
        _start_new_trajectory!(unmatched)
        _update_cols_to_match!(unmatched, matched_pairs)
        _update_cols_to_match!(trajectories, matched_pairs)

        # Attach new matches and unmatched floes to trajectories
        trajectories = vcat(trajectories, matched_pairs, unmatched)
    end
    trajectories = _drop_short_trajectories(trajectories, :trajectory_uuid)
    DataFrames.sort!(trajectories, [:trajectory_uuid, :passtime])
    _add_integer_id!(trajectories, :trajectory_uuid, :ID)
    # Move ID columns to the front
    select!(trajectories, :ID, :trajectory_uuid, :head_uuid, :uuid, :)
    return trajectories
end

# helper functions: all these should start with _ and should be defined in this file
"""
    _start_new_trajectory!(floes)

Initialize trajectory by adding empty `head_uuid` column and adding uuids to each row for `trajectory_uuid`.
"""
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
    heads[:, :elapsed_time] = current_time_step .- heads[:, order_col]
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

function _update_cols_to_match!(target::DataFrame, source::DataFrame; fill_value=missing)
    missing_cols = [c for c ∈ names(source) if c ∉ names(target)]
    for c in missing_cols
        target[!, c] .= fill_value
    end
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
