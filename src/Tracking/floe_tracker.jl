import DataFrames: DataFrames, DataFrame, AbstractDataFrame, eachrow, select!, subset!
import Images: SegmentedImage
import Dates: DateTime, Period, Day
import ..Segmentation: regionprops_table

abstract type AbstractTracker end

"""
    FloeTracker(    
        filter_function::AbstractFloeFilterFunction
        matching_function::AbstractFloeMatchingFunction
        minimum_area::Real = 100
        maximum_area::Real = 90e3
        maximum_time_step::Period = Day(2)
    )
    FloeTracker()(
        segmented_images::Vector{<:Union{SegmentedImage,Matrix{Int64}}},
        image_times::Vector{DateTime}
    )
    
Track ice floes over multiple observations.

The FloeTracker functor initializes the floe tracking function by setting the `filter_function` 
(see  [`FilterFunction`](@ref)) and the `matching_function` (e.g., ['MinimumWeightMatchingFunction'](@ref)),
and basic filter parameters for the area range and maximum time step. 

Trajectories are built as follows:  
- Assume the floes detected in observation 1 are trajectories of length 1.
- For each subsequent observation at time `t`:
- Determine the latest observation for each trajectory -- these are the "current trajectory heads".
- Select the subset of trajectory heads observed within the window `maximum_time_step, t`
- Apply the filter function in order to determine possible floe pairings 
- Apply the matching function to produce unique pairs of floes
- Update the trajectories to include the newly paired floes
- Add all unmatched floes as heads for new trajectories.

## Arguments
- `filter_function`: Function that reduces the number of rows in a DataFrame and adds columns for tests.
- `matching_function`: Function that decides between conflicting matches to form a unique matching.
- `minimum_area`: Minimum object size in pixels.
- `maximum_area`: Maximum object size in pixels.
- `maximum_time_step`: Maximum time step between observations using Dates period type.
- `segmented_images`: Vector of SegmentedImage objects or labeled indexmaps
- `image_times`: Vector of DateTimes of the same length as `segmented_images`

## Examples
Using the default functions, initialize as:
```jldoctest; setup = :(using IceFloeTracker)
julia> tracker = FloeTracker(filter_function=FilterFunction(), matching_function=MinimumWeightMatchingFunction())
```

Once the tracker is defined, it can be run on a list of either SegmentedImages (or labeled image indexmaps) and
a list of corresponding observation times. As a simple toy example, we can use labeled blocks. We can't expect 
shapes to have consistent labels across images before tracking, so we'll intentionally mislabel them. We also need
to provide observation times. The default thresholds are time-step dependent, so we choose a short time step. We also
need to set the minimum size to accommodate the toy example.

```jldoctest; setup = :(using IceFloeTracker, Dates)
julia> tracker = FloeTracker(
        filter_function=FilterFunction(),
        matching_function=MinimumWeightMatchingFunction(),
        minimum_area=1)

julia> A = zeros(Int, 13, 16); A[2:6, 2:6] .= 1; A[4:8, 7:10] .= 2; A[10:12,13:15] .= 3; A[10:12,3:6] .= 4;
julia> B = zeros(Int, 13, 16); B[1:5, 1:5] .= 2; B[5:9, 7:10] .= 3; B[10:12,12:14] .= 4; B[10:12,2:5] .= 1;
julia> times = [DateTime("2025-05-01T11:00"), DateTime("2025-05-01T13:00")]

julia> tracked_floes = tracker([A, B], times)

julia> tracked_floes[:, ["ID", "label", "passtime"]]
8×3 DataFrame
Row │ ID     label  passtime            
    │ Int64  Int64  DateTime            
─────┼───────────────────────────────────
1 │     1      2  2025-05-01T11:00:00
2 │     1      3  2025-05-01T13:00:00
3 │     2      4  2025-05-01T11:00:00
4 │     2      1  2025-05-01T13:00:00
5 │     3      1  2025-05-01T11:00:00
6 │     3      2  2025-05-01T13:00:00
7 │     4      3  2025-05-01T11:00:00
8 │     4      4  2025-05-01T13:00:00
```

Note that the tracker has assigned each object a unique ID, and that the objects are linked correctly:
1=>2, 2=>3, 3=>4, and 4=>1.
"""
@kwdef struct FloeTracker <: AbstractTracker
    filter_function::AbstractFloeFilterFunction
    matching_function::AbstractFloeMatchingFunction
    minimum_area::Real = 100
    maximum_area::Real = 90e3
    maximum_time_step::Period = Day(2)
end
# TODO: Add the minimum area and maximum area to the FilterFunction
# TODO: Add method to functor to get list of needed columns from the filter functions or struct (e.g., if doing cross correlation)

function (t::FloeTracker)(
    segmented_images::Vector{<:Union{SegmentedImage,Matrix{Int64}}},
    image_times::Vector{DateTime},
)
    props = regionprops_table.(segmented_images)
    add_uuids!.(props)
    !issorted(image_times) && @warn "Passtimes are not in ascending order."
    add_passtimes!.(props, image_times) # TODO: Change function name to image_times
    add_floemasks!.(props, segmented_images)
    add_ψs!.(props)

    tracking_results = floe_tracker(
        props,
        t.filter_function,
        t.matching_function;
        minimum_area=t.minimum_area,
        maximum_area=t.maximum_area,
        maximum_time_step=t.maximum_time_step,
    )

    return tracking_results
end

# TODO: Make this one an internal function, but extend FloeTracker functor to allow list of props
"""
    floe_tracker(props; filter_function, matching_function, minimum_floe_size, maximum_floe_size, maximum_time_step)

Lower-level function for tracking from an already-existing property table. See the  [`FloeTracker`](@ref) function for comparison.

## Arguments
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

## Returns
A DataFrame with the above columns, plus extra columns:
- columns added by the filter function, such as similarity measures
- `head_uuid`, the floe which was best matched by this floe.
- Trajectories are identified by: 
  - a unique identifier `ID` and the 
  - UUID of the trajectory, `trajectory_uuid`.

Note: the props dataframes are modified in place.
"""
function floe_tracker(
    props::Vector{DataFrame},
    filter_function,
    matching_function;
    minimum_area=100,
    maximum_area=90e3,
    maximum_time_step=Day(2),
)

    # dmw: give users option to copy props rather than modify in place?
    floe_size_filter = filter(r -> r.area >= minimum_area && r.area <= maximum_area)
    props .= floe_size_filter.(props)

    # Start_new_trajectory adds head_uuid and trajectory_uuid columns to props
    # The starting trajectories are just the floes visible and large enough on day 1.
    init_idx = 1
    nrow(props[init_idx]) == 0 && begin
        while nrow(props[init_idx]) == 0 && init_idx < length(props)
            init_idx += 1
        end
    end
    trajectories = props[init_idx]
    _start_new_trajectory!(trajectories)

    for candidates in props[2:end]
        # Note: assumes each property table comes from a single observation time!
        nrow(candidates) > 0 && begin
            trajectory_heads = _get_trajectory_heads(
                trajectories, candidates[1, :passtime], maximum_time_step
            )

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
    pairs::T,
    current_time_step,
    maximum_time_step;
    group_col=:trajectory_uuid,
    order_col=:passtime,
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

## Arguments
- `trajectories`: dataframe containing floe trajectories.
- `col`: column name for the floe ID.
"""
function _drop_short_trajectories(trajectories::DataFrame, col::Symbol=:ID; min_length=2)
    trajectories = filter(
        :count => x -> x >= min_length,
        transform(groupby(trajectories, col), nrow => :count),
    )
    select!(trajectories, Not("count"))
    return trajectories
end

function _update_cols_to_match!(target::DataFrame, source::DataFrame; fill_value=missing)
    missing_cols = [c for c in names(source) if c ∉ names(target)]
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
