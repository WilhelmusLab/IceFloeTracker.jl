"""
    long_tracker(props, condition_thresholds, mc_thresholds)

Track ice floes over multiple observations.

Trajectories are built as follows:
- Assume the floes detected in observation 1 are trajectories of length 1.
- For each subsequent observation:
  - Determine the newest observation for each trajectory – these are the "current trajectory heads".
  - Find matches between the the current trajectory heads and the new observed floes, extending those trajectories.
  - Any unmatched floe in an observation is added as a new trajectory starting point.

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
    - "mask": 2D array of booleans
    - "passtime": A timestamp for the floe
    - "psi": the psi-s curve for the floe
    - "uuid": a universally unique identifier for each segmented floe
- `condition_thresholds`: namedtuple of thresholds for deciding whether to match floe `i` from day `k` to floe j from day `k+1`. See `IceFloeTracker.condition_thresholds` for sample values.
- `mc_thresholds`: thresholds for area mismatch and psi-s shape correlation. See `IceFloeTracker.mc_thresholds` for sample values.

# Returns
A DataFrame with the above columns, plus extra columns:
- `area_mismatch` and `corr`, which are the area mismatch and correlation between a floe and the one that preceeds it in the trajectory. 
- `head_uuid`, the floe which was best matched by this floe.
- Trajectories are identified by: 
  - a unique identifier `ID` and the 
  - UUID of the trajectory, `trajectory_uuid`.
"""
function long_tracker(props::Vector{DataFrame}, condition_thresholds, mc_thresholds)
    filter_out_small_floes = filter(
        r -> r.area >= condition_thresholds.small_floe_settings.minimumarea
    )

    # The starting trajectories are just the floes visible and large enough on day 1.
    trajectories = props[1]
    trajectories = filter_out_small_floes(trajectories)
    # Order by largest first
    sort!(trajectories, :area; rev=true)

    _new_trajectory!(trajectories)

    for prop in props[2:end]
        trajectory_heads = get_trajectory_heads(trajectories)

        prop = filter_out_small_floes(prop)
        sort!(prop, :area; rev=true)

        new_matches = find_floe_matches(
            trajectory_heads, prop, condition_thresholds, mc_thresholds
        )

        # Get unmatched floes in day 2 (iterations > 2)
        matched_uuids = new_matches.uuid
        unmatched = filter((f) -> !(f.uuid in matched_uuids), prop)
        _new_trajectory!(unmatched)

        # Attach new matches and unmatched floes to trajectories
        trajectories = vcat(trajectories, new_matches, unmatched)
    end
    trajectories = drop_trajectories_length1(trajectories, :trajectory_uuid)
    DataFrames.sort!(trajectories, [:trajectory_uuid, :passtime])
    _add_integer_id!(trajectories, :trajectory_uuid, :ID)
    # Move ID columns to the front
    select!(trajectories, :ID, :trajectory_uuid, :head_uuid, :uuid, :)
    return trajectories
end

function _new_trajectory!(floes::DataFrame)
    floes[!, :head_uuid] .= missing
    floes[!, :trajectory_uuid] .= [_uuid() for _ in eachrow(floes)]
    floes[!, :area_mismatch] .= missing
    floes[!, :corr] .= missing
    return floes
end

"""
    find_floe_matches(
    tracked,
    candidate_props,
    condition_thresholds,
    mc_thresholds
)

Find matches for floes in `tracked` from floes in  `candidate_props`.

# Arguments
- `tracked`: dataframe containing floe trajectories.
- `candidate_props`: dataframe containing floe candidate properties.
- `condition_thresholds`: thresholds for deciding whether to match floe `i` from tracked to floe j from `candidate_props`
- `mc_thresholds`: thresholds for area mismatch and psi-s shape correlation
"""
function find_floe_matches(
    tracked::T, candidate_props::T, condition_thresholds, mc_thresholds
) where {T<:AbstractDataFrame}
    matches = []
    for floe1 in eachrow(tracked), floe2 in eachrow(candidate_props)
        Δt = get_dt(floe1, floe2)
        ratios, conditions, dist = compute_ratios_conditions(
            floe1, floe2, Δt, condition_thresholds
        )

        if callmatchcorr(conditions)
            (area_mismatch, corr) = matchcorr(
                floe1.mask, floe2.mask, Δt; mc_thresholds.comp...
            )
            if isfloegoodmatch(conditions, mc_thresholds.goodness, area_mismatch, corr)
                @debug "** Found a good match for ", floe1.uuid, "<=", floe2.uuid
                push!(
                    matches,
                    (;
                        Δt,
                        measures=(; ratios..., area_mismatch, corr),
                        conditions,
                        dist,
                        floe1,
                        floe2,
                    ),
                )
            end
        end
    end

    remaining_matches_df = DataFrame(matches)
    best_matches = []

    for floe2 in eachrow(candidate_props)  # leave at most one match for each 
        matches_involving_floe2_df = filter((r) -> r.floe2 == floe2, remaining_matches_df)
        if nrow(matches_involving_floe2_df) == 0
            continue
        end
        best_match = matches_involving_floe2_df[1, :]
        measures_df = DataFrame(matches_involving_floe2_df.measures)
        best_match_idx = getidxmostminimumeverything(measures_df)
        best_match = matches_involving_floe2_df[best_match_idx, :]
        push!(
            best_matches,
            (;
                best_match.floe2...,
                area_mismatch=best_match.measures.area_mismatch,
                corr=best_match.measures.corr,
                trajectory_uuid=best_match.floe1.trajectory_uuid,
                head_uuid=best_match.floe1.uuid,
            ),
        )
        # Filter out from the remaining matches any cases involving the matched floe1
        remaining_matches_df = filter(
            (r) -> !(r.floe1 === best_match.floe1), remaining_matches_df
        )
        # Filter out from the remaining matches any cases involving the matched floe2
        remaining_matches_df = filter(
            (r) -> !(r.floe2 === best_match.floe2), remaining_matches_df
        )
    end

    best_matches_df = similar(tracked, 0)
    append!(best_matches_df, best_matches; promote=true)
    return best_matches_df
end

# Sample values for condition_thresholds
small_floe_minimum_area = 400
large_floe_minimum_area = 1200
_dt = (30.0, 100.0, 1300.0)
_dist = (200, 250, 300)
search_thresholds = (dt=_dt, dist=_dist)

large_floe_settings = (
    minimumarea=large_floe_minimum_area,
    arearatio=0.28,
    majaxisratio=0.10,
    minaxisratio=0.12,
    convexarearatio=0.14,
)

small_floe_settings = (
    minimumarea=small_floe_minimum_area,
    arearatio=0.18,
    majaxisratio=0.1,
    minaxisratio=0.15,
    convexarearatio=0.2,
)
condition_thresholds = (
    search_thresholds=search_thresholds,
    small_floe_settings=small_floe_settings,
    large_floe_settings=large_floe_settings,
)

mc_thresholds = (
    goodness=(small_floe_area=0.18, large_floe_area=0.236, corr=0.68),
    comp=(mxrot=10, sz=16),
)
