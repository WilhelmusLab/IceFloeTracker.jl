"""
    long_tracker(props, condition_thresholds, mc_thresholds)

Track ice floes over multiple days.

Trajectories are built in two steps:
0. Get pairs of floes in day 1 and day 2. Any unmatched floes, in both day 1 and day 2, become the "heads" of their respective trajectories.
1. For each subsequent day, find pairs of floes for the current trajectory heads. Again, any unmatched floe in the new prop table starts a new trajectory.

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
A DataFrame with the above columns, plus two extra columns, "area_mismatch" and "corr", which are the area mismatch and correlation between a floe and the one that follows it in the trajectory. Trajectories are identified by a unique identifier.
"""
function long_tracker(props::Vector{DataFrame}, condition_thresholds, mc_thresholds)
    begin # Filter out floes with area less than `small_floe_minimum_area` pixels
        small_floe_minimum_area = condition_thresholds.small_floe_settings.minimumarea
        for (i, prop) in enumerate(props)
            props[i] = prop[prop[:, :area] .>= small_floe_minimum_area, :]
            DataFrames.sort!(props[i], :area; rev=true)
        end
    end

    # The starting trajectories are just the floes visible on day 1.
    trajectories = props[1]
    trajectories[!, :head_uuid] .= trajectories[:, :uuid]
    trajectories[!, :area_mismatch] .= missing
    trajectories[!, :corr] .= missing

    for prop in props[2:end]
        trajectory_heads = get_trajectory_heads(trajectories)
        new_matches = find_floe_matches(
            trajectory_heads, prop, condition_thresholds, mc_thresholds
        )

        # Get unmatched floes in day 2 (iterations > 2)
        matched_uuids = new_matches.uuid
        unmatched = filter((f) -> !(f.uuid in matched_uuids), prop)
        unmatched[!, :head_uuid] = unmatched[:, :uuid]  # unmatched floes start new trajectories
        unmatched[!, :area_mismatch] .= missing
        unmatched[!, :corr] .= missing

        # Attach new matches and unmatched floes to trajectories
        trajectories = vcat(trajectories, new_matches, unmatched)
    end
    trajectories = drop_trajectories_length1(trajectories, :head_uuid)
    DataFrames.sort!(trajectories, [:head_uuid, :passtime])
    _add_integer_id!(trajectories, :head_uuid, :ID)
    # Move ID columns to the front
    select!(trajectories, :ID, :head_uuid, :uuid, :)
    return trajectories
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
                head_uuid=best_match.floe1.head_uuid,
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
