"""
    long_tracker(props, condition_thresholds, mc_thresholds)

Track ice floes over multiple observations.

Trajectories are built as follows:
- Assume the floes detected in observation 1 are trajectories of length 1.
- For each subsequent observation:
  - Determine the latest observation for each trajectory – these are the "current trajectory heads".
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
- `candidate_filter_settings`: namedtuple of settings and functions for reducing the number of possible matches. See `IceFloeTracker.candidate_filter_settings` for sample values.
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
function long_tracker(props::Vector{DataFrame}, candidate_filter_settings, candidate_matching_settings)
    # dmw: this can be generalized to have a global area filter setting in the candidate_filter_settings
    filter_out_small_floes = filter(
        r -> r.area >= candidate_filter_settings.small_floe_settings.minimumarea
    )
    # dmw: we should warn users that props gets mutated in place, or give an option to avoid mutating it
    props .= filter_out_small_floes.(props)

    # The starting trajectories are just the floes visible and large enough on day 1.
    trajectories = props[1]
    _start_new_trajectory!(trajectories)

    for prop in props[2:end]
        trajectory_heads = get_trajectory_heads(trajectories)


        new_matches = find_floe_matches(
            trajectory_heads, prop, candidate_filter_settings, candidate_matching_settings
        )

        # Get unmatched floes in day 2 (iterations > 2)
        matched_uuids = new_matches.uuid
        unmatched = filter((f) -> !(f.uuid in matched_uuids), prop)
        _start_new_trajectory!(unmatched)

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

function _start_new_trajectory!(floes::DataFrame)
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
- `candidate_filter_settings`: thresholds for deciding whether to match floe `i` from tracked to floe j from `candidate_props`
- `candidate_matching_settings`: thresholds for area mismatch and psi-s shape correlation
"""
function find_floe_matches(
    tracked::T, candidate_props::T, candidate_filter_settings, candidate_matching_settings
) where {T<:AbstractDataFrame}
    matches = []
    
    for floe1 in eachrow(tracked)
        floe1_candidate_props = deepcopy(candidate_props)

        #### Step 1: Restrict to plausible distances
        # TODO: Update to use x/y when we've added that to the property tables.
        # TODO: Rename dist to euclidean_distance to avoid confusion with built-ins
        floe1_candidate_props[!, :Δt] = floe1_candidate_props.passtime .- floe1.passtime
        floe1_candidate_props[!, :Δx] = [dist(getcentroid(floe1), getcentroid(floe2); candidate_filter_settings.resolution) 
                                            for floe2 in eachrow(floe1_candidate_props)]        
        floe1_candidate_props .= filter(r -> candidate_filter_settings.time_distance_threshold_function(r.Δx, r.Δt), floe1_candidate_props)
        
    
        #### Step 2: Compute and check geometric ratios
        # TODO: Add geometric_thresholds function to the candidate filter settings
        if nrow(floe1_candidate_props) > 0
           
            # need to get ratio names in order to build dataframe columns
            # better to merge? or to loop?
            ratios = compute_ratios(floe1, floe1_candidate_props) # won't work yet!
           
           
            # Add function (or functor) to compute whether the geometric ratios are acceptable.
            # large_floe_condition = get_large_floe_condition(floe1.area, ratios, candidate_filter_settings)
            # small_floe_condition = get_small_floe_condition(floe1.area, ratios, candidate_filter_settings)
           floe1_candidate_props .= filter(r -> candidate_filter_settings.geometric_thresholds(r), floe1_candidate_props)
        end

        #### Step 3: Psi-s correlation
        if nrow(flow1_candidate_props) > 0

            # broadcast this function across the psi-s curves
            r, r_low, r_high = round(normalized_maximum_crosscorrelation(_psi...), digits=3)

            # filter the candidates -- transform this from the pseudocode to real code
            floe1_candidate_props .= filter(r -> candidate_filter_settings.minimum_psi_s_corr in (r_low, r_high), floe1_candidate_props)
        end
            
        if nrow(flow1_candidate_props) > 0

            # Named tuple with results. Do we want to do psi-s and correlation and registration all at once?
            area_mismatch, corr, rot, corr_ci, mismatch_ci, rotation_ci = matchcorr(
                floe1.mask, floe2.mask, Δt; candidate_matching_settings.comp...
            )

           end

        end

            
            if isfloegoodmatch(conditions, candidate_matching_settings.goodness, area_mismatch, corr)
                @debug "** Found a good match for ", floe1.uuid, "<=", floe2.uuid
                push!(
                    matches,
                    (;
                        Δt,
                        measures=(; ratios..., area_mismatch, corr), # TODO: Add confidence intervals here
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

    for floe2 in eachrow(
        sort(candidate_props, :area; rev=true),  # prioritize larger floes
    )
        matches_involving_floe2_df = filter((r) -> r.floe2 == floe2, remaining_matches_df)
        nrow(matches_involving_floe2_df) == 0 && continue
        best_match = matches_involving_floe2_df[1, :]
        measures_df = DataFrame(matches_involving_floe2_df.measures)
        best_match_idx = getidxmostminimumeverything(measures_df) # TODO: Update this to take travel distance into account; potentially add weighing function
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
