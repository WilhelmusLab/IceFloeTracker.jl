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
- `condition_thresholds`: 3-tuple of thresholds (each a named tuple) for deciding whether to match floe `i` from day `k` to floe j from day `k+1`
- `mc_thresholds`: thresholds for area mismatch and psi-s shape correlation

# Returns
A DataFrame with the above columns, plus two extra columns, "area_mismatch" and "corr", which are the area mismatch and correlation between a floe and the one that follows it in the trajectory. Trajectories are identified by a unique identifier, "uuid".
"""
function long_tracker(props::Vector{DataFrame}, condition_thresholds, mc_thresholds)
    begin # 0th iteration: pair floes in day 1 and day 2 and add unmatched floes to _pairs
        props1, props2 = props[1:2]
        matched_pairs0 = find_floe_matches(props1, props2, condition_thresholds, mc_thresholds)

        # Get unmatched floes from day 1/2
        unmatched1 = get_unmatched(props1, matched_pairs0.props1)
        unmatched2 = get_unmatched(props2, matched_pairs0.props2)
        unmatched = vcat(unmatched1, unmatched2)
        consolidated_matched_pairs = consolidate_matched_pairs(matched_pairs0)

        # Get _pairs: preliminary matched and unmatched floes
        trajectories = vcat(consolidated_matched_pairs, unmatched)
        trajectories[:, [:uuid, :passtime, :area_mismatch, :corr]]
    end

    begin # Start 3:end iterations
        for i in 3:length(props)
            trajectory_heads = get_trajectory_heads(trajectories)
            new_pairs = IceFloeTracker.find_floe_matches(trajectory_heads, props[i], condition_thresholds, mc_thresholds)
            # Get unmatched floes in day 2 (iterations > 2)
            unmatched2 = get_unmatched(props[i], new_pairs.props2)
            new_pairs = IceFloeTracker.get_matches(new_pairs)

            # Attach new matches and unmatched floes to trajectories
            trajectories = vcat(trajectories, new_pairs, unmatched2)
            DataFrames.sort!(trajectories, [:uuid, :passtime])
            _swap_last_values!(trajectories)
        end
    end
    IceFloeTracker.reset_id!(trajectories)
    trajectories.ID = trajectories.uuid
    # list the uuid in the leftmost column
    cols = [col for col in names(trajectories) if col ∉ ["ID", "uuid"]]
    return trajectories[!, ["ID", cols...]]
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
    tracked::T,
    candidate_props::T,
    condition_thresholds,
    mc_thresholds
) where {T<:AbstractDataFrame}
    props1 = deepcopy(tracked)
    props2 = deepcopy(candidate_props)
    match_total = MatchedPairs(props2)
    while true # there are no more floes to match in props1

        # This routine mutates both props1 and props2.
        # Get preliminary matches for floe r in props1 from floes in props2

        # Container for props of matched floe pairs and their similarity ratios. Matches will be updated and added to match_total
        matched_pairs = MatchedPairs(props2)
        for r in 1:nrow(props1) # TODO: consider using eachrow(props1) to iterate over rows
            # 1. Collect preliminary matches for floe r in matching_floes
            matching_floes = makeemptydffrom(props2)

            for s in 1:nrow(props2) # TODO: consider using eachrow(props2) to iterate over rows
                Δt = get_dt(props1, r, props2, s)
                @debug "Considering floe 2:$s for floe 1:$r"
                ratios, conditions, dist = compute_ratios_conditions(
                    (props1, r), (props2, s), Δt, condition_thresholds
                )

                if callmatchcorr(conditions)
                    @debug "Getting mismatch and correlation for floe 1:$r and floe 2:$s"
                    (area_mismatch, corr) = matchcorr(
                        props1.mask[r], props2.mask[s], Δt; mc_thresholds.comp...
                    )

                    if isfloegoodmatch(
                        conditions, mc_thresholds.goodness, area_mismatch, corr
                    )
                        @debug "** Found a good match for floe 1:$r => 2:$s"
                        appendrows!(
                            matching_floes,
                            props2[s, :],
                            (ratios..., area_mismatch, corr),
                            s,
                            dist,
                        )
                        @debug "Matching floes" matching_floes
                    end
                end
            end # of s for loop

            # 2. Find the best match for floe r
            @debug "Finding best match for floe 1:$r"
            best_match_idx = getidxmostminimumeverything(matching_floes.ratios)
            @debug "Best match index for floe 1:$r: $best_match_idx"
            if isnotnan(best_match_idx)
                bestmatchdata = getbestmatchdata(
                    best_match_idx, r, props1, matching_floes
                ) # might be copying data unnecessarily
                addmatch!(matched_pairs, bestmatchdata)
                @debug "Matched pairs" matched_pairs
            end
        end # of for r = 1:nrow(props1)

        # exit while loop if there are no more floes to match
        @debug "Matched pairs" matched_pairs
        isempty(matched_pairs) && break

        #= Resolve collisions:
        Are there floes in day k+1 paired with more than one
        floe in day k? If so, keep the best matching pair and remove all others. =#

        matched_pairs = remove_collisions(matched_pairs)
        deletematched!((props1, props2), matched_pairs)
        update!(match_total, matched_pairs)
    end # of while loop
    return match_total
end
