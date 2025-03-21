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
    begin # 0th iteration: pair floes in day 1 and day 2 and add unmatched floes to _pairs
        props1, props2 = props[1:2]
        matched_pairs0 = find_floe_matches(
            props1, props2, condition_thresholds, mc_thresholds, 2
        )
        @show props1[!, Main.cols]
        @show props2[!, Main.cols]
        @show matched_pairs0
        # @assert false
        Main.foo = (p1=matched_pairs0.props1[!, Main.cols], p2=matched_pairs0.props2[!, Main.cols])

        # Get unmatched floes from day 1/2
        unmatched1 = get_unmatched(props1, matched_pairs0.props1)
        @info "Unmatched floes in day 1:"
        @show unmatched1[!, Main.cols]
        unmatched2 = get_unmatched(props2, matched_pairs0.props2)
        @info "Unmatched floes in day 2:"
        @show unmatched2[!, Main.cols]
        unmatched = vcat(unmatched1, unmatched2)
        @info "Unmatched floes in day 1 and 2:"
        @show unmatched[!, Main.cols]
        consolidated_matched_pairs = consolidate_matched_pairs(matched_pairs0)

        # Get _pairs: preliminary matched and unmatched floes
        trajectories = vcat(consolidated_matched_pairs, unmatched)
        trajectories[:, [:uuid, :passtime, :area_mismatch, :corr]]
    end
    Main.tafter2 = trajectories

    begin # Start 3:end iterations
        for i in 3:length(props)
            @info "Processing Day $i"
            trajectory_heads = get_trajectory_heads(trajectories)
            candidate_props = props[i]
            Main.foo = (t=trajectories, h=trajectory_heads)
            # @assert false
            begin # Check trajectories heads
                tg = groupby(trajectory_heads, :uuid)
                tcounts = [nrow(g) for g in tg]
                mxcounts = maximum(tcounts)
                @assert mxcounts == 1
                @info "Trajectory heads"
                @show sort(trajectory_heads[!, Main.cols], "uuid")
            end

            floes_in_day_i = nrow(candidate_props)
            new_matches = IceFloeTracker.find_floe_matches(
                trajectory_heads, candidate_props, condition_thresholds, mc_thresholds, i
            )

            # Get unmatched floes in "day-2" (props[i]) (iterations > 2)
            unmatched2 = get_unmatched(props[i], new_matches.props2)
            @info "Unmatched floes in day $i:"
            @show unmatched2[!, Main.cols]
            new_pairs = IceFloeTracker.get_matches(new_matches)
            @show new_pairs[!, Main.cols]
            @assert nrow(new_pairs) + nrow(unmatched2) == floes_in_day_i

            # Attach new matches and unmatched floes to trajectories
            Main.foo = (
                trajectories=trajectories[!, Main.cols],
                new_pairs=new_pairs[!, Main.cols],
                unmatched2=unmatched2[!, Main.cols],
                new_matches=new_matches,
            )

            trajectories = vcat(trajectories, new_pairs, unmatched2)

            DataFrames.sort!(trajectories, [:uuid, :passtime])
            _swap_last_values!(trajectories)

            # Check trajectories counts
            tg = groupby(trajectories, :uuid)
            tcounts = [nrow(g) for g in tg]
            @show tcounts
            mxcounts = maximum(tcounts)
            argmaxcounts = argmax(tcounts)
            @show mxcounts argmaxcounts
            @show tg[argmaxcounts][!, Main.cols]
            @assert mxcounts <= i

            if i == 4
                @info "Terminating after Day $i"
                # @assert false
                break
            end

            # @show groupby(new_pairs[!, Main.cols], :uuid)
            # if i == 4
            #     @assert false
            # end
        end
    end
    trajectories = IceFloeTracker.drop_trajectories_length1(trajectories, :uuid)
    IceFloeTracker.reset_id!(trajectories, :uuid)
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
    heads::T, candidate_props, condition_thresholds, mc_thresholds, i
) where {T<:AbstractDataFrame}
    props1 = deepcopy(heads)
    props2 = deepcopy(candidate_props)
    # Main.foo = (props1=props1[!, Main.cols], props2=props2[!, Main.cols])
    # i == 4 && @assert false
    i == 4 && @show heads[!, Main.cols]
    # @assert false
    match_total = MatchedPairs(props2)

    while_rounds = 1
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
                # @assert false
                bestmatchdata = getbestmatchdata(best_match_idx, r, props1, matching_floes) # might be copying data unnecessarily
                i == 4 && @info "Processing day $i"
                # @assert false
                if i == 4 && r > 2
                    @info "From find_floe_matches, bestmatchdata"
                    p1 = matched_pairs.props1[!, Main.cols]
                    p2 = matched_pairs.props2[!, Main.cols]
                    @show p1
                    @show p2
                    @show bestmatchdata
                    Main.foo = (p1=p1, p2=p2, bestmatchdata=bestmatchdata)
                    # @assert false
                end
                addmatch!(matched_pairs, bestmatchdata, i, r)
                @debug "Matched pairs" matched_pairs

                # @show matched_pairs.props1[!, Main.cols]
                # @show matched_pairs.props2[!, Main.cols]
                Main.foo = (
                    props1=bestmatchdata.props1[Main.cols],
                    props2=bestmatchdata.props2[Main.cols],
                )
                # @assert false
                # @show matched_pairs.ratios[!, Main.cols]
                # @show matched_pairs.dist[!, Main.cols]
            end
        end # of for r = 1:nrow(props1)

        # exit while loop if there are no more floes to match
        @debug "Matched pairs" matched_pairs
        isempty(matched_pairs) && break

        #= Resolve collisions:
        Are there floes in day k+1 paired with more than one
        floe in day k? If so, keep the best matching pair and remove all others. =#
        if i == 4 && while_rounds == 1
            @info "Matched pairs before removing collisions"
            @info "while round: $while_rounds day: $i"
            bar = sort(
                hcat(matched_pairs.props1, matched_pairs.props2; makeunique=true), :uuid
            )
            Main.mp_prior_collision_removal = bar
            println(sort(bar[!, [:uuid, :_label, :uuid_1, :_label_1]], :uuid))
            # @assert false
        end

        matched_pairs = remove_collisions(matched_pairs)
        # if i == 4
        #     Main.mp_post_collision_removal = matched_pairs
        Main.check_matched_pairs(matched_pairs)
        # @assert false
        # end
        @info "Getting ready to delete matched pairs from props1 and props2"
        countofmatches = nrow(matched_pairs.props1)
        @info "There are $countofmatches matched pairs"
        # @assert false
        @show nrow(heads)
        @show nrow(props1)
        @show nrow(props2)
        # @assert false

        deletematched!((props1, props2), matched_pairs)

        if i == 4 && while_rounds == 1
            @info "After deleting matched pairs from props1 and props2"
            @show nrow(props1)
            @show nrow(props2)
            @info "Unmatched floes in day $i so far:"
            @show props1[!, Main.cols]
            @show matched_pairs.props1[!, Main.cols]
            @show matched_pairs.props2[!, Main.cols]
            # @assert false
            update!(match_total, matched_pairs)
        end
        while_rounds += 1
    end # of while loop
    Main.foo = (match_total=match_total,)
    # i == 4 && @assert false
    @info "total while rounds: $while_rounds"
    # @assert false
    return match_total
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
