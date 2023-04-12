"""
    pair_floes(indir::String, condition_thresholds, mc_thresholds=(area3=0.18, area2=0.236, corr=0.68))

$(include("pair_floes_docstring.jl"))
"""
function pair_floes(
    indir::String, condition_thresholds, mc_thresholds=(area3=0.18, area2=0.236, corr=0.68)
)
    input_data = deserialize(joinpath(indir, "MB_tracker_inputs.dat"))
    # Keys:
    #     "FLOE_LIBRARY"  [x]
    #     "old_data_expected" [x]
    #     "prop"     [x]
    #     "delta_t"  [x]

    properties = input_data["prop"]
    floe_library = input_data["FLOE_LIBRARY"] # used in matchcorr
    numdays = size(properties, 2) - 1

    # Initialize container for props of matched pairs of floes, their similarity ratios, and their distances between their centroids
    tracked = Tracked()

    # Traverse each pair of succesive "days" of floe properties. Grab the properties of floe r in day k. Find all preliminary matches for floe r in day k+1 and keep the best match. If floe s in day k+1 is paired with more than one floe in day k, keep the best matching pair and remove all others. 
    for dayi in 1:numdays
        props_day1, props_day2 = getpropsday1day2(properties, dayi)
        delta_time = input_data["delta_t"][dayi]

        # Container for matches in dayi which will be used to populate tracked
        match_total = MatchedPairs(props_day1)

        while true # there are no more floes to match in props_day1
            # This rutine mutates props_day1 and props_day2.

            # Container for props of matched floe pairs and their similarity ratios. Matches will be updated and added to match_total
            matched_pairs = MatchedPairs(props_day1)

            for r in 1:nrow(props_day1) # TODO: consider using eachrow(props_day1) to iterate over rows

                # 1. Collect preliminary matches for floe r in matching_floes
                matching_floes = makeemptydffrom(props_day1)

                Threads.@threads for s in 1:nrow(props_day2) # TODO: consider using eachrow(props_day2) to iterate over rows
                    ratios, conditions, dist = compute_ratios_conditions(
                        (props_day1, r), (props_day2, s), delta_time, condition_thresholds
                    )

                    if callmatchcorr(conditions)
                        (area_under, corr) = matchcorr(r, s, dayi, delta_time, floe_library) # TODO: build matchcorr

                        if isfloegoodmatch(conditions, mc_thresholds, area_under, corr)
                            # collect collect data for matching floe s 
                            appendrow!(
                                matching_floes,
                                props_day2[s, :],
                                (ratios..., area_under, 1 - corr),
                                s, # is this really needed?
                                dist,
                            )
                        end
                    end
                end # of s for loop

                # 2. Find the best match for floe r
                best_match_idx = getidxmostminimumeverything(matching_floes.ratios)
                if isnotnan(best_match_idx)
                    bestmatchdata = getbestmatchdata(
                        best_match_idx, r, props_day1, matching_floes
                    ) # might be copying data unnecessarily
                    addmatch!(matched_pairs, bestmatchdata)
                end
            end # of for r = 1:nrow(props_day1)

            # exit while loop if there are no more floes to match
            isempty(matched_pairs) && break

            #= Resolve collisions:
            Are there floes in day k+1 paired with more than one
            floe in day k? If so, keep the best matching pair and remove all others. =#
            resolvecollisions!(matched_pairs)

            deletematched!((props_day1, props_day2), matched_pairs)
            update!(match_total, matched_pairs)
        end # of while loop
        update!(tracked, match_total)
    end
    sort!(tracked)
    return tracked
end