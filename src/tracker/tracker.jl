"""
    pair_floes(indir::String, condition_thresholds, mc_thresholds=(area3=0.18, area2=0.236, corr=0.68))

$(include("pair_floes_docstring.jl"))
"""
function pair_floes(
    indir::String, condition_thresholds, mc_thresholds=(area3=0.18, area2=0.236, corr=0.68)
)
    input_data = deserialize(joinpath(indir, "MB_tracker_inputs.dat"))
    properties = input_data["prop"]
    segmented_imgs = input_data["FLOE_LIBRARY"] # used in matchcorr
    delta_time = input_data["delta_t"]
    return newpairfloes(
        segmented_imgs, properties, delta_time, condition_thresholds, mc_thresholds
    )
end

function pairfloes(
    segmented_imgs::Vector{BitMatrix},
    props::Vector{DataFrame},
    dt::Vector{Int64},
    condition_thresholds,
    mc_thresholds,
)
    # Initialize container for props of matched pairs of floes, their similarity ratios, and their distances between their centroids
    tracked = Tracked()

    # Crop floes from the images using the bounding box data in `props`.
    addfloemasks!(props, segmented_imgs)

    numdays = length(imgs) - 1

    for dayi in 1:numdays
        props1, props2 = getpropsday1day2(props, dayi)
        Δt = dt[dayi]

        # Container for matches in dayi which will be used to populate tracked
        match_total = MatchedPairs(props1)
        while true # there are no more floes to match in props1
            # This rutine mutates both props1 and props2.

            # Container for props of matched floe pairs and their similarity ratios. Matches will be updated and added to match_total
            matched_pairs = MatchedPairs(props1)
            for r in 1:nrow(props1) # TODO: consider using eachrow(props1) to iterate over rows
                # 1. Collect preliminary matches for floe r in matching_floes
                matching_floes = makeemptydffrom(props1)
                # Threads.@threads 
                for s in 1:nrow(props2) # TODO: consider using eachrow(props2) to iterate over rows
                    ratios, conditions, dist = compute_ratios_conditions(
                        (props1, r), (props2, s), Δt, condition_thresholds
                    )

                    if callmatchcorr(conditions)
                        (area_under, corr) = matchcorr((r, props1), (s, props2), Δt) # TODO: build matchcorr

                        if isfloegoodmatch(conditions, mc_thresholds, area_under, corr)
                            # collect data for matching floe s 
                            appendrows!(
                                matching_floes,
                                props2[s, :],
                                (ratios..., area_under, 1 - corr),
                                s,
                                dist,
                            )
                        end
                    end
                end # of s for loop

                # 2. Find the best match for floe r
                best_match_idx = getidxmostminimumeverything(matching_floes.ratios)
                if isnotnan(best_match_idx)
                    bestmatchdata = getbestmatchdata(
                        best_match_idx, r, props1, matching_floes
                    ) # might be copying data unnecessarily

                    addmatch!(matched_pairs, bestmatchdata)
                end
            end # of for r = 1:nrow(props1)

            # exit while loop if there are no more floes to match
            isempty(matched_pairs) && break

            #= Resolve collisions:
            Are there floes in day k+1 paired with more than one
            floe in day k? If so, keep the best matching pair and remove all others. =#
            resolvecollisions!(matched_pairs)
            deletematched!((props1, props2), matched_pairs)
            update!(match_total, matched_pairs)
        end # of while loop
        update!(tracked, match_total)
    end
    sort!(tracked)
    return tracked
end
