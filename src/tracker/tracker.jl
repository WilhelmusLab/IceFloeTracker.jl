"""
    pairfloes(
    segmented_imgs::Vector{BitMatrix},
    props::Vector{DataFrame},
    dt::Vector{Int64},
    condition_thresholds,
    mc_thresholds,
)

Pair floes in `props[k]` to floes in `props[k+1]` for `k=1:length(props)-1` into an array of `MatchedPairs`.

The main steps of the algorithm are as follows:

1. Crop floes from `segmented_imgs` using bounding box data in `props`.
2. For each floe_k_r in `props[k]`, compare to floe_k+1_s in `props[k+1]` by computing similarity ratios, set of `conditions`, and drift distance `dist`. If the conditions are met, compute the area mismatch `mm` and psi-s correlation `c` for this pair of floes. Pair these two floes if `mm` and `c` satisfy the thresholds in `mc_thresholds`.
3. If there are collisions (i.e. floe `s` in `props[k+1]` is paired to more than one floe in `props[k]`), then the floe in `props[k]` with the best match is paired to floe `s` in `props[k+1]`.
4. Drop paired floes from `props[k]` and `props[k+1]` and repeat steps 2 and 3 until there are no more floes to match in `props[k]`.
5. Repeat steps 2-4 for `k=2:length(props)-1`.

# Arguments
- `segmented_imgs`: array of images with segmented floes.
- `props`: array of dataframes containing floe properties.
- `dt`: array of time elapsed between images in `segmented_imgs`.
- `condition_thresholds`: 3-tuple of thresholds (each a named tuple) for deciding whether to match floe `i` from day `k` to floe j from day `k+1`.
- `mc_thresholds`: thresholds for area mismatch and psi-s shape correlation.

Returns an array of `MatchedPairs` containing the properties of matched floe pairs, their similarity ratios, and their distances between their centroids.
"""
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
    addψs!(props)

    numdays = length(segmented_imgs) - 1

    for dayi in 1:numdays
        props1, props2 = getpropsday1day2(props, dayi)
        Δt = dt[dayi]

        # Container for matches in dayi which will be used to populate tracked
        match_total = MatchedPairs(props1)
        while true # there are no more floes to match in props1
            # This routine mutates both props1 and props2.

            # Container for props of matched floe pairs and their similarity ratios. Matches will be updated and added to match_total
            matched_pairs = MatchedPairs(props1)
            for r in 1:nrow(props1) # TODO: consider using eachrow(props1) to iterate over rows
                # 1. Collect preliminary matches for floe r in matching_floes
                matching_floes = makeemptydffrom(props1)
                for s in 1:nrow(props2) # TODO: consider using eachrow(props2) to iterate over rows
                    ratios, conditions, dist = compute_ratios_conditions(
                        (props1, r), (props2, s), Δt, condition_thresholds
                    )

                    if callmatchcorr(conditions)
                        (area_under, corr) = matchcorr(
                            props1.mask[r], props2.mask[s], Δt; mc_thresholds.comp...
                        )

                        if isfloegoodmatch(
                            conditions, mc_thresholds.goodness, area_under, corr
                        )
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
    return tracked.data
end
