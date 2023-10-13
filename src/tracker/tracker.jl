"""
    pairfloes(
    segmented_imgs::Vector{BitMatrix},
    props::Vector{DataFrame},
    passtimes::Vector{DateTime},
    dt::Vector{Float64},
    condition_thresholds,
    mc_thresholds,
)

Pair floes in `props[k]` to floes in `props[k+1]` for `k=1:length(props)-1`.

The main steps of the algorithm are as follows:

1. Crop floes from `segmented_imgs` using bounding box data in `props`. Floes in the edges are removed.
2. For each floe_k_r in `props[k]`, compare to floe_k+1_s in `props[k+1]` by computing similarity ratios, set of `conditions`, and drift distance `dist`. If the conditions are met, compute the area mismatch `mm` and psi-s correlation `c` for this pair of floes. Pair these two floes if `mm` and `c` satisfy the thresholds in `mc_thresholds`.
3. If there are collisions (i.e. floe `s` in `props[k+1]` is paired to more than one floe in `props[k]`), then the floe in `props[k]` with the best match is paired to floe `s` in `props[k+1]`.
4. Drop paired floes from `props[k]` and `props[k+1]` and repeat steps 2 and 3 until there are no more floes to match in `props[k]`.
5. Repeat steps 2-4 for `k=2:length(props)-1`.

# Arguments
- `segmented_imgs`: array of images with segmented floes.
- `props`: array of dataframes containing floe properties.
- `passtimes`: array of `DateTime` objects containing the time of the image in which the floes were captured.
- `dt`: array of time elapsed between images in `segmented_imgs`.
- `condition_thresholds`: 3-tuple of thresholds (each a named tuple) for deciding whether to match floe `i` from day `k` to floe j from day `k+1`.
- `mc_thresholds`: thresholds for area mismatch and psi-s shape correlation.

Returns a tuple `(props, trackdata)` where `props` is a long dataframe containing floe ID's, passtimes, the original set of physical properties, and their masks and `trackdata` is a dataframe containing the floe tracking data.
"""
function pairfloes(
    segmented_imgs::Vector{BitMatrix},
    props::Vector{DataFrame},
    passtimes::Vector{DateTime},
    dt::Vector{Float64},
    condition_thresholds,
    mc_thresholds,
)
    sort_floes_by_area!(props)

    # Assign a unique ID to each floe in each image
    for (i, prop) in enumerate(props)
        props[i].uuid = [randstring(12) for _ in 1:nrow(prop)]
    end

    add_passtimes!(props, passtimes)

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
    _pairs = tracked.data

    # Make a dict with keys in _pairs[i].props2.uuid and values in _pairs[i-1].props1.uuid
    mappings = [Dict(pair.props2.uuid .=> pair.props1.uuid) for pair in _pairs]

    # Convert mappings to functions
    funcsfrommappings = [x -> get(mapping, x, x) for mapping in mappings]
    
    # Compose functions in reverse order to push uuids forward
    mapuuid = foldr((f, g) -> x -> f(g(x)), funcsfrommappings)

    for prop in props[2:end]
        prop.uuid = mapuuid.(prop.uuid)
    end

    # Collect all unique uuids in props[i] to label as simple ints starting from 1
    uuids = unique([uuid for prop in props for uuid in prop.uuid])
    
    # create mapping from uuids to index
    uuid2index = Dict(uuid => i for (i, uuid) in enumerate(uuids))

    # apply the uuid2index mapping to props
    for prop in props
        prop.uuid .= [uuid2index[uuid] for uuid in prop.uuid]
    end

    # Merge all props into one long DataFrame
    propsvert = vcat(props...)

    # rename uuid to ID
    rename!(propsvert, :uuid => :ID)

    # 2. Sort propsvert by uuid and then by passtime
    DataFrames.sort!(propsvert, [:ID, :passtime])

    # 3. Move ID, passtime columns to the front
    propsvert = propsvert[:, unique(["ID", "passtime", names(propsvert)...])]

    return (props = propsvert[:, names(propsvert)[1:15]], trackdata = _pairs)
end


"""
    add_passtimes!(props, passtimes)

Add a column `passtime` to each DataFrame in `props` containing the time of the image in which the floes were captured.

# Arguments
- `props`: array of DataFrames containing floe properties.
- `passtimes`: array of `DateTime` objects containing the time of the image in which the floes were captured.

"""
function add_passtimes!(props, passtimes)
    for (i, passtime) in enumerate(passtimes)
        props[i].passtime .= passtime
    end
    nothing
end

"""
    sort_floes_by_area!(props)

Sort floes in `props` by area in descending order.
"""
function sort_floes_by_area!(props)
    for prop in props
        # sort by area in descending order
        DataFrames.sort!(prop, :area; rev=true)
        nothing
    end
end
