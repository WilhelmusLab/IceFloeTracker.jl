"""
    addlatlon(pairedfloesdf::DataFrame, refimage::AbstractString)

Add columns `latitude`, `longitude`, and pixel coordinates `x`, `y` to `pairedfloesdf`.

# Arguments
- `pairedfloesdf`: dataframe containing floe tracking data.
- `refimage`: path to reference image.
"""
function addlatlon!(pairedfloesdf::DataFrame, refimage::AbstractString)
    latlondata = latlon(refimage)
    colstodrop = [:row_centroid, :col_centroid, :min_row, :min_col, :max_row, :max_col]
    converttounits!(pairedfloesdf, latlondata, colstodrop)
    return nothing
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
    return nothing
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

function _pairfloes(
    segmented_imgs::Vector{<:FloeLabelsImage},
    props::Vector{DataFrame},
    passtimes::Vector{DateTime},
    condition_thresholds,
    mc_thresholds,
)
    dt = diff(passtimes) ./ Minute(1)

    sort_floes_by_area!(props)

    # Assign a unique ID to each floe in each image
    adduuid!(props)

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
                        props1[r, :], props2[s, :], Δt, condition_thresholds
                    )

                    if callmatchcorr(conditions)
                        (area_mismatch, corr) = matchcorr(
                            props1.mask[r], props2.mask[s], Δt; mc_thresholds.comp...
                        )

                        if isfloegoodmatch(
                            conditions, mc_thresholds.goodness, area_mismatch, corr
                        )
                            appendrows!(
                                matching_floes,
                                props2[s, :],
                                (ratios..., area_mismatch, corr),
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

    # Concatenate horizontally props1, props2, tracked, and add dist as the last column for each item in _pairs
    _pairs = [
        hcat(
            hcat(p.props1, p.props2; makeunique=true),
            p.ratios[:, ["area_mismatch", "corr"]],
        ) for p in _pairs
    ]

    # Make a dict with keys in _pairs[i].props2.uuid and values in _pairs[i-1].props1.uuid
    mappings = [Dict(zip(p.uuid_1, p.uuid)) for p in _pairs]

    # Convert mappings to functions
    funcsfrommappings = [x -> get(mapping, x, x) for mapping in mappings]

    # Compose functions in reverse order to push uuids forward
    mapuuid = foldr((f, g) -> x -> f(g(x)), funcsfrommappings)

    # Apply mapuuid to uuid_1 in each set of props in _pairs => get consolidated uuids
    [prop.uuid_0 = mapuuid.(prop.uuid_1) for prop in _pairs]

    # Reshape _pairs to a long df
    propsvert = vcat(_pairs...)
    DataFrames.sort!(propsvert, [:uuid_0, :passtime])
    rightcolnames = vcat(
        [
            name for name in names(propsvert) if
            all([!(name in ["uuid_1", "psi_1", "mask_1"]), endswith(name, "_1")])
        ],
        ["uuid_0"],
    )
    leftcolnames = [split(name, "_1")[1] for name in rightcolnames]
    matchcolnames = ["area_mismatch", "corr", "uuid_0", "passtime", "passtime_1"]

    leftdf = propsvert[:, leftcolnames]
    rightdf = propsvert[:, rightcolnames]
    matchdf = propsvert[:, matchcolnames]
    rename!(rightdf, Dict(zip(rightcolnames, leftcolnames)))

    _pairs = vcat(leftdf, rightdf)

    # sort by uuid_0, passtime and keep unique rows
    _pairs = unique(DataFrames.sort!(_pairs, [:uuid_0, :passtime]))

    _pairs = leftjoin(_pairs, matchdf; on=[:uuid_0, :passtime])
    DataFrames.sort!(_pairs, [:uuid_0, :passtime])

    # create mapping from uuids to index as ID
    uuids = unique(_pairs.uuid_0)
    uuid2index = Dict(uuid => i for (i, uuid) in enumerate(uuids))
    _pairs.ID = [uuid2index[uuid] for uuid in _pairs.uuid_0]
    _pairs = _pairs[:, [name for name in names(_pairs) if name != "uuid_0"]]
    return _pairs
end

"""
    pairfloes(
    segmented_imgs::Vector{BitMatrix},
    props::Vector{DataFrame},
    passtimes::Vector{DateTime},
    latlonrefimage::AbstractString,
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

### Arguments
- `segmented_imgs`: array of images with segmented floes
- `props`: array of dataframes containing floe properties
- `passtimes`: array of `DateTime` objects containing the time of the image in which the floes were captured
- `latlonrefimage`: path to geotiff reference image for getting latitude and longitude of floe centroids
- `condition_thresholds`: namedtuple of thresholds (each a named tuple) for deciding whether to match floe `i` from day `k` to floe j from day `k+1`
   (see sample values `IceFloeTracker.condition_thresholds`)


- `mc_thresholds`: thresholds for area mismatch and psi-s shape correlation

Returns a dataframe containing the following columns:
- `ID`: unique ID for each floe pairing.
- `passtime`: time of the image in which the floes were captured.
- `area`: area of the floe in sq. kilometers
- `convex_area`: area of the convex hull of the floe in sq. kilometers
- `major_axis_length`: length of the major axis of the floe in kilometers
- `minor_axis_length`: length of the minor axis of the floe in kilometers
- `orientation`: angle between the major axis and the x-axis in radians
- `perimeter`: perimeter of the floe in kilometers
- `latitude`: latitude of the floe centroid
- `longitude`: longitude of the floe centroid
- `x`: x-coordinate of the floe centroid
- `y`: y-coordinate of the floe centroid
- `area_mismatch`: area mismatch between the two floes in row_i and row_i+1 after registration
- `corr`: psi-s shape correlation between the two floes in row_i and row_i+1
"""
function pairfloes(
    segmented_imgs::Vector{<:FloeLabelsImage},
    props::Vector{DataFrame},
    passtimes::Vector{DateTime},
    latlonrefimage::AbstractString,
    condition_thresholds,
    mc_thresholds,
)
    _pairs = _pairfloes(
        segmented_imgs, props, passtimes, condition_thresholds, mc_thresholds
    )
    addlatlon!(_pairs, latlonrefimage)

    cols = [
        :ID,
        :passtime,
        :area,
        :convex_area,
        :major_axis_length,
        :minor_axis_length,
        :orientation,
        :perimeter,
        :latitude,
        :longitude,
        :x,
        :y,
        :area_mismatch,
        :corr,
    ]
    _pairs = _pairs[:, cols]
    return _pairs
end
