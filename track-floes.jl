
HOME = "." # path to the root of the project two levels up

# Activate the environment
using Pkg
Pkg.activate(HOME)
Pkg.precompile()

using IceFloeTracker: pairfloes, deserialize, PaddedView, float64, mosaicview, Gray, sort_floes_by_area!, add_passtimes!, addfloemasks!, addψs!, getpropsday1day2, MatchedPairs, makeemptydffrom, compute_ratios_conditions, callmatchcorr, isfloegoodmatch, appendrows!, getidxmostminimumeverything, isnotnan, getbestmatchdata, addmatch!, resolvecollisions!, deletematched!, update!, Tracked, get_trajectory_heads, get_unmatched, _pairfloes, _swap_last_values!, get_dt, find_floe_matches
using DataFrames
using Dates
using Random
using StatsBase
using Test
imshow(x) = Gray.(x);

# begin

# Set thresholds
begin
    t1 = (dt=(30.0, 100.0, 1300.0), dist=(200, 250, 300))
    t2 = (
        area=1200,
        arearatio=0.28,
        majaxisratio=0.10,
        minaxisratio=0.12,
        convexarearatio=0.14,
    )
    t3 = (
        area=10_000,
        arearatio=0.18,
        majaxisratio=0.1,
        minaxisratio=0.15,
        convexarearatio=0.2,
    )
    condition_thresholds = (t1, t2, t3)
    mc_thresholds = (
        goodness=(area3=0.18, area2=0.236, corr=0.68), comp=(mxrot=10, sz=16)
    )
end

# Load data
begin
    pth = joinpath(HOME, "test", "test_inputs", "tracker")
    _floedata = deserialize(joinpath(pth, "tracker_test_data.dat"))
    _passtimes = deserialize(joinpath(pth, "passtimes.dat"))
    latlonimgpth = joinpath(HOME, "test", "test_inputs", "NE_Greenland_truecolor.2020162.aqua.250m.tiff")
    _props, _imgs = deepcopy(_floedata.props), deepcopy(_floedata.imgs)
    IceFloeTracker.addfloemasks!(_props, _imgs)
end

# Filter out floes with area less than `floe_area_threshold` pixels
floe_area_threshold = 400
for (i, prop) in enumerate(_props)
    _props[i] = prop[prop[:, :area].>=floe_area_threshold, :] # 500 working good
    sort!(_props[i], :area, rev=true)
end

# Prep data for _pairfloes. TODO: use as test case
begin
    # Unmatched floe in day 1, unmatched floe in day 2, and matches for every floe starting in day 3
    props_test_case2 = deepcopy(_props)
    delete!(props_test_case2[1], 1)
    delete!(props_test_case2[2], 5)
    segmented_imgs = _imgs
    passtimes = _passtimes
    # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2
end

begin # 0th iteration: pair floes in day 1 and day 2 and add unmatched floes to _pairs
    foo = _pairfloes(segmented_imgs, props_test_case2, passtimes, condition_thresholds, mc_thresholds)

    # Get unmatched floes from day 1
    unmatched1 = get_unmatched(props_test_case2[1], foo[1].props1)
    unmatched2 = get_unmatched(props_test_case2[2], foo[1].props2)
    unmatched = vcat(unmatched1, unmatched2)
end

begin # Pairs consolidation
    # Consolidate (horizontally) props1, props2, and ratios into a single data structure
    missingcols = [:area_mismatch, :corr]
    _pairs = [hcat(hcat(p.props1, p.props2, makeunique=true), p.ratios[:, missingcols]) for p in foo]
    # Reshape _pairs to a long df
    propsvert = vcat(_pairs...) # same as _pairs[1] as _pairs is a vector of DataFrames with one element
    DataFrames.sort!(propsvert, [:passtime])

    # Update uuid_1 to uuid to sort later by uuid and then by passtime
    propsvert.uuid_1 = propsvert.uuid
    rightcolnames = vcat([name for name in names(propsvert) if endswith(name, "_1")])
    leftcolnames = [split(name, "_1")[1] for name in rightcolnames]
    matchcolnames = ["uuid", "passtime", "area_mismatch", "corr"]
    leftdf = propsvert[:, leftcolnames]
    rightdf = propsvert[:, rightcolnames]
    rename!(rightdf, Dict(zip(rightcolnames, leftcolnames)))
    matchdf = propsvert[:, matchcolnames]
    _pairs = vcat(leftdf, rightdf)
    # sort by uuid, passtime and keep unique rows
    _pairs = DataFrames.sort!(_pairs, [:uuid, :passtime]) |> unique
    _pairs = leftjoin(_pairs, matchdf, on=[:uuid, :passtime])
    DataFrames.sort!(_pairs, [:uuid, :passtime])
end

# Update _pairs with unmatched floes
_pairs = vcat(_pairs, unmatched[:, names(_pairs)])


# Start 2:end iterations
trajectory_heads = get_trajectory_heads(_pairs)

begin # Set up next i+1 iteration with trajectory heads
    tracked = Tracked()
    props = props_test_case2
    dt = diff(passtimes) ./ Minute(1)
    dayi = 2
end

# TODO: Instead of using props1, use trajectory_heads in getpropsday1day2. Perhaps rename function.
_, props2 = getpropsday1day2(props, dayi)
props1 = trajectory_heads

# IceFloeTracker.addfloemasks!(prop, img)

# Container for matches in dayi which will be used to populate tracked
match_total = MatchedPairs(props2)
while true # there are no more floes to match in props1
    # This routine mutates both props1 and props2.

    # Container for props of matched floe pairs and their similarity ratios. Matches will be updated and added to match_total
    matched_pairs = MatchedPairs(props2)
    for r in 1:nrow(props1) # TODO: consider using eachrow(props1) to iterate over rows
        @info "Matching floe $r in day $dayi"
        # 1. Collect preliminary matches for floe r in matching_floes
        matching_floes = makeemptydffrom(props2)

        for s in 1:nrow(props2) # TODO: consider using eachrow(props2) to iterate over rows
            Δt = get_dt(props1, r, props2, s)
            @info "Considering floe $s in day $(dayi+1) for floe $r in day $dayi"
            ratios, conditions, dist = compute_ratios_conditions(
                (props1, r), (props2, s), Δt, condition_thresholds
            )

            if callmatchcorr(conditions)
                @info "Getting mismatch and correlation for floe $r in day $dayi and floe $s in day $(dayi+1)"
                (area_mismatch, corr) = matchcorr(
                    props1.mask[r], props2.mask[s], Δt; mc_thresholds.comp...
                )

                if isfloegoodmatch(
                    conditions, mc_thresholds.goodness, area_mismatch, corr
                )
                    @info "names in matching_floes" matching_floes |> names
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
        @info "Finding best match for floe $r in day $dayi"
        best_match_idx = getidxmostminimumeverything(matching_floes.ratios)
        @info "Best match index for floe $r in day $dayi: $best_match_idx"
        if isnotnan(best_match_idx)
            bestmatchdata = getbestmatchdata(
                best_match_idx, r, props1, matching_floes
            ) # might be copying data unnecessarily
            @info "names of bestmatchdata.props1" names(bestmatchdata.props1)
            @info "names of bestmatchdata.props2" names(bestmatchdata.props2)
            @info "names of matched_pairs.props1" matched_pairs.props1 |> names
            @assert names(bestmatchdata.props1) == names(matched_pairs.props1)
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

# end
sort!(tracked)
_foo = tracked.data
foo = IceFloeTracker.find_floe_matches(trajectory_heads, props[3], condition_thresholds, mc_thresholds)

begin # Get unmatched floes in day 2 (iterations > 2)
    unmatched2 = get_unmatched(props[3], foo.props2)
    @assert isempty(unmatched2)
end

# Apply flattening workflow to the tracked floes and update pairs
# begin
# 1. Get unmatched floes and add missing columns for join compatibility

# Test adding an unmatched floe in props2
begin
    fake_props2 = deepcopy(props[3])
    newrow = deepcopy(fake_props2[end, :])
    newrow.area = 5
    uuidfake = "zfakefloe123"
    newrow.uuid = uuidfake
    push!(fake_props2, newrow)
    unmatched_df = get_unmatched(fake_props2, foo.props2)
    @assert unmatched_df[1, :uuid] == uuidfake
    @assert "corr" in names(unmatched_df)
end


begin # Consolidate (horizontally) props1, props2, and ratios into a single data structure
    consolidated = [hcat(hcat(p.props1, p.props2, makeunique=true), p.ratios[:, missingcols]) for p in [foo]]
    propsvert = vcat(consolidated...) # same as newpairs[1] as newpairs is a vector of DataFrames with one element
    DataFrames.sort!(propsvert, [:area, :passtime]) # TODO: might not be necessary
    # Update uuid_1 to uuid to sort later by uuid and then by passtime
    propsvert.uuid_1 = propsvert.uuid
    # Fix column names in rightdf
    leftcolnames = names(propsvert)[1:16]
    rightcolnames = names(propsvert)[17:end]
    rightdf = propsvert[:, vcat(rightcolnames)]
    rename!(rightdf, Dict(zip(rightcolnames, leftcolnames)))

    _pairs_rightdf = vcat(_pairs, rightdf, unmatched2)
    DataFrames.sort!(_pairs_rightdf, [:uuid, :passtime])
    _swap_last_values!(_pairs_rightdf)
    # IceFloeTracker.reset_ids!(_pairs_rightdf)
end

# Repeat

# At the end make :uuid => :ID where ID is a unique identifier for each floe from 1 to length(unique(:uuid))