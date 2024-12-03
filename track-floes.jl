HOME = "." # path to the root of the project two levels up

# Activate the environment
using Pkg
Pkg.activate(HOME)
Pkg.precompile()

using IceFloeTracker: pairfloes, deserialize, PaddedView, float64, mosaicview, Gray, sort_floes_by_area!, add_passtimes!, addfloemasks!, addψs!, getpropsday1day2, MatchedPairs, makeemptydffrom, compute_ratios_conditions, callmatchcorr, isfloegoodmatch, appendrows!, getidxmostminimumeverything, isnotnan, getbestmatchdata, addmatch!, resolvecollisions!, deletematched!, update!, Tracked
using DataFrames
using Dates
using Random
using StatsBase
using Test
imshow(x) = Gray.(x);


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

# Filter out floes with area less than 350 pixels
for (i, prop) in enumerate(_props)
    _props[i] = prop[prop[:, :area].>=500, :]
    sort!(_props[i], :area, rev=true)
end

pf(p) = pairfloes(_imgs, p, passtimes, latlonimgpth, condition_thresholds, mc_thresholds)
pf(p, _imgs) = pairfloes(_imgs, p, passtimes, latlonimgpth, condition_thresholds, mc_thresholds)

# Test case 1: New floes after the first day
# props_test_case1 = deepcopy(_props);
# delete!(props_test_case1[1], 4); # delete the smallest floe in day 1
# passtimes = _passtimes
# pairs_test_case1 = pf(props_test_case1)
# Four trajectories, three of length 3 and one of length 2
# Set([3 => 3, 2 => 1])
# @test countmap(pairs_test_case1[:, :ID]) |> values |> countmap |> Set == Set([3 => 3, 2 => 1])

begin
    # Test case 2: No new floes after the first day but some floes are missing
    props_test_case2 = deepcopy(_props)
    delete!(props_test_case2[2], 4)
    segmented_imgs = _imgs
    passtimes = _passtimes
    # delete!(props_test_case2[3], 4);
    # Expect four trajectories, three of length 3 and one of length 1
    # Outcome: three trajectories of length 3. Unpaired floe in day 1 left behind
    pairs_test_case2 = pairfloes(segmented_imgs, (props_test_case2), passtimes, latlonimgpth, condition_thresholds, mc_thresholds)
end

# begin
    # Get unmatched floes
    unmatched = IceFloeTracker.get_unmatched(props_test_case2[1], pairs_test_case2[1].props1)


    # Consolidate (horizontally) props1, props2, and ratios into a single data structure
    _pairs = [hcat(hcat(p.props1, p.props2, makeunique=true), p.ratios[:, ["area_mismatch", "corr"]]) for p in pairs_test_case2]
    # names(_pairs[1])



    # Reshape _pairs to a long df
    propsvert = vcat(_pairs...) # same as _pairs[1] as _pairs is a vector of DataFrames with one element
    DataFrames.sort!(propsvert, [:passtime])

    # Update uuid_1 to uuid to sort later by uuid and then by passtime
    propsvert.uuid_1 = propsvert.uuid;
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

    add_missing = ["area_mismatch", "corr"]
    [unmatched[!, n] = [missing for _ in 1:nrow(unmatched)] for n in add_missing]

    _pairs = vcat(_pairs, unmatched[:, names(_pairs)])

    gdf = groupby(_pairs, :uuid)
    baseprops = combine(gdf, last)[:, names(_pairs)]
    # baseprops[:, [:uuid, :passtime]]
# end
# vcat(_pairs, baseprops)

tracked = Tracked()
props = props_test_case2
dt = diff(passtimes) ./ Minute(1)
dayi = 2
props1, props2 = getpropsday1day2(props, dayi)

props1 = baseprops

# IceFloeTracker.addfloemasks!(prop, img)
Δt = dt[dayi]

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
foo = tracked.data

# begin
    # Get unmatched floes
    @info baseprops |> names
    @info foo[1].props1 |> names
    unmatched = IceFloeTracker.get_unmatched(baseprops, foo[1].props1)


    # Consolidate (horizontally) props1, props2, and ratios into a single data structure
    _pairs = [hcat(hcat(p.props1, p.props2, makeunique=true), p.ratios[:, ["area_mismatch", "corr"]]) for p in pairs_test_case2]
    # names(_pairs[1])



    # Reshape _pairs to a long df
    propsvert = vcat(_pairs...) # same as _pairs[1] as _pairs is a vector of DataFrames with one element
    DataFrames.sort!(propsvert, [:passtime])

    # Update uuid_1 to uuid to sort later by uuid and then by passtime
    propsvert.uuid_1 = propsvert.uuid;
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

    add_missing = ["area_mismatch", "corr"]
    [unmatched[!, n] = [missing for _ in 1:nrow(unmatched)] for n in add_missing]

    _pairs = vcat(_pairs, unmatched[:, names(_pairs)])

    gdf = groupby(_pairs, :uuid)
    baseprops = combine(gdf, last)[:, names(_pairs)]
    # baseprops[:, [:uuid, :passtime]]
# end