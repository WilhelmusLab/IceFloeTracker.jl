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

# Load data
pth = joinpath(HOME, "test", "test_inputs", "tracker")
_floedata = deserialize(joinpath(pth, "tracker_test_data.dat"))
_passtimes = deserialize(joinpath(pth, "passtimes.dat"))
latlonimgpth = joinpath(HOME, "test", "test_inputs", "NE_Greenland_truecolor.2020162.aqua.250m.tiff")
_props, _imgs = deepcopy(_floedata.props), deepcopy(_floedata.imgs);


# Filter out floes with area less than 350 pixels
for (i, prop) in enumerate(_props)
    _props[i] = prop[prop[:, :area].>=500, :]
    sort!(_props[i], :area, rev=true)
end

pf(p) = pairfloes(_imgs, p, passtimes, latlonimgpth, condition_thresholds, mc_thresholds)
pf(p, _imgs) = pairfloes(_imgs, p, passtimes, latlonimgpth, condition_thresholds, mc_thresholds)

# Test case 1: New floes after the first day
props_test_case1 = deepcopy(_props);
delete!(props_test_case1[1], 4); # delete the smallest floe in day 1
passtimes = _passtimes
pairs_test_case1 = pf(props_test_case1)
# Four trajectories, three of length 3 and one of length 2
# Set([3 => 3, 2 => 1])
# @test countmap(pairs_test_case1[:, :ID]) |> values |> countmap |> Set == Set([3 => 3, 2 => 1])



# Test case 2: No new floes after the first day but some floes are missing
props_test_case2 = deepcopy(_props)[1:end-1];
delete!(props_test_case2[2], 4);
segmented_imgs = _imgs[1:end-1]
passtimes = _passtimes[1:end-1]
# delete!(props_test_case2[3], 4);
# Expect four trajectories, three of length 3 and one of length 1
# Outcome: three trajectories of length 3. Unpaired floe in day 1 left behind
pairs_test_case2 = pairfloes(segmented_imgs, (props_test_case2), passtimes, latlonimgpth, condition_thresholds, mc_thresholds);

# Get unmatched floes
unmatched = IceFloeTracker.get_unmatched(props_test_case2[1], pairs_test_case2[1].props1);


# Consolidate (horizontally) props1, props2, and ratios into a single data structure
_pairs = [hcat(hcat(p.props1, p.props2, makeunique=true), p.ratios[:, ["area_mismatch", "corr"]]) for p in pairs_test_case2]
# names(_pairs[1])



# Reshape _pairs to a long df
propsvert = vcat(_pairs...) # same as _pairs[1] as _pairs is a vector of DataFrames with one element
DataFrames.sort!(propsvert, [:passtime])

# Update uuid_1 to uuid to sort later by uuid and then by passtime
propsvert.uuid_1 = propsvert.uuid
rightcolnames = vcat([name for name in names(propsvert) if all([!(name in ["psi_1", "mask_1"]), endswith(name, "_1")])])
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




