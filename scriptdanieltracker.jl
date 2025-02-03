# using Pkg
# Pkg.activate(".")
# Pkg.instantiate()
using IceFloeTracker: imshow, getidxmostminimumeverything, MatchedPairs, buildψs, corr, deleteallbut!
using Images
using Revise

# Optional for viewing psi-s curves:
using Plots
using Random

Random.seed!(123)
# Image timestamps can be found via the Satellite Overpass Information Tool (SOIT)
# or acquired manually

begin
    im1 = load("data/006-baffin_bay-20220530-aqua-labeled_floes.png")
    im2 = load("data/006-baffin_bay-20220530-terra-labeled_floes.png")

    # im1 = load("data/081-barents_kara_seas-20046708-terra-labeled_floes.png");
    # im2 = load("data/081-barents_kara_seas-20046708-aqua-labeled_floes.png");
    # passtimes = Vector([DateTime(2020, 6, 28, 22, 24), DateTime(2020, 6, 28, 22, 42)])




    passtimes = Vector([DateTime(2022, 5, 30, 15, 28), DateTime(2022, 5, 30, 16, 44)])
    init_images = [im1, im2]
    # Images.mosaicview(init_images, ncol=2, npad=10, fillvalue=1)

    labeled_images = [IceFloeTracker.label_components(image) for image in init_images]
    props = [IceFloeTracker.regionprops_table(image) for image in labeled_images]

    for p in props
        p.idx = 1:nrow(p)
    end

    props[1].idx
    props[2].idx

    # Filter out floes with area less than `floe_area_threshold` pixels
    # Too-small of objects can cause issues
    begin
        floe_area_threshold = 50
        for (i, prop) in enumerate(props)
            props[i] = prop[prop[:, :area].>=floe_area_threshold, :]
            sort!(props[i], :area, rev=true)
        end
    end


    IceFloeTracker.addfloemasks!(props, labeled_images) # Note the "!" indicates that the objects are modified in place
    IceFloeTracker.addψs!(props)

    # the tracker uses the elapsed time to determine the thr
    IceFloeTracker.add_passtimes!(props, passtimes)

    # adduuid gives each item in the props dataframes a unique identifier
    IceFloeTracker.adduuid!(props)

    # Set thresholds
    search_thresholds = (dt=(30.0, 100.0, 1300.0), dist=(200, 250, 300))
    large_floe_settings = (
        area=1000,
        arearatio=0.28,
        majaxisratio=0.10,
        minaxisratio=0.12,
        convexarearatio=0.15,
    )
    small_floe_settings = (
        area=1000,
        arearatio=0.18,
        majaxisratio=0.1,
        minaxisratio=0.15,
        convexarearatio=0.2,
    )
    condition_thresholds = (search_thresholds, small_floe_settings, large_floe_settings)
    mc_thresholds = (
        goodness=(small_floe_area=0.18, large_floe_area=0.236, corr=0.68), comp=(mxrot=10, sz=16)
    )
end




foo = Ref{Any}()
bar = Ref{Any}()
zoo = Ref{Any}()
tracker_results = IceFloeTracker.long_tracker(props, condition_thresholds, mc_thresholds)


ratp2 = hcat(rat, p2, makeunique=true)



# g = groupby(ratp2, :uuid)

# result = combine(groupby(ratp2, :uuid),
#                  g -> @view g[getidxmostminimumeverything(g[!, rationames]), :])
# result[!, [:uuid, :idx, :area, :area_mismatch, :corr]]



original_matched_pairs = bar[]
mp = deepcopy(original_matched_pairs)
perm = sortperm(mp.props2.uuid)
p2 = mp.props2[perm, :]
p1 = mp.props1[perm, :]
rat = mp.ratios[perm, :]
dist = mp.dist[perm]
matched_pairs = IceFloeTracker.MatchedPairs(p1, p2, rat, dist)

p1 = original_matched_pairs.props1
p1.idx_matches = 1:nrow(p1)
p2.idx_matches = 1:nrow(p1)
p2 = original_matched_pairs.props2
names(p2)
p2[!,[:idx, :idx_matches]]
rat = original_matched_pairs.ratios
rationames = names(rat)

matches = hcat(p1,p2, makeunique=true)[!, [:idx_matches, :uuid, :uuid_1]]
matches.cnt = [mp[uuid] for uuid in matches.uuid_1]
g = groupby(matches, :uuid_1)
mp = Dict(g.uuid_1[1] => nrow(g) for g in g)
matches = sort(matches, [:cnt, :uuid_1], rev=true)
matches[matches.cnt .> 1, :][!, [:uuid, :uuid_1, :cnt]]
matches[!, [:uuid, :uuid_1, :cnt]]

p2[!, [:uuid, :idx, :idx_matches,:area]]

ratp2 = hcat(rat, p2, makeunique=true)

g = groupby(ratp2, :uuid)

result = combine(groupby(ratp2, :uuid),
                 g -> @view g[getidxmostminimumeverything(g[!, rationames]), :])
result[!, [:uuid, :idx, :idx_matches, :area, :area_mismatch, :corr]]

filtered_df_other = filter(row -> row.idx_matches in result.idx_matches, p1)

filtered_df_other[!, [:uuid, :idx, :idx_matches, ]]



rng = 1:end
hcat(p1, p2, makeunique=true)[!, [:idx_matches, :uuid, :uuid_1]]
hcat(filtered_df_other, result, makeunique=true)[!, [:idx_matches,:idx_matches_1, :uuid, :uuid_1]]
























begin
    m1, m2, ratios, dist = foo[].matched_pairs.props1, foo[].matched_pairs.props2, foo[].matched_pairs.ratios, foo[].matched_pairs.dist
    m2ratios = hcat(m2, ratios, makeunique=true)
    m2ratios.idx = 1:nrow(m2ratios)
end

f1, f2 = (m1.mask[1], m2.mask[1])
_psi = buildψs.([f1, f2])
c = corr(_psi...)

begin
    # Group by uuid and select the row with max corr
    result = combine(groupby(m2ratios, :uuid),
        x -> @view x[argmax(x.corr), :])
end


idx = 6
idx1 = result[idx, :idx]
result[idx, :mask] |> imshow
m1[idx1, :mask] |> imshow

getidxmostminimumeverything(ratios)





toupdate = m2ratios[:, [:uuid, :idx, :area, :area_mismatch, :corr]]
get_best_match(g) = g.idx[argmax(g.corr)]

collisions = groupby(m2ratios, :uuid)
# collisions = [g for g in collisions if nrow(g) > 1]

g = collisions[2]
g[!, [:uuid, :idx, :area, :area_mismatch, :corr]]

bestmatches = [get_best_match(g) for g in collisions]

m1[bestmatches, :]
m2[bestmatches, :]
ratios[bestmatches, :]
dist[bestmatches, :] |> vec
mymatches = IceFloeTracker.MatchedPairs(m1[bestmatches, :], m2[bestmatches, :], ratios[bestmatches, :], dist[bestmatches, :] |> vec)

props1, props2 = props

IceFloeTracker.deletematched!((props1, props2), mymatches)