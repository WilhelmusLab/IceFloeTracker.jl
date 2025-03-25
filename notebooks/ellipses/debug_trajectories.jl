begin
    using Revise
    using Test
    using IceFloeTracker
    using Serialization
    using DataFrames
    using IceFloeTracker: imshow, compute_ratios_conditions, get_dt

    function get_day_and_floe_num(label)
        # Use a regular expression to match the desired parts
        match_result = match(r"day_(\d+)-#(\d+)", label)

        day_number = match_result[1]  # The number between "day_" and "-"
        hash_number = match_result[2]  # The number after "#"

        return (daynum=day_number, hashnum=hash_number)
    end

    function Counter(collection)
        counts = Dict{eltype(collection),Int}()  # Create a dictionary to store counts
        for item in collection
            counts[item] = get(counts, item, 0) + 1  # Increment the count for each item
        end
        return counts
    end

    function check_matched_pairs(mp)
        @assert sort(collect(Main.Counter(mp.props1.uuid)); by=x -> x[2], rev=true)[1][2] ==
            1
    end

    wait_for_key(prompt=nothing) = (print(stdout, prompt); read(stdin, 1); nothing)
    cols = [:uuid, :_label]
end

begin # Load the data / set config
    tracker_data = deserialize("notebooks/ellipses/tracker-trajectories-bug.jls")
    IceFloeTracker.adduuid!(tracker_data.props, :_label)
    _show(v) =
        for i in v
            println(i)
        end
    # trajectories_ = IceFloeTracker.long_tracker(
    #     tracker_data.props, tracker_data.ct, tracker_data.thresholds
    # )

    ct = (
        search_thresholds=(dt=(30.0, 100.0, 1300.0), dist=(200, 250, 300)),
        small_floe_settings=(
            minimumarea=100,
            arearatio=0.18,
            majaxisratio=0.1,
            minaxisratio=0.15,
            convexarearatio=0.2,
        ),
        large_floe_settings=(
            minimumarea=1200,
            arearatio=0.28,
            majaxisratio=0.1,
            minaxisratio=0.12,
            convexarearatio=0.14,
        ),
    )
end

begin
    cols = [:uuid, :_label]
    tafter2 = Ref{Any}()
    foo = Ref{Any}()
    mp_prior_collision_removal = Ref{Any}()
    mp_post_collision_removal = Ref{Any}()
    trajectories = IceFloeTracker.long_tracker(
        tracker_data.props, ct, tracker_data.thresholds
    )
end

foo

foo[foo.uuid .== "loVjHZ7UHVGv", cols] |> println
foo.new_pairs[foo.new_pairs.uuid .== "loVjHZ7UHVGv", :] |> println
sort!(foo.new_matches.props1[!, cols], :uuid) |> println
foo.new_matches.props1[!, cols] |> println
# investigate foo.new_matches which is used to obtain new_pairs


trajectories[!, [:ID, :_label]] |> println

# begin
    println(sort(foo.trajectories[!, cols], cols))
    grps = groupby(foo.trajectories, "uuid")
    [println((i, nrow(g))) for (i, g) in enumerate(grps) if nrow(g) > 4]
# end
for i in 1:5
    println(grps[i][!, [:uuid, :_label]])
end
grps[1][!, [:uuid, :_label]]

# remove collisions workflow
pairs = mp_prior_collision_removal
nm1 = names(pairs.props1)
nm2 = ["$(n)_1" for n in nm1]
old_nmratios = names(pairs.ratios)
nmratios = ["$(n)_ratio" for n in old_nmratios]
rename!(pairs.ratios, nmratios)
pairsdf = hcat(
    pairs.props1, pairs.props2, pairs.ratios, DataFrame(; dist=pairs.dist); makeunique=true
)

pairs_post = hcat(
    mp_post_collision_removal.props1,
    mp_post_collision_removal.props2,
    mp_post_collision_removal.ratios,
    DataFrame(; dist=mp_post_collision_removal.dist);
    makeunique=true,
)
pairs_post = pairs_post[!, [:_label, :_label_1]]
println(sort(pairs_post[!, [:_label, :_label_1]], :_label))
println(sort(pairs_post[!, [:_label, :_label_1]], :_label_1))
check_matched_pairs(mp_post_collision_removal)

println(trajectories[!, [:ID, :_label]])
grps = groupby(trajectories, "ID")
[println((i, nrow(g))) for (i, g) in enumerate(grps) if nrow(g) > 5];
grps[15][!, [:ID, :_label]]

foo.props1[cols] #, foo.props2[!, cols]

foo.props1[!, cols]
lpairs = sort(collect(Main.Counter(foo.props1.uuid)); by=x -> x[2], rev=true)
println.(lpairs)

uuid = foo.trajectories.uuid[1]

foo[foo.uuid .== uuid, cols]

[foo.props1._label foo.props2._label]

r, s = 1, 2
imshow(foo.props1.mask[r])
imshow(foo.props2.mask[s])

foo.props1[r, :]
foo.props2[s, :]
delta_time = get_dt(foo.props1, r, foo.props2, s)
compute_ratios_conditions((foo.props1, r), (foo.props2, s), delta_time, ct)

imshow(foo.props1.mask[1])
imshow(foo.props2.mask[4])

meta = get_day_and_floe_num.(foo.uuid)
days = [m.daynum for m in meta]
floe_counts = Counter(days)

_show(Counter(foo.uuid))

trajectory_groups = groupby(foo, "uuid")
print(trajectory_groups[1])

# assert matched pairs has one floe per uuid
