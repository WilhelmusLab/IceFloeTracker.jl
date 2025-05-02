begin
    using Test
    using IceFloeTracker
    using Serialization
    using DataFrames
    using IceFloeTracker: imshow, compute_ratios_conditions, get_dt, long_tracker
    using CSV

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

function load_props_from_csv(path; eval_cols=[:mask, :psi])
    df = DataFrame(CSV.File(path))
    for column in eval_cols
        df[!, column] = eval.(Meta.parse.(df[:, column]))
    end
    return df
end

function check_tracker(
    path;
    ct=(
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
    ),
    thresholds=(
        goodness=(small_floe_area=0.18, large_floe_area=0.236, corr=0.68),
        comp=(mxrot=10, sz=16),
    ),
)
    props = [
        load_props_from_csv(p) for p in readdir(path; join=true) if endswith(p, ".csv")
    ]
    trajectories_ = long_tracker(props, ct, thresholds)
    @show trajectories_

    counts = combine(groupby(trajectories_, :head_uuid), nrow)

    counts[!, :fine] .= counts.nrow .<= length(props)
    @show counts
    @test all(counts.fine)
    return println("\n\n\n-----\n\n\n")
end

begin # Load the data / set config
    @info "obs 1–2"
    check_tracker("notebooks/ellipses/example-31-25-obs1-2")
    # @info "obs 1–3"
    # check_tracker("notebooks/ellipses/example-31-25-obs1-3")
    # @info "obs 5–7"
    # check_tracker("notebooks/ellipses/example-31-25-obs5-7")
end
