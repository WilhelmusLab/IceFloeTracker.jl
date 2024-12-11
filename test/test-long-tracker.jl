using IceFloeTracker: long_tracker, _imhist

"""
addgaps(props)

Add gaps to the props array after the first and before the last day.
"""
function addgaps(props)
    blank_props = fill(similar(props[1], 0), rand(1:5))

    # add gap after first day
    props = vcat(props[1:1], blank_props, props[2:end])
    # add gap before last day
    props = vcat(props[1:end-1], blank_props, [props[end]])
    return props
end


begin # Set thresholds
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


begin # Load data
    pth = joinpath("test_inputs", "tracker")
    _floedata = deserialize(joinpath(pth, "tracker_test_data.dat"))
    _passtimes = deserialize(joinpath(pth, "passtimes.dat"))
    _props, _imgs = deepcopy.([_floedata.props, _floedata.imgs])

    # This order is important: masks, uuids, passtimes, ψs
    IceFloeTracker.addfloemasks!(_props, _imgs)
    IceFloeTracker.addψs!(_props)
    IceFloeTracker.add_passtimes!(_props, _passtimes)
    IceFloeTracker.adduuid!(_props)
end

begin # Filter out floes with area less than `floe_area_threshold` pixels
    floe_area_threshold = 400
    for (i, prop) in enumerate(_props)
        _props[i] = prop[prop[:, :area].>=floe_area_threshold, :] # 500 working good
        sort!(_props[i], :area, rev=true)
    end
end

@ntestset "$(@__FILE__)" begin
    @ntestset "Case 1" begin
        # Every floe is matched in every day
        props_test_case1 = deepcopy(_props)
        trajectories = IceFloeTracker.long_tracker(props_test_case1, condition_thresholds, mc_thresholds)

        # Expected: 5 trajectories, all of which have length 3
        IDs = trajectories[!, :ID]
        ids, counts = _imhist(IDs, unique(IDs))
        @test maximum(ids) == 5

        ids, counts = _imhist(counts, unique(counts))
        @test ids == [3]
        @test counts == [5]
    end

    @ntestset "Case 2" begin
        # Unmatched floe in day 1, unmatched floe in day 2, and matches for every floe starting in day 3
        props_test_case2 = deepcopy(_props)
        delete!(props_test_case2[1], 1)
        delete!(props_test_case2[2], 5)
        trajectories = IceFloeTracker.long_tracker(props_test_case2, condition_thresholds, mc_thresholds)

        # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2
        IDs = trajectories[!, :ID]
        ids, counts = _imhist(IDs, unique(IDs))
        @test maximum(ids) == 5

        ids, counts = _imhist(counts, unique(counts))
        @test all(ids .== counts)
    end

    @ntestset "Test gaps" begin
        @ntestset "Case 3" begin
            # Every floe is matched in every day for which there is data

            props = addgaps(props_test_case1)

            trajectories = IceFloeTracker.long_tracker(props, condition_thresholds, mc_thresholds)

            # Expected: 5 trajectories, all of which have length 3 as in test case 1
            IDs = trajectories[!, :ID]
            ids, counts = _imhist(IDs, unique(IDs))
            @test maximum(ids) == 5

            ids, counts = _imhist(counts, unique(counts))
            @test ids == [3]
            @test counts == [5]
        end

        @ntestset "Case 4" begin
            # Every floe is matched in every day for which there is data

            props = addgaps(props_test_case2)

            trajectories = IceFloeTracker.long_tracker(props, condition_thresholds, mc_thresholds)

            # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2 as in test case 2
            IDs = trajectories[!, :ID]
            ids, counts = _imhist(IDs, unique(IDs))
            @test maximum(ids) == 5

            ids, counts = _imhist(counts, unique(counts))
            @test all(ids .== counts)
        end
    end

end
