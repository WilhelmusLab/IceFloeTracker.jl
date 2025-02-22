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
    props = vcat(props[1:(end - 1)], blank_props, [props[end]])
    return props
end

begin # Set thresholds
    search_thresholds = (dt=(30.0, 100.0, 1300.0), dist=(200, 250, 300))
    large_floe_settings = (
        area=1200,
        arearatio=0.28,
        majaxisratio=0.10,
        minaxisratio=0.12,
        convexarearatio=0.14,
    )
    small_floe_settings = (
        area=1200, arearatio=0.18, majaxisratio=0.1, minaxisratio=0.15, convexarearatio=0.2
    )
    condition_thresholds = (search_thresholds, small_floe_settings, large_floe_settings)
    mc_thresholds = (
        goodness=(small_floe_area=0.18, large_floe_area=0.236, corr=0.68),
        comp=(mxrot=10, sz=16),
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
    Random.seed!(123)
    IceFloeTracker.adduuid!(_props)
end

begin # Filter out floes with area less than `floe_area_threshold` pixels
    floe_area_threshold = 400
    for (i, prop) in enumerate(_props)
        _props[i] = prop[prop[:, :area] .>= floe_area_threshold, :] # 500 working good
        sort!(_props[i], :area; rev=true)
    end
end

@ntestset "$(@__FILE__)" begin
    @ntestset "Case 1" begin
        # Every floe is matched in every day
        props_test_case1 = deepcopy(_props)
        trajectories = IceFloeTracker.long_tracker(
            props_test_case1, condition_thresholds, mc_thresholds
        )

        # Expected: 5 trajectories, all of which have length 3
        IDs = trajectories[!, :ID]
        ids, counts = _imhist(IDs, unique(IDs))
        @test maximum(ids) == 5

        ids, counts = _imhist(counts, unique(counts))
        @test ids == [3]
        @test counts == [5]
    end

    begin # Unmatched floe in day 1, unmatched floe in day 2, and matches for every floe starting in day 3
        props_test_case2 = deepcopy(_props)
        deleteat!(props_test_case2[1], 1)
        deleteat!(props_test_case2[2], 5)
    end

    @ntestset "Case 2" begin
        trajectories = IceFloeTracker.long_tracker(
            props_test_case2, condition_thresholds, mc_thresholds
        )

        # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2
        IDs = trajectories[!, :ID]
        @test IDs == [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5]
    end

    @ntestset "Test gaps" begin
        @ntestset "Case 3" begin
            # Every floe is matched in every day for which there is data

            props = addgaps(_props)

            trajectories = IceFloeTracker.long_tracker(
                props, condition_thresholds, mc_thresholds
            )

            # Expected: 5 trajectories, all of which have length 3 as in test case 1
            IDs = trajectories[!, :ID]
            @test IDs == [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5]
        end

        @ntestset "Case 4" begin
            # Add gaps to props_test_case2
            props = addgaps(props_test_case2)
            trajectories = IceFloeTracker.long_tracker(
                props, condition_thresholds, mc_thresholds
            )

            # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2 as in test case 2
            IDs = trajectories[!, :ID]
            @test IDs == [1, 1, 1, 2, 2, 2, 3, 3, 4, 4, 4, 5, 5]
        end
    end

    @ntestset "Test Case 5: drop_trajectories_length1" begin
        df = DataFrame(; ID=[1, 1, 2, 3, 3, 3, 4], foo=1:7)
        df = IceFloeTracker.drop_trajectories_length1(df, :ID)
        @test df == DataFrame(; ID=[1, 1, 3, 3, 3], foo=[1, 2, 4, 5, 6])

        # Test empty dataframe
        df = similar(df, 0)
        df_ = IceFloeTracker.drop_trajectories_length1(df, :ID)
        @test df_ == df

        # Test no trajectories of length 1
        df = DataFrame(; ID=[1, 1, 2, 2, 3, 3], foo=1:6)
        df_ = IceFloeTracker.drop_trajectories_length1(df, :ID)
        @test df_ == df
    end
end
