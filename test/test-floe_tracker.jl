
@testitem "Utilities" begin
    using DataFrames
    import IceFloeTracker.Tracking: _get_trajectory_heads
    using Dates

    include("config.jl")

    @testset "get_trajectory_heads" begin
        @testset "basic case" begin
            
            test_time0 = DateTime("2020-01-01 0:00", "y-m-d H:M")
            test_time1 = DateTime("2020-01-01 1:00", "y-m-d H:M")
            test_time2 = DateTime("2020-01-01 2:00", "y-m-d H:M")
            test_time3 = DateTime("2020-01-01 3:00", "y-m-d H:M")
            test_time5 = DateTime("2020-01-01 5:00", "y-m-d H:M")
            df = DataFrame([
                (; floe_id=12, group_id=11, passtime=test_time0),
                (; floe_id=27, group_id=16, passtime=test_time5),  # this is the newest entry of the head_uuid=16 trajectory
                (; floe_id=13, group_id=11, passtime=test_time2),
                (; floe_id=14, group_id=11, passtime=test_time3),  # this is the newest entry of the head_uuid=11 trajectory
                (; floe_id=17, group_id=16, passtime=test_time2),
                (; floe_id=11, group_id=11, passtime=test_time0),
                (; floe_id=15, group_id=16, passtime=test_time1),
                (; floe_id=16, group_id=16, passtime=test_time0),
            ])
            # Check that we only get two heads
            current_time_step = test_time5
            maximum_time_step = Day(5)
            heads = _get_trajectory_heads(df, current_time_step, maximum_time_step; group_col=:group_id)
            @test nrow(heads) == 2

            # Check that the heads we get are the ones we want, 
            # despite the fact that the dataframe is unsorted
            sorted_heads = sort(heads, :group_id)
            @test sorted_heads[1, :] == (; group_id=11, floe_id=14, passtime=test_time3)
            @test sorted_heads[2, :] == (; group_id=16, floe_id=27, passtime=test_time5)
        end
        @testset "no existing trajectories" begin
            test_time1 = DateTime("2020-01-01 1:00", "y-m-d H:M")
            df = DataFrame([
                (; floe_id=12, group_id=12, passtime=test_time1),
                (; floe_id=27, group_id=27, passtime=test_time1),
                (; floe_id=13, group_id=13, passtime=test_time1),
                (; floe_id=14, group_id=14, passtime=test_time1),
                (; floe_id=17, group_id=17, passtime=test_time1),
                (; floe_id=11, group_id=11, passtime=test_time1),
                (; floe_id=15, group_id=15, passtime=test_time1),
                (; floe_id=16, group_id=16, passtime=test_time1),
            ])
            # Check that we get a head for every row
            current_time_step = test_time1
            maximum_time_step = Day(1)
            heads = _get_trajectory_heads(df, current_time_step, maximum_time_step; group_col=:group_id)
            @test nrow(heads) == 8

            # Check that each head appears once
            @test length(Set(heads.group_id)) == 8
        end
    end
end
# TODO: Include tests for other floe_tracker utilities


@testitem "Basic cases" begin
    using Random
    using DataFrames
    using IceFloeTracker: floe_tracker, FilterFunction, MinimumWeightMatchingFunction
    using Serialization: deserialize

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

    begin # Load data
        pth = joinpath("test_inputs", "tracker")
        _floedata = deserialize(joinpath(pth, "tracker_test_data.dat"))
        _passtimes = deserialize(joinpath(pth, "passtimes.dat"))
        _props, _imgs = deepcopy.([_floedata.props, _floedata.imgs])

        # This order is important: masks, uuids, passtimes, ψs
        add_floemasks!.(_props, _imgs)
        add_ψs!(_props)
        add_passtimes!(_props, _passtimes)
        Random.seed!(123)
        add_uuids!(_props)
    end

    begin # Filter out floes with area less than `floe_area_threshold` pixels
        floe_area_threshold = 400
        for (i, prop) in enumerate(_props)
            _props[i] = prop[prop[:, :area] .>= floe_area_threshold, :] # 500 working good
            sort!(_props[i], :area; rev=true)
        end
    end

    @testset "Case 1" begin
        # Every floe is matched in every day
        props_test_case1 = deepcopy(_props)
        trajectories = floe_tracker(
            props_test_case1, FilterFunction, MinimumWeightMatchingFunction()
        )

        # Expected: 5 trajectories, all of which have length 3
        counts = combine(groupby(trajectories, [:ID]), nrow => :count)
        @test nrow(counts) == 5
        @test all(counts[!, :count] .== 3)
    end

    begin # Unmatched floe in day 1, unmatched floe in day 2, and matches for every floe starting in day 3
        props_test_case2 = deepcopy(_props)
        deleteat!(props_test_case2[1], 1)
        deleteat!(props_test_case2[2], 5)
    end

    @testset "Case 2" begin
        trajectories = IceFloeTracker.floe_tracker(
            props_test_case2, FilterFunction, IceFloeTracker.MinimumWeightMatchingFunction()
        )

        # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2
        
        counts = combine(groupby(trajectories, [:ID]), nrow => :count)
        @test sum(counts[:, :count] .== 3) == 3 && sum(counts[:, :count] .== 2) == 2
    end

    @testset "Test gaps" begin
        @testset "Case 3" begin
            # Every floe is matched in every day for which there is data
            Random.seed!(123)
            props = addgaps(_props)

            trajectories = IceFloeTracker.floe_tracker(
                props, FilterFunction, MinimumWeightMatchingFunction()
            )

            # Expected: 5 trajectories, all of which have length 3 as in test case 1
            IDs = trajectories[!, :ID]
            counts = combine(groupby(trajectories, [:ID]), nrow => :count)
            @test nrow(counts) == 5
            @test all(counts[!, :count] .== 3)
        end

        @testset "Case 4" begin
            # Add gaps to props_test_case2
            Random.seed!(123)
            props = addgaps(props_test_case2)
            trajectories = IceFloeTracker.floe_tracker(
                props, FilterFunction, IceFloeTracker.MinimumWeightMatchingFunction()
            )

            # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2 as in test case 2
            counts = combine(groupby(trajectories, [:ID]), nrow => :count)
            @test sum(counts[:, :count] .== 3) == 3 && sum(counts[:, :count] .== 2) == 2
        end
    end
end

@testitem "Ellipses" begin
    using CSVFiles
    using DataFrames
    using IceFloeTracker:
        floe_tracker, FilterFunction, MinimumWeightMatchingFunction

    function load_props_from_csv(path; eval_cols=[:mask, :psi])
        df = DataFrame(load(path))
        for column in eval_cols
            df[!, column] = eval.(Meta.parse.(df[:, column]))
        end
        return df
    end

    function check_tracker_results(path)
        props = [
            load_props_from_csv(p) for p in readdir(path; join=true) if endswith(p, ".csv")
        ]
        # TODO: Check types for the ShapeDifference function. What's different about these props tables?
        trajectories_ = floe_tracker(
            props, FilterFunction, MinimumWeightMatchingFunction()
        )

        trajectory_lengths = combine(groupby(trajectories_, :trajectory_uuid), nrow)

        # Each trajectory is at most the legnth of the dataset
        # Weak test for a regression where a trajectory would have more than one element for a particular day
        trajectory_lengths[!, :not_longer_than_dataset] .=
            trajectory_lengths.nrow .<= length(props)

        @test all(trajectory_lengths.not_longer_than_dataset)

        # Each trajectory is at least two rows long – all single-match trajectories are removed.
        trajectory_lengths[!, :longer_than_one] .= trajectory_lengths.nrow .>= 2
        @test all(trajectory_lengths.longer_than_one)

        # Each UUID appears at most once
        # Weak test for a regression where a trajectory would have more than one element for a particular day, 
        # and one floe might be matched multiple times
        uuid_counts = combine(groupby(trajectories_, :uuid), nrow)
        @test all(uuid_counts.nrow .== 1)
    end

    @testset "10 observations of 2 floes" begin
        check_tracker_results(
            joinpath("test_inputs", "tracker", "ellipses", "example-2floes-10obs")
        )
    end
    @testset "10 observations of 40 floes" begin
        check_tracker_results(
            joinpath("test_inputs", "tracker", "ellipses", "example-40floes-10obs")
        )
    end
    @testset "some observations missing" begin
        check_tracker_results(
            joinpath("test_inputs", "tracker", "ellipses", "example-floes-missing-10obs")
        )
    end
    @testset "exclude small floes" begin
        path = joinpath("test_inputs", "tracker", "ellipses", "example-40floes-10obs")
        props = [
            load_props_from_csv(p) for p in readdir(path; join=true) if endswith(p, ".csv")
        ]

        trajectories_ = floe_tracker(
            props, FilterFunction, MinimumWeightMatchingFunction(); minimum_area=1200
        )

        @test all(
            1200 .<= trajectories_.area,
        )
    end
end
