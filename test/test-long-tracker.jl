
@testitem "Utilities" begin
    using IceFloeTracker:
        long_tracker, _imhist, condition_thresholds, mc_thresholds, get_trajectory_heads
    using CSV

    include("config.jl")

    @testset "get_trajectory_heads" begin
        @testset "basic case" begin
            df = DataFrame([
                (; floe_id=12, group_id=11, passtime=1),
                (; floe_id=27, group_id=16, passtime=5),  # this is the newest entry of the head_uuid=16 trajectory
                (; floe_id=13, group_id=11, passtime=2),
                (; floe_id=14, group_id=11, passtime=3),  # this is the newest entry of the head_uuid=11 trajectory
                (; floe_id=17, group_id=16, passtime=2),
                (; floe_id=11, group_id=11, passtime=0),
                (; floe_id=15, group_id=16, passtime=1),
                (; floe_id=16, group_id=16, passtime=0),
            ])
            # Check that we only get two heads
            heads = get_trajectory_heads(df; group_col=:group_id)
            @test nrow(heads) == 2

            # Check that the heads we get are the ones we want, 
            # despite the fact that the dataframe is unsorted
            sorted_heads = sort(heads, :group_id)
            @test sorted_heads[1, :] == (; group_id=11, floe_id=14, passtime=3)
            @test sorted_heads[2, :] == (; group_id=16, floe_id=27, passtime=5)
        end
        @testset "no existing trajectories" begin
            df = DataFrame([
                (; floe_id=12, group_id=12, passtime=1),
                (; floe_id=27, group_id=27, passtime=1),
                (; floe_id=13, group_id=13, passtime=1),
                (; floe_id=14, group_id=14, passtime=1),
                (; floe_id=17, group_id=17, passtime=1),
                (; floe_id=11, group_id=11, passtime=1),
                (; floe_id=15, group_id=15, passtime=1),
                (; floe_id=16, group_id=16, passtime=1),
            ])
            # Check that we get a head for every row
            heads = get_trajectory_heads(df; group_col=:group_id)
            @test nrow(heads) == 8

            # Check that each head appears once
            @test length(Set(heads.group_id)) == 8
        end
        @testset "wider range of numbers for ranking" begin
            df = DataFrame([
                (; id=13, rank=1),
                (; id=13, rank=2),
                (; id=24, rank=300),
                (; id=24, rank=4),
                (; id=24, rank=0),
                (; id=32, rank=1),
                (; id=32, rank=-1),
                (; id=32, rank=-100),
            ])
            # Check that we only get three heads
            heads = get_trajectory_heads(df; group_col=:id, order_col=:rank)
            @test nrow(heads) == 3

            # Check that each head appears once
            @test length(Set(heads.id)) == 3

            # Check that the heads are the ones we care about
            sorted_heads = sort(heads, :id)
            @test sorted_heads[1, :] == (; id=13, rank=2)
            @test sorted_heads[2, :] == (; id=24, rank=300)
            @test sorted_heads[3, :] == (; id=32, rank=1)
        end
    end
end

@testitem "Basic cases" begin
    using Random
    using DataFrames
    using IceFloeTracker: _imhist, condition_thresholds, mc_thresholds

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

    @testset "Case 1" begin
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

    @testset "Case 2" begin
        trajectories = IceFloeTracker.long_tracker(
            props_test_case2, condition_thresholds, mc_thresholds
        )

        # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2
        IDs = trajectories[!, :ID]
        @test IDs == [1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5, 5]
    end

    @testset "Test gaps" begin
        @testset "Case 3" begin
            # Every floe is matched in every day for which there is data
            Random.seed!(123)
            props = addgaps(_props)

            trajectories = IceFloeTracker.long_tracker(
                props, condition_thresholds, mc_thresholds
            )

            # Expected: 5 trajectories, all of which have length 3 as in test case 1
            IDs = trajectories[!, :ID]
            @test IDs == [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5]
        end

        @testset "Case 4" begin
            # Add gaps to props_test_case2
            Random.seed!(123)
            props = addgaps(props_test_case2)
            trajectories = IceFloeTracker.long_tracker(
                props, condition_thresholds, mc_thresholds
            )

            # Expected: 5 trajectories, 3 of which have length 3 and 2 of which have length 2 as in test case 2
            IDs = trajectories[!, :ID]
            @test IDs == [1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 5, 5, 5]
        end
    end
end

@testitem "Ellipses" begin
    using CSV
    using DataFrames
    using IceFloeTracker: long_tracker, condition_thresholds, mc_thresholds

    function load_props_from_csv(path; eval_cols=[:mask, :psi])
        df = DataFrame(CSV.File(path))
        for column in eval_cols
            df[!, column] = eval.(Meta.parse.(df[:, column]))
        end
        return df
    end

    function check_tracker_results(path)
        props = [
            load_props_from_csv(p) for p in readdir(path; join=true) if endswith(p, ".csv")
        ]
        trajectories_ = long_tracker(props, condition_thresholds, mc_thresholds)

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

        modified_condition_thresholds = (
            search_thresholds=condition_thresholds.search_thresholds,
            small_floe_settings=(;
                condition_thresholds.small_floe_settings..., minimumarea=1200
            ),
            large_floe_settings=condition_thresholds.large_floe_settings,
        )
        trajectories_ = long_tracker(props, modified_condition_thresholds, mc_thresholds)

        @test all(
            modified_condition_thresholds.small_floe_settings.minimumarea .<=
            trajectories_.area,
        )
    end
end
