@testitem "FloeTracker – utilities" begin
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
            heads = _get_trajectory_heads(
                df, current_time_step, maximum_time_step; group_col=:group_id
            )
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
            heads = _get_trajectory_heads(
                df, current_time_step, maximum_time_step; group_col=:group_id
            )
            @test nrow(heads) == 8

            # Check that each head appears once
            @test length(Set(heads.group_id)) == 8
        end
    end
end
# TODO: Include tests for other floe_tracker utilities

@testitem "FloeTracker – basic cases" begin
    using Random
    using DataFrames
    using IceFloeTracker: floe_tracker, FilterFunction, MinimumWeightMatchingFunction
    using Serialization: deserialize
    import Images: label_components
    using Dates

    begin # Load data
        pth = joinpath("test_inputs", "tracker")
        _floedata = deserialize(joinpath(pth, "tracker_test_data.dat"))
        _passtimes = deserialize(joinpath(pth, "passtimes.dat"))
        # The serialized data has props tables, but we'll generate those ourselves
        _props, _imgs = deepcopy.([_floedata.props, _floedata.imgs])
        labeled_imgs = label_components.(_imgs)
    end

    floe_area_threshold = 400

    @testset "Case 1" begin
        # Every floe is matched in every day
        tracker = FloeTracker(;
            filter_function=FilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
            minimum_area=floe_area_threshold,
        )
        trajectories = tracker(labeled_imgs, _passtimes)

        # Expected: 5 trajectories, all of which have length 3
        # (other floes are below the area threshold)
        counts = combine(groupby(trajectories, [:ID]), nrow => :count)
        @test nrow(counts) == 5
        @test all(counts[!, :count] .== 3)
    end

    @testset "Case 2" begin
        # Add single floe gaps
        labeled_imgs_gaps = deepcopy(labeled_imgs)
        labeled_imgs_gaps[2][labeled_imgs_gaps[2] .== 36] .= 0
        labeled_imgs_gaps[3][labeled_imgs_gaps[3] .== 33] .= 0

        tracker = FloeTracker(;
            filter_function=FilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
            minimum_area=floe_area_threshold,
        )

        trajectories = tracker(labeled_imgs_gaps, _passtimes)

        # Expected: 5 trajectories, 4 of which have length 3 and 1 of which have length 2
        counts = combine(groupby(trajectories, [:ID]), nrow => :count)
        @test sum(counts[:, :count] .== 3) == 4 && sum(counts[:, :count] .== 2) == 1
    end

    @testset "Test gaps" begin
        @testset "Case 3" begin
            # Every floe is matched in every day for which there is data
            # Here we insert a blank image into the series
            labeled_imgs_gaps = [
                labeled_imgs[1], labeled_imgs[2], labeled_imgs[2] * 0, labeled_imgs[3]
            ]
            tracker = FloeTracker(;
                filter_function=FilterFunction(),
                matching_function=MinimumWeightMatchingFunction(),
                minimum_area=floe_area_threshold,
            )
            # Add an extra pass-time to simulate a longer time series
            passtimes_gaps = [
                _passtimes[1], _passtimes[2], _passtimes[3], DateTime("2022-09-16T12:44:49")
            ]

            trajectories = tracker(labeled_imgs_gaps, passtimes_gaps)

            # Expected: 5 trajectories, all of which have length 3 as in test case 1
            IDs = trajectories[!, :ID]
            counts = combine(groupby(trajectories, [:ID]), nrow => :count)
            @test nrow(counts) == 5
            @test all(counts[!, :count] .== 3)
        end

        @testset "Case 4" begin
            tracker = FloeTracker(;
                filter_function=FilterFunction(),
                matching_function=MinimumWeightMatchingFunction(),
                minimum_area=floe_area_threshold,
            )

            # Add full image gap
            labeled_imgs_gaps = [
                labeled_imgs[1], labeled_imgs[2], labeled_imgs[2] * 0, labeled_imgs[3]
            ]

            # Add single floe gaps
            labeled_imgs_gaps[2][labeled_imgs_gaps[2] .== 36] .= 0
            labeled_imgs_gaps[4][labeled_imgs_gaps[4] .== 33] .= 0

            # Extend passtimes
            passtimes_gaps = [
                _passtimes[1], _passtimes[2], _passtimes[3], DateTime("2022-09-16T12:44:49")
            ]

            tracker = FloeTracker(;
                filter_function=FilterFunction(),
                matching_function=MinimumWeightMatchingFunction(),
                minimum_area=floe_area_threshold,
            )

            trajectories = tracker(labeled_imgs_gaps, passtimes_gaps)
            counts = combine(groupby(trajectories, [:ID]), nrow => :count)
            @test sum(counts[:, :count] .== 3) == 4 && sum(counts[:, :count] .== 2) == 1
        end
    end
end

@testitem "FloeTracker – ellipses" begin
    using CSVFiles
    using DataFrames
    using IceFloeTracker: floe_tracker, FilterFunction, MinimumWeightMatchingFunction

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

        # This test uses the inner floe tracker function, which takes a props argument instead of 
        # the list of images
        trajectories_ = floe_tracker(
            props, FilterFunction(), MinimumWeightMatchingFunction()
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
            props, FilterFunction(), MinimumWeightMatchingFunction(); minimum_area=1200
        )

        @test all(1200 .<= trajectories_.area)
    end
end

@testitem "FloeTracker – pipeline" begin
    import Dates: DateTime
    import DataFrames: nrow, DataFrame
    import IceFloeTracker.Tracking:
        FilterFunction,
        MinimumWeightMatchingFunction,
        ChainedFilterFunction,
        DistanceThresholdFilter,
        RelativeErrorThresholdFilter

    dataset = Watkins2026Dataset(; ref="v0.1")

    @testset "Basic functionality" begin
        filter!(c -> c.case_number == 6, dataset)
        sort!([:pass_time], dataset)
        segmenter = LopezAcosta2019Tiling.Segment()
        segmentation_results =
            segmenter.(
                modis_truecolor.(dataset),
                modis_falsecolor.(dataset),
                modis_landmask.(dataset),
            )
        tracker = FloeTracker(;
            filter_function=FilterFunction(),
            matching_function=MinimumWeightMatchingFunction(),
        )
        tracking_results = tracker(segmentation_results, info(dataset).pass_time)
        @test isa(tracking_results, DataFrame)
        @test "trajectory_uuid" in names(tracking_results)
        @test nrow(tracking_results) == 116 # Just a test of consistency, not indicator of correctness
    end

    @testset "Simpler filter function" begin

        # Note: Only tests that the function runs, does not check the error in the results!
        filter!(c -> c.case_number == 6, dataset)
        sort!([:pass_time], dataset)
        segmenter = LopezAcosta2019Tiling.Segment()
        segmentation_results =
            segmenter.(
                modis_truecolor.(dataset),
                modis_falsecolor.(dataset),
                modis_landmask.(dataset),
            )

        tracker = FloeTracker(;
            filter_function=ChainedFilterFunction(;
                filters=[
                    DistanceThresholdFilter(),
                    RelativeErrorThresholdFilter(; variable=:area),
                ],
            ),
            matching_function=MinimumWeightMatchingFunction(;
                columns=[:scaled_distance, :relative_error_area],
                weights=ones(2), # Not yet used
            ),
        )
        tracking_results = tracker(segmentation_results, info(dataset).pass_time)
        @test isa(tracking_results, DataFrame)
        @test "trajectory_uuid" in names(tracking_results)
        @test nrow(tracking_results) == 150 # Just checks consistency -- not a direct indicator of quality!
    end

    # TODO: LopezAcosta filter function
    # TODO: LopezAcosta matching function
    # TODO: Update to use error metric for matching
end

@testsnippet TrackerValidation begin
    using DataFrames: DataFrame, nrow
    using Dates: DateTime
    function tracker_runs_without_error(
        img1::Matrix{Int},
        time1::DateTime,
        img2::Matrix{Int},
        time2::DateTime;
        tracker::AbstractTracker,
    )
        result = tracker([img1, img2], [time1, time2])
        return is_wellformed_tracker_result(result)
    end

    """
        is_wellformed_tracker_result(result) -> Bool

    Check that the result of running the FloeTracker is a well-formed DataFrame.
    A well-formed result is a DataFrame that, when non-empty, contains the expected
    trajectory and floe property columns.
    """
    function is_wellformed_tracker_result(result)
        result isa DataFrame || return false
        if nrow(result) > 0
            expected_columns = ["ID", "trajectory_uuid", "uuid", "passtime"]
            all(col in names(result) for col in expected_columns) || return false
        end
        return true
    end
end

@testitem "FloeTracker – validated cases" setup = [TrackerValidation] begin
    import DataFrames: DataFrame, nrow

    """
        load_tracker_pair(dataset, case_number) -> (labeled_imgs, pass_times)

    Load the validated labeled floes and pass times for a pair of satellite observations
    with the given `case_number` from the validation `dataset`.

    Filters to observations that have validated floe labels (`fl_analyst != ""`),
    then sorts by pass time so the earliest observation comes first.

    Returns a tuple `(labeled_imgs, pass_times)` where:
    - `labeled_imgs`: vector of SegmentedImages with validated floe labels
    - `pass_times`: vector of DateTime observation times
    """
    function load_tracker_pair(dataset, case_number)
        pair = filter(c -> c.case_number == case_number && c.fl_analyst != "", dataset)
        sort!([:pass_time], pair)
        labeled_imgs = validated_labeled_floes.(pair)
        pass_times = info(pair).pass_time
        return labeled_imgs, pass_times
    end

    function tracker_runs_without_error(dataset, case_number, tracker)
        labeled_imgs, pass_times = load_tracker_pair(dataset, case_number)
        result = tracker(labeled_imgs, pass_times)
        return is_wellformed_tracker_result(result)
    end

    # Iterates through a sample of common pairs from the Watkins 2026 validation dataset,
    # loading the validated labeled floes and pass times, and running the default tracker.
    dataset = Watkins2026Dataset(; ref="v0.1")
    tracker = FloeTracker(;
        filter_function=FilterFunction(), matching_function=MinimumWeightMatchingFunction()
    )

    all_cases_with_validated_floes = case_number.(filter(c -> c.fl_analyst != "", dataset))
    known_broken_cases = [
        53, # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
        84, # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
        105, # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
        141, # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
        142, # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
        188, # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    ]
    working_cases = setdiff(all_cases_with_validated_floes, known_broken_cases)

    # Some cases are known to be broken due to issues in the tracker.
    # Each of these cases should be fixed before being removed from this list.
    @testset "Known broken cases" begin
        for case_number in known_broken_cases
            @test tracker_runs_without_error(dataset, case_number, tracker) broken = true
        end
    end

    # The rest of the cases should run without error.
    @testset "Working cases" begin
        for case_number in working_cases
            @test tracker_runs_without_error(dataset, case_number, tracker)
        end
    end
end

@testitem "FloeTracker – synthetic shapes" setup = [TrackerValidation] begin
    using Dates: DateTime
    function tracker_runs_without_error(
        tracker::AbstractTracker,
        img::Matrix{Int};
        start=DateTime("2025-01-01T00:00:00"),
        end_=DateTime("2025-01-01T00:00:01"),
    )
        return is_wellformed_tracker_result(tracker([img, img], [start, end_]))
    end

    tracker = FloeTracker(;
        filter_function=FilterFunction(),
        matching_function=MinimumWeightMatchingFunction(),
        minimum_area=1,
    )

    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0
            0 0 0
            0 0 0
        ],
    )

    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0
            0 1 0
            0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/911

    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0
            0 1 0
            0 1 0
            0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 0 0
            0 0 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 0 1 0
            0 1 0 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/912
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 1 0
            0 0 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 0 0
            0 1 1 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 0 1 0 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/913
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 1 0
            0 1 1 0
            0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0
            0 1 0 0
            0 1 1 0
            0 1 0 0
            0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 1 1 1 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 0 0
            0 1 1 0 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 0 1 0 0
            0 1 1 1 1 0
            0 0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0 0
            0 1 1 1 1 1 0
            0 0 0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 0 0
            0 1 0 0 0
            0 0 0 0 0
        ],
    ) broken = true # https://github.com/WilhelmusLab/IceFloeTracker.jl/issues/919
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 1 1 1 0
            0 1 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 0
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 0 0
            0 0 1 0 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 0 0
            0 1 1 1 0
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 0 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 0 1 0 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 0 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 0 1 1 1 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 0 1 1 0
            0 1 1 1 1
            0 0 1 1 0
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0
            0 1 1 1 0
            0 1 1 1 0
            0 0 1 1 1
            0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 1 1 1 0 0
            0 0 1 0 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 0 0
            0 1 1 1 0 0
            0 1 1 1 1 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 1 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 1 1 1 1 0
            0 1 1 1 1 0
            0 0 1 1 0 0
            0 0 0 0 0 0
        ],
    )
    @test tracker_runs_without_error(
        tracker,
        Int[
            0 0 0 0 0 0
            0 0 1 1 1 0
            0 1 1 1 1 0
            0 1 1 1 0 0
            0 0 0 0 0 0
        ],
    )
end