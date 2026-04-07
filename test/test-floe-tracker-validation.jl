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

@testitem "FloeTracker – sample of cases" setup = [TrackerValidation] tags = [:e2e] begin
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
    known_broken_cases = [53, 84, 105, 141, 142, 188]
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
