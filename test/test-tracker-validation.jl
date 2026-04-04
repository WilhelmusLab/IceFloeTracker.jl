@testsnippet TrackerValidation begin
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

    function tracker_runs_without_error(dataset, case_number, tracker)
        labeled_imgs, pass_times = load_tracker_pair(dataset, case_number)
        result = tracker(labeled_imgs, pass_times)
        return is_wellformed_tracker_result(result)
    end

    function tracker_runs_without_error(case_number; dataset, tracker)
        return tracker_runs_without_error(dataset, case_number, tracker)
    end
end

@testitem "FloeTracker – validated pairs" setup = [TrackerValidation] tags = [:e2e] begin
    import DataFrames: DataFrame
    import IceFloeTracker:
        FloeTracker, FilterFunction, MinimumWeightMatchingFunction, Watkins2026Dataset

    # Each pair is a set of two satellite passes (aqua and terra) for the same
    # location and date in the Watkins 2026 validation dataset.
    # The pairs listed here are those case_numbers for which both passes have
    # validated labeled floe data (fl_analyst != "").
    dataset = Watkins2026Dataset(; ref="v0.1")
    tracker = FloeTracker(;
        filter_function=FilterFunction(), matching_function=MinimumWeightMatchingFunction()
    )

    @test tracker_runs_without_error(dataset, 1, tracker)
    @test tracker_runs_without_error(dataset, 4, tracker)
    @test tracker_runs_without_error(dataset, 5, tracker)
    @test tracker_runs_without_error(dataset, 6, tracker)
    @test tracker_runs_without_error(dataset, 7, tracker)
    @test tracker_runs_without_error(dataset, 8, tracker)
    @test tracker_runs_without_error(dataset, 9, tracker)
    @test tracker_runs_without_error(dataset, 10, tracker)
    @test tracker_runs_without_error(dataset, 11, tracker)
    @test tracker_runs_without_error(dataset, 12, tracker)
    @test tracker_runs_without_error(dataset, 13, tracker)
    @test tracker_runs_without_error(dataset, 14, tracker)
    @test tracker_runs_without_error(dataset, 15, tracker)
    @test tracker_runs_without_error(dataset, 16, tracker)
    @test tracker_runs_without_error(dataset, 17, tracker)
    @test tracker_runs_without_error(dataset, 18, tracker)
    @test tracker_runs_without_error(dataset, 19, tracker)
    @test tracker_runs_without_error(dataset, 21, tracker)
    @test tracker_runs_without_error(dataset, 22, tracker)
    @test tracker_runs_without_error(dataset, 23, tracker)
    @test tracker_runs_without_error(dataset, 25, tracker)
    @test tracker_runs_without_error(dataset, 29, tracker)
    @test tracker_runs_without_error(dataset, 32, tracker)
    @test tracker_runs_without_error(dataset, 33, tracker)
    @test tracker_runs_without_error(dataset, 34, tracker)
    @test tracker_runs_without_error(dataset, 36, tracker)
    @test tracker_runs_without_error(dataset, 37, tracker)
    @test tracker_runs_without_error(dataset, 39, tracker)
    @test tracker_runs_without_error(dataset, 43, tracker)
    @test tracker_runs_without_error(dataset, 44, tracker)
    @test tracker_runs_without_error(dataset, 46, tracker)
    @test tracker_runs_without_error(dataset, 47, tracker)
    @test tracker_runs_without_error(dataset, 48, tracker)
    @test tracker_runs_without_error(dataset, 49, tracker)
    @test tracker_runs_without_error(dataset, 50, tracker)
    @test tracker_runs_without_error(dataset, 51, tracker)
    @test tracker_runs_without_error(dataset, 52, tracker)
    @test tracker_runs_without_error(dataset, 53, tracker) broken = true
    @test tracker_runs_without_error(dataset, 54, tracker)
    @test tracker_runs_without_error(dataset, 55, tracker)
    @test tracker_runs_without_error(dataset, 56, tracker)
    @test tracker_runs_without_error(dataset, 57, tracker)
    @test tracker_runs_without_error(dataset, 58, tracker)
    @test tracker_runs_without_error(dataset, 61, tracker)
    @test tracker_runs_without_error(dataset, 62, tracker)
    @test tracker_runs_without_error(dataset, 63, tracker)
    @test tracker_runs_without_error(dataset, 65, tracker)
    @test tracker_runs_without_error(dataset, 67, tracker)
    @test tracker_runs_without_error(dataset, 68, tracker)
    @test tracker_runs_without_error(dataset, 69, tracker)
    @test tracker_runs_without_error(dataset, 70, tracker)
    @test tracker_runs_without_error(dataset, 71, tracker)
    @test tracker_runs_without_error(dataset, 73, tracker)
    @test tracker_runs_without_error(dataset, 75, tracker)
    @test tracker_runs_without_error(dataset, 77, tracker)
    @test tracker_runs_without_error(dataset, 80, tracker)
    @test tracker_runs_without_error(dataset, 81, tracker)
    @test tracker_runs_without_error(dataset, 83, tracker)
    @test tracker_runs_without_error(dataset, 84, tracker) broken = true
    @test tracker_runs_without_error(dataset, 86, tracker)
    @test tracker_runs_without_error(dataset, 87, tracker)
    @test tracker_runs_without_error(dataset, 93, tracker)
    @test tracker_runs_without_error(dataset, 95, tracker)
    @test tracker_runs_without_error(dataset, 97, tracker)
    @test tracker_runs_without_error(dataset, 98, tracker)
    @test tracker_runs_without_error(dataset, 99, tracker)
    @test tracker_runs_without_error(dataset, 100, tracker)
    @test tracker_runs_without_error(dataset, 103, tracker)
    @test tracker_runs_without_error(dataset, 104, tracker)
    @test tracker_runs_without_error(dataset, 105, tracker) broken = true
    @test tracker_runs_without_error(dataset, 106, tracker)
    @test tracker_runs_without_error(dataset, 107, tracker)
    @test tracker_runs_without_error(dataset, 108, tracker)
    @test tracker_runs_without_error(dataset, 109, tracker)
    @test tracker_runs_without_error(dataset, 110, tracker)
    @test tracker_runs_without_error(dataset, 111, tracker)
    @test tracker_runs_without_error(dataset, 112, tracker)
    @test tracker_runs_without_error(dataset, 115, tracker)
    @test tracker_runs_without_error(dataset, 116, tracker)
    @test tracker_runs_without_error(dataset, 117, tracker)
    @test tracker_runs_without_error(dataset, 118, tracker)
    @test tracker_runs_without_error(dataset, 119, tracker)
    @test tracker_runs_without_error(dataset, 121, tracker)
    @test tracker_runs_without_error(dataset, 122, tracker)
    @test tracker_runs_without_error(dataset, 128, tracker)
    @test tracker_runs_without_error(dataset, 129, tracker)
    @test tracker_runs_without_error(dataset, 130, tracker)
    @test tracker_runs_without_error(dataset, 132, tracker)
    @test tracker_runs_without_error(dataset, 133, tracker)
    @test tracker_runs_without_error(dataset, 135, tracker)
    @test tracker_runs_without_error(dataset, 136, tracker)
    @test tracker_runs_without_error(dataset, 138, tracker)
    @test tracker_runs_without_error(dataset, 140, tracker)
    @test tracker_runs_without_error(dataset, 141, tracker) broken = true
    @test tracker_runs_without_error(dataset, 142, tracker) broken = true
    @test tracker_runs_without_error(dataset, 144, tracker)
    @test tracker_runs_without_error(dataset, 146, tracker)
    @test tracker_runs_without_error(dataset, 148, tracker)
    @test tracker_runs_without_error(dataset, 150, tracker)
    @test tracker_runs_without_error(dataset, 152, tracker)
    @test tracker_runs_without_error(dataset, 155, tracker)
    @test tracker_runs_without_error(dataset, 156, tracker)
    @test tracker_runs_without_error(dataset, 157, tracker)
    @test tracker_runs_without_error(dataset, 158, tracker)
    @test tracker_runs_without_error(dataset, 160, tracker)
    @test tracker_runs_without_error(dataset, 161, tracker)
    @test tracker_runs_without_error(dataset, 164, tracker)
    @test tracker_runs_without_error(dataset, 166, tracker)
    @test tracker_runs_without_error(dataset, 168, tracker)
    @test tracker_runs_without_error(dataset, 171, tracker)
    @test tracker_runs_without_error(dataset, 175, tracker)
    @test tracker_runs_without_error(dataset, 186, tracker)
    @test tracker_runs_without_error(dataset, 188, tracker) broken = true
    @test tracker_runs_without_error(dataset, 189, tracker)
end

@testitem "FloeTracker – sample of cases" setup = [TrackerValidation] tags = [:e2e] begin
    import IceFloeTracker:
        FloeTracker, FilterFunction, MinimumWeightMatchingFunction, Watkins2026Dataset

    # Iterates through a sample of common pairs from the Watkins 2026 validation dataset,
    # loading the validated labeled floes and pass times, and running the default tracker.
    # Selects every 7th case number from among the set of common pairs (case_numbers for
    # which both satellite passes have validated labeled floe data).
    dataset = Watkins2026Dataset(; ref="v0.1")
    tracker = FloeTracker(;
        filter_function=FilterFunction(), matching_function=MinimumWeightMatchingFunction()
    )

    for case_number in
        [7, 14, 21, 49, 56, 63, 70, 77, 98, 112, 119, 133, 140, 161, 168, 175, 189]
        @test tracker_runs_without_error(dataset, case_number, tracker)
    end
    for case_number in [53, 84, 105, 141, 142, 188]
        @test tracker_runs_without_error(dataset, case_number, tracker) broken = true
    end
end
