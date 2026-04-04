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
