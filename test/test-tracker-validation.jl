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
        filter_function=FilterFunction(),
        matching_function=MinimumWeightMatchingFunction(),
    )

    @testset "Pair case 1" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 1)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 4" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 4)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 5" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 5)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 6" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 6)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 7" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 7)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 8" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 8)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 9" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 9)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 10" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 10)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 11" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 11)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 12" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 12)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 13" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 13)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 14" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 14)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 15" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 15)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 16" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 16)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 17" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 17)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 18" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 18)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 19" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 19)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 21" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 21)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 22" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 22)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 23" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 23)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 25" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 25)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 29" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 29)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 32" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 32)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 33" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 33)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 34" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 34)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 36" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 36)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 37" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 37)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 39" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 39)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 43" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 43)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 44" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 44)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 46" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 46)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 47" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 47)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 48" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 48)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 49" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 49)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 50" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 50)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 51" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 51)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 52" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 52)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 53" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 53)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 54" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 54)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 55" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 55)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 56" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 56)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 57" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 57)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 58" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 58)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 61" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 61)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 62" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 62)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 63" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 63)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 65" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 65)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 67" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 67)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 68" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 68)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 69" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 69)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 70" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 70)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 71" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 71)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 73" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 73)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 75" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 75)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 77" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 77)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 80" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 80)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 81" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 81)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 83" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 83)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 84" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 84)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 86" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 86)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 87" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 87)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 93" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 93)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 95" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 95)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 97" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 97)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 98" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 98)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 99" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 99)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 100" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 100)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 103" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 103)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 104" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 104)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 105" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 105)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 106" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 106)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 107" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 107)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 108" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 108)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 109" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 109)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 110" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 110)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 111" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 111)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 112" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 112)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 115" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 115)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 116" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 116)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 117" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 117)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 118" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 118)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 119" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 119)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 121" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 121)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 122" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 122)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 128" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 128)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 129" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 129)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 130" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 130)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 132" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 132)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 133" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 133)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 135" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 135)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 136" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 136)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 138" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 138)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 140" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 140)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 141" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 141)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 142" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 142)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 144" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 144)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 146" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 146)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 148" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 148)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 150" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 150)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 152" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 152)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 155" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 155)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 156" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 156)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 157" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 157)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 158" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 158)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 160" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 160)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 161" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 161)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 164" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 164)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 166" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 166)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 168" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 168)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 171" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 171)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 175" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 175)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 186" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 186)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 188" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 188)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
    @testset "Pair case 189" begin
        labeled_imgs, pass_times = load_tracker_pair(dataset, 189)
        result = tracker(labeled_imgs, pass_times)
        @test is_wellformed_tracker_result(result)
    end
end

@testitem "FloeTracker – sample of cases" setup = [TrackerValidation] tags = [:e2e] begin
    import DataFrames: DataFrame
    import IceFloeTracker:
        FloeTracker, FilterFunction, MinimumWeightMatchingFunction, Watkins2026Dataset

    # Iterates through a sample of common pairs from the Watkins 2026 validation dataset,
    # loading the validated labeled floes and pass times, and running the default tracker.
    # Selects every 7th case number from among the set of common pairs (case_numbers for
    # which both satellite passes have validated labeled floe data).
    dataset = Watkins2026Dataset(; ref="v0.1")
    tracker = FloeTracker(;
        filter_function=FilterFunction(),
        matching_function=MinimumWeightMatchingFunction(),
    )

    sample_case_numbers = [7, 14, 21, 49, 56, 63, 70, 77, 84, 98, 105, 112, 119, 133, 140, 161, 168, 175, 189]

    successes = Bool[]
    for case_number in sample_case_numbers
        labeled_imgs, pass_times = load_tracker_pair(dataset, case_number)
        result = tracker(labeled_imgs, pass_times)
        push!(successes, is_wellformed_tracker_result(result))
    end
    @test all(successes)
end
