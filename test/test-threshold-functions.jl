@testitem "Time-Distance Thresholds" begin
    using Dates: Minute, Hour, Day
    # Test short and long times

    lopez_acosta_dist = LopezAcostaTimeDistanceFunction()
    loglog_quadratic = LogLogQuadraticTimeDistanceFunction()
    linear = LinearTimeDistanceFunction()

    @test distance_threshold(1, Minute(19), lopez_acosta_dist)
    @test !distance_threshold(7e5, Minute(19), lopez_acosta_dist)
    @test distance_threshold(1, Minute(19), loglog_quadratic)
    @test !distance_threshold(3e3, Minute(19), loglog_quadratic)
    @test distance_threshold(1, Minute(1), linear)
    @test !distance_threshold(1e3, Minute(1), linear)

    @test distance_threshold(1e3, Hour(15), lopez_acosta_dist)
    @test !distance_threshold(1e5, Hour(15), lopez_acosta_dist)
    @test distance_threshold(1e3, Hour(15), loglog_quadratic)
    @test !distance_threshold(1e5, Hour(15), loglog_quadratic)

    # LopezAcosta allows any length of time but limits the maximum distance
    @test distance_threshold(1, Day(10), lopez_acosta_dist)
    @test !distance_threshold(1e5, Day(10), lopez_acosta_dist)

    # LogLogQuadratic has a default maximum time of 7 days
    @test !distance_threshold(1, Day(10), loglog_quadratic)
    @test !distance_threshold(1e3, Day(10), loglog_quadratic)
end

@testitem "Geometric thresholds" begin
    # LopezAcosta2019 threshold functions
    stepwise = StepwiseLinearThresholdFunction(700, 0.5, 1.0)
    @test stepwise(500, 0.2) && stepwise(800, 0.2)
    @test !stepwise(500, 1.2) && !stepwise(800, 1.2)
    @test !stepwise(200, 0.7) && stepwise(800, 0.7)
end

@testitem "Filter function tests" begin
    using IceFloeTracker
    using DataFrames

    dataset = Watkins2026Dataset(; ref="v0.2")
    dataset = filter(c -> c.case_number == 9, dataset)
    cases = [x for x in dataset]
    if occursin("terra", name(cases[1]))
        terra = cases[1]
        aqua = cases[2]
    else
        aqua = cases[1]
        terra = cases[2]
    end

    if info(aqua)[:pass_time] < info(terra)[:pass_time]
        order = ["aqua", "terra"]
        labeled_images = [
            validated_labeled_floes(aqua).image_indexmap,
            validated_labeled_floes(terra).image_indexmap,
        ]
        passtimes = [info(aqua)[:pass_time], info(terra)[:pass_time]] # True time delta

    else
        order = ["terra", "aqua"]
        labeled_images = [
            validated_labeled_floes(terra).image_indexmap,
            validated_labeled_floes(aqua).image_indexmap,
        ]
        passtimes = [info(terra)[:pass_time], info(aqua)[:pass_time]] # True time delta
    end
    props = IceFloeTracker.regionprops_table.(labeled_images)

    # Adding floe masks: it may be that we need a step in the shape difference and the
    # psi-s curve test to check if the mask exists already, and only add it if it isn't there.
    add_floemasks!.(props, labeled_images)
    add_passtimes!.(props, passtimes)
    add_ψs!.(props)
    floe = props[1][1, :]
    candidates = props[1][2:end, :]

    n = nrow(candidates)
    candidates_after_map = map(DistanceThresholdFilter(), floe, candidates)
    n_after_map = nrow(candidates_after_map)
    candidates_after_filter = filter(DistanceThresholdFilter(), floe, candidates)
    n_after_filter = nrow(candidates_after_filter)
    # and that using it with just (floe, candidates) does the subsetting
    @test (n == n_after_map)
    @test (n >= n_after_filter)
    @test ("time_distance_test" ∉ names(candidates_after_filter))

    candidates = props[1][2:end, :]
    candidates = map(RelativeErrorThresholdFilter(; variable=:area), floe, candidates)
    candidates = map(
        RelativeErrorThresholdFilter(; variable=:convex_area), floe, candidates
    )
    # Check that variable names are being passed through correctly
    @test ("relative_error_area" ∈ names(candidates)) &&
        ("relative_error_convex_area" ∈ names(candidates))

    candidates = map(ShapeDifferenceThresholdFilter(), floe, candidates)
    @test "shape_difference_test" ∈ names(candidates)
    @test candidates[1, :shape_difference] == 268

    candidates = map(PsiSCorrelationThresholdFilter(), floe, candidates)
    @test "psi_s_correlation" ∈ names(candidates)
    @test candidates[1, :psi_s_correlation] == 0.914
end
