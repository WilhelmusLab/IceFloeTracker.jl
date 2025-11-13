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
@test !stepwise(200 , 0.7) && stepwise(800, 0.7)

end

@testitem "Filter function tests" begin
    using IceFloeTracker
    using DataFrames

    data_loader = Watkins2025GitHub(; ref="25ba4d46814a5423b65ad675aaec05633d17a37e")
    dataset = data_loader(c-> c.case_number == 9)
    cases = [x for x in dataset.data]
    if occursin("terra", cases[1].name)
        terra = cases[1]
        aqua = cases[2]
    else
        aqua = cases[1]
        terra = cases[2]
    end

    if aqua.metadata[:pass_time] < terra.metadata[:pass_time]
        order = ["aqua", "terra"]
        labeled_images = [aqua.validated_labeled_floes.image_indexmap,
                    terra.validated_labeled_floes.image_indexmap]
        passtimes = [aqua.metadata[:pass_time], terra.metadata[:pass_time]] # True time delta

    else
        order = ["terra", "aqua"]
        labeled_images = [terra.validated_labeled_floes.image_indexmap,
                    aqua.validated_labeled_floes.image_indexmap]
        passtimes = [terra.metadata[:pass_time], aqua.metadata[:pass_time]] # True time delta

    end
    props = IceFloeTracker.regionprops_table.(labeled_images);

    # Adding floe masks: it may be that we need a step in the shape difference and the
    # psi-s curve test to check if the mask exists already, and only add it if it isn't there.
    greaterthan0(x) = x .> 0
    addfloemasks!(props, greaterthan0.(labeled_images))
    add_passtimes!(props, passtimes)
    floe = props[1][1, :]
    candidates = props[1][2:end, :]
    

    n = nrow(candidates)
    dt_test = DistanceThresholdFilter()
    dt_test(floe, candidates, Val(:raw))
    n2 = nrow(candidates)
    dt_test(floe, candidates)
    n3 = nrow(candidates)
    # Check that using the Val(:raw) forces the test to skip the subset step
    # and that using it with just (floe, candidates) does the subsetting
    @test (n == n2) && (n >= n3) && ("time_distance_test" ∉ names(candidates))

    candidates = props[1][2:end, :]
    re_test_area = RelativeErrorThresholdFilter(variable=:area)
    re_test_area(floe, candidates, Val(:raw))
    re_test_convex_area = RelativeErrorThresholdFilter(variable=:convex_area)
    re_test_convex_area(floe, candidates, Val(:raw))
    # Check that variable names are being passed through correctly
    @test ("relative_error_area" ∈ names(candidates)) && ("relative_error_convex_area" ∈ names(candidates))

    sd_test = ShapeDifferenceThresholdFilter()
    sd_test(floe, candidates, Val(:raw))
    @test "shape_difference_test" ∈ names(candidates)
    @test candidates[1, :shape_difference] == 268

    ps_test = PsiSCorrelationThresholdFunction()
    ps_test(floe, candidates, Val(:raw))
    @test "psi_s_correlation" ∈ names(candidates)
    @test candidates[1, :psi_s_correlation] == 0.914
end 
