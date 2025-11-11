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

@testitem "Extend regionprops threshold tests" begin
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
    floe = props[1][1, :]
    candidates = props[1][2:end, :]
    add_passtimes!(props, passtimes)

    shape_difference_test!(
        floe,
        candidates;
        threshold_function = PiecewiseLinearThresholdFunction(100, 700, 0.4, 0.2),
        threshold_column=:shape_difference_test,
        scale_by=:area,
        area_column=:area
    )
    @test "shape_difference_test" ∈ names(candidates)
    @test candidates[1, :shape_difference] == 268

    psi_s_correlation_test!(
        floe, 
        candidates;
        threshold_function=PiecewiseLinearThresholdFunction(100, 700, 0.8, 0.9),
        threshold_column=:psi_s_correlation_test,
        area_column=:area
    )
    @test "psi_s_correlation" ∈ names(candidates)
    @test candidates[1, :psi_s_correlation] == 0.914

    floe = props[1][1, :]
    candidates = props[2]

    time_distance_test!(
            floe,
            candidates;
            threshold_function=LopezAcostaTimeDistanceFunction(),
            threshold_column=:time_distance_test)
    @test "time_distance_test" ∈ names(candidates)
    @test candidates[1, :time_distance_test]

    relative_error_test!(
        floe,
        candidates;
        variable=:convex_area,
        threshold_column=:relative_error_convex_area,
        area_variable=:area,
        threshold_function=PiecewiseLinearThresholdFunction(100, 700, 0.4, 0.2)
        )

    @test "relative_error_convex_area" ∈ names(candidates)
end 
