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