@testset "tracker/matchcorr" begin
    path = joinpath(test_data_dir, "tracker")

    @testset "matchcorr" begin
        floes =
            deserialize.([
                joinpath(path, f) for f in ["f1.dat", "f2.dat", "f3.dat", "f4.dat"]
            ])
        mm, c = matchcorr(floes[1], floes[2], 400.0)
        @test isapprox(mm, 0.0; atol=0.05) && isapprox(c, 0.99; atol=0.05)
        @test all(isnan.(collect(matchcorr(floes[3], floes[4], 400.0))))
    end

    @testset "process_soit" begin
        soitdir = dirname("test/test_inputs/pipeline/tracker/passtimes_lat.csv")
        datetimes = process_soit(soitdir)
        @test names(datetimes) = ["sat", "pass_time"]
        @test nrow(datetimes) == 2
        @test typeof(datetimes.sat[1]) == String
        @test typeof(datetimes.pass_time[1]) == DateTime
    end

    @testset "tracker" begin
        # Set thresholds
        t1 = (dt=(30.0, 100.0, 1300.0), dist=(200, 250, 300))
        t2 = (
            area=1200,
            arearatio=0.28,
            majaxisratio=0.10,
            minaxisratio=0.12,
            convexarearatio=0.14,
        )
        t3 = (
            area=10_000,
            arearatio=0.18,
            majaxisratio=0.1,
            minaxisratio=0.15,
            convexarearatio=0.2,
        )
        condition_thresholds = (t1, t2, t3)
        mc_thresholds = (
            goodness=(area3=0.18, area2=0.236, corr=0.68), comp=(mxrot=10, sz=16)
        )
        dt = [15.0, 20.0]

        # Load data
        data = deserialize(joinpath(path, "tracker_test_data.dat"))

        # sort floe data by area
        for i in 1:3
            sort!(data.props[i], :area; rev=true)
        end

        pairs = IceFloeTracker.pairfloes(
            data.imgs, data.props, dt, condition_thresholds, mc_thresholds
        )
        @test length(pairs) == 2

        r = rand(1:(length(pairs[1].props1[:, :area])))
        f1 = pairs[1].props1[r, :]
        f2 = pairs[1].props2[r, :]
        ratios = pairs[1].ratios[r, :]
        @test sqrt(
            (f1.row_centroid - f2.row_centroid)^2 + (f1.col_centroid - f2.col_centroid)^2
        ) == pairs[1].dist[r]
        @test all(collect(ratios)[1:4] .< collect(t3)[2:end])
    end
end
