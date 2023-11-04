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
        passtimes = deserialize(joinpath(path, "passtimes.dat"))
        latlondata = deserialize(joinpath(path, "latlondatatest.dat"))

        # Filtering out small floes. Algorithm performs poorly on small, amorphous floes as they seem to look similar (too `blobby`) to each other
        for (i, prop) in enumerate(data.props)
            data.props[i] = prop[prop[:, :area].>=350, :]
        end

        _pairs = IceFloeTracker.pairfloes(
            data.imgs, data.props, passtimes, dt, condition_thresholds, mc_thresholds
        )

        IceFloeTracker.addlatlon!(_pairs, latlondata)

        @test maximum(_pairs.ID) == 6
        @test names(_pairs) == ["ID",
            "passtime",
            "area",
            "convex_area",
            "major_axis_length",
            "minor_axis_length",
            "orientation",
            "perimeter",
            "area_under",
            "corr",
            "latitude",
            "longitude",
            "x",
            "y"]
        @test issorted(_pairs, :ID)
        @test all([issorted(grp, :passtime) for grp in groupby(_pairs, :ID)])
    end
end
