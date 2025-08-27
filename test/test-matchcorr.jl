
@testitem "matchcorr" setup = [Paths] begin

    path = joinpath(test_data_dir, "tracker")

    floes =
        deserialize.([joinpath(path, f) for f in ["f1.dat", "f2.dat", "f3.dat", "f4.dat"]])
    mm, c = matchcorr(floes[1], floes[2], 400.0)
    @test isapprox(mm, 0.0; atol=0.05) && isapprox(c, 0.99; atol=0.05)
    @test all(isnan.(collect(matchcorr(floes[3], floes[4], 400.0))))
end
