@testset "matchcorr" begin
    pth = joinpath(test_data_dir, "tracker")
    floes =
        deserialize.([joinpath(pth, f) for f in ["f1.dat", "f2.dat", "f3.dat", "f4.dat"]])
    mm, c = matchcorr(floes[1], floes[2], 400)
    @test isapprox(mm, 0.11; atol=0.05) && isapprox(c, 0.99; atol=0.05)
    @test all(isnan.(collect(matchcorr(floes[3], floes[4], 400))))
end
