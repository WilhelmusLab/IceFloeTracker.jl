
@testitem "matchcorr" begin
    include("config.jl")
    path = joinpath(test_data_dir, "tracker")

    floes =
        deserialize.([
            joinpath(path, f) for f in ["f1.dat", "f2.dat", "f3.dat", "f4.dat"]
        ])
    mm, c, rot, corr_ci, sd_ci, rot_ci = matchcorr(floes[1], floes[2], 400.0)
    # Should use cases where the mm is known! I set the test to 0.6 just to see when things change.

    @test isapprox(mm, 0.6; atol=0.05) && isapprox(c, 0.99; atol=0.05)
   
    # floes 3 and 4 are too dissimilar in shape (corr = 0.91).
    # matchcorr returns NaN for the shape difference and rotation due to exiting after
    # the correlation test.
    mm, c, rot, corr_ci, sd_ci, rot_ci = matchcorr(floes[3], floes[4], 400.0)
    @test all([isnan(mm); c < 0.95; isnan(rot)])
end
