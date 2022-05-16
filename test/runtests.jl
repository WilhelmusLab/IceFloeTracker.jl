using IceFloeTracker
using Images
using Test
using DelimitedFiles

@testset "IceFloeTracker.jl" begin
    lm_image = load("./data/landmask.tiff")
    matlab_landmask = load("./data/matlab_landmask.png")
    struct_elem = readdlm("../data/input/se.csv", ',', Bool)
    masked_image = IceFloeTracker.create_landmask(lm_image, struct_elem; num_pixels_closing=50)
    @test (@test_approx_eq_sigma_eps masked_image matlab_landmask [0,0] 0.001) === nothing 
end
