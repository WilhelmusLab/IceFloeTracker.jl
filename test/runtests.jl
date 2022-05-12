using IceFloeTracker
using ImageIO
using Images
using Test
using TiffImages

lm_image = load("/Users/tdivoll/IceFloeTracker.jl/data/input/landmask.tiff")

function create_landmask(landmask_image; num_pixels_dilate::Int=50, num_pixels_closing::Int=15)
    # Drop third dimension if it exists (test image had 3 dims: height x width x 1)
    landmask_image = dropdims(landmask_image, dims = 3)
    landmask_binary = Gray.(landmask_image) .== 0
    landmask_binary = LocalFilters.dilate(.!landmask_binary, num_pixels_dilate)
    landmask_binary = LocalFilters.closing(landmask_binary, num_pixels_closing)
    return landmask_binary
    # update to process inline
end

masked_image = create_landmask(lm_image, 50, 15)


@testset "IceFloeTracker.jl" begin
    @test created_landmask == matlab_landmask
    ## Next test will return a percent difference
    result = (@test_approx_eq_sigma_eps(created_landmask, matlab_landmask, 0, 0))
    #expected landmask #calc 95% similar
    result <= 0.95   
    
end
