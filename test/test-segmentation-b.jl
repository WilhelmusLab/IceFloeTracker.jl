@testitem "Segmentation-B" begin
    using Images: @test_approx_eq_sigma_eps, load, float64, strel_diamond, closing, clamp01nan
    import StatsBase: mean
    include("config.jl")
    include("test_error_rate.jl")

    sharpened_image = float64.(load(sharpened_test_image_file))
    segmented_a_ice_mask = convert(BitMatrix, load(segmented_a_ice_mask_file))
    cloudmask = .!convert(BitMatrix, load(cloudmask_test_file))
    matlab_ice_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_c.png")
    )
    matlab_brightened = float64.(load("$(test_data_dir)/matlab_I.png"))


    # TBD: Update the segB function to match this exactly.
    not_ice_mask = deepcopy(sharpened_image)
    not_ice_mask[not_ice_mask .< 0.4] .= 0
    not_ice_bit = not_ice_mask .* 0.3
    not_ice_mask .= not_ice_bit .+ sharpened_image
    brightened_image = map(clamp01nan, not_ice_mask)


    # threshold_mask = sharpened_image .> 0.4
    # brightened_image = (sharpened_image .* 1.3) .* threshold_mask
    # clamp!(brightened_image, 0, 1)
    segB = LopezAcosta2019.segB_binarize(sharpened_image, brightened_image, cloudmask)
    
    @info "Merging segmentation results"
    segAB_intersect = closing(segmented_a_ice_mask, strel_diamond((5,5))) .* segB


    @persist brightened_image "./test_outputs/brightened_image.png" true
    @persist segAB_intersect "./test_outputs/segAB_intersect.png" true
  
    @test typeof(brightened_image) == typeof(matlab_brightened)
    
    # Temporary replacement of the test approx eq. The brightened image does look different than
    # the matlab image, and it's mainly in ocean areas where the clamped image is darker. 
    # It's worth checking how much differs with the choice of how to clamp / rescale intensity.
    @test (@test_approx_eq_sigma_eps brightened_image matlab_brightened [0, 0] 0.02) ===
        nothing # Images differ in the 

    @test mean(abs.(vec(Float64.(brightened_image) .- Float64.(matlab_brightened)))) < 0.04
    @test typeof(segAB_intersect) == typeof(matlab_ice_intersect)
    @test test_similarity(segAB_intersect, matlab_ice_intersect, 0.03) 
end
