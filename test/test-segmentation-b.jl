@testitem "Segmentation-B" begin
    using Images: @test_approx_eq_sigma_eps, load, float64, strel_diamond, closing

    include("config.jl")
    include("test_error_rate.jl")

    sharpened_image = float64.(load(sharpened_test_image_file))
    segmented_a_ice_mask = convert(BitMatrix, load(segmented_a_ice_mask_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    matlab_ice_intersect = convert(
        BitMatrix, load("$(test_data_dir)/matlab_segmented_c.png")
    )
    matlab_brightened = float64.(load("$(test_data_dir)/matlab_I.png"))

    # @time segB = LopezAcosta2019.segmentation_B(
    #     sharpened_image, .!cloudmask, segmented_a_ice_mask
    # )

    threshold_mask = sharpened_image .> 0.4
    brightened_image = (sharpened_image .* 1.3) .* threshold_mask
    clamp!(brightened_image, 0, 1)
  
    segB = LopezAcosta2019.segB_binarize(sharpened_image, brightened_image, cloudmask)
    
    @info "Merging segmentation results"
    segAB_intersect = closing(segmented_a_ice_mask, strel_diamond((5,5))) .* segB



    @persist brightened_image "./test_outputs/brightened_image.png" true
    @persist segAB_intersect "./test_outputs/segAB_intersect.png" true
  
    @test typeof(brightened_image) == typeof(matlab_brightened)
    
    # no longer identical - why and where?
    @test (@test_approx_eq_sigma_eps brightened_image matlab_brightened [0, 0] 0.002) ===
        nothing

    @test typeof(segAB_intersect) == typeof(matlab_ice_intersect)
    @test test_similarity(segAB_intersect, matlab_ice_intersect, 0.005)
end
