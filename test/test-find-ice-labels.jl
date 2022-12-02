@testset "Find_Ice_Labels" begin
    println("------------------------------------------------")
    println("------------ Create Ice Labels Test --------------")

    reflectance_image = float64.(load(reflectance_test_image_file)[test_region...])
    landmask = convert(BitMatrix, load(current_landmask_file))
    ice_labels_matlab = DelimitedFiles.readdlm(
        "$(test_data_dir)/ice_labels_matlab.csv", ','
    )
    ice_labels_matlab = vec(ice_labels_matlab)

    @time ice_labels_julia = IceFloeTracker.find_ice_labels(reflectance_image, landmask)

    DelimitedFiles.writedlm("ice_labels_julia.csv", ice_labels_julia, ',')

    @test ice_labels_matlab == ice_labels_julia

    @time ice_labels_ice_floe_region = IceFloeTracker.find_ice_labels(
        reflectance_image[ice_floe_test_region...], landmask[ice_floe_test_region...]
    )

    DelimitedFiles.writedlm("ice_labels_floe_region.csv", ice_labels_ice_floe_region, ',')
end
