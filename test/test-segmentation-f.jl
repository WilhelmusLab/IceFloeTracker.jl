@testitem "Segmentation-F" begin
    using DelimitedFiles
    using Images: complement
    
    include("config.jl")
    include("test_error_rate.jl")

    ## Load inputs for comparison
    segmentation_B_not_ice_mask = float64.(load("$(test_data_dir)/matlab_I.png"))
    segmentation_B_ice_intersect = convert(BitMatrix, load(segmented_c_test_file))
    matlab_BW7 = load("$(test_data_dir)/matlab_BW7.png") .> 0.499

    ## Load function arg files

    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    # convert ocean mask into land mask, so land=1
    landmask = .!convert(BitMatrix, load(current_landmask_file))
    watershed_intersect = load(watershed_test_file) .> 0.499
    ice_labels =
        Int64.(
            vec(DelimitedFiles.readdlm("$(test_data_dir)/ice_labels_floe_region.csv", ','))
        )

    ## Run function with Matlab inputs

    @time isolated_floes = IceFloeTracker.segmentation_F(
        segmentation_B_not_ice_mask[ice_floe_test_region...],
        segmentation_B_ice_intersect[ice_floe_test_region...],
        watershed_intersect[ice_floe_test_region...],
        ice_labels,
        .!cloudmask[ice_floe_test_region...],
        landmask[ice_floe_test_region...],
    )

    IceFloeTracker.@persist isolated_floes "./test_outputs/isolated_floes.png" true
    IceFloeTracker.@persist matlab_BW7[ice_floe_test_region...] "./test_outputs/matlab_isolated_floes.png" true

    @test test_similarity(matlab_BW7[ice_floe_test_region...], isolated_floes, 0.013)
end
