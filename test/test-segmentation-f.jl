@testset "Segmentation-F" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-F Test --------------")

    ## Load inputs for comparison
    segmentation_B_not_ice_mask = float64.(load("$(test_data_dir)/matlab_not_ice_mask.png"))
    segmentation_C_ice_mask = load(segmented_c_test_file) .> 0.5
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    landmask = convert(BitMatrix, load(current_landmask_file))
    watershed_intersect = load(watershed_test_file) .> 0.5
    ice_labels =
        Int64.(
            vec(DelimitedFiles.readdlm("$(test_data_dir)/ice_labels_floe_region.csv", ','))
        )
    matlab_isolated_floes = convert(
        BitMatrix, load("$(test_data_dir)/matlab_isolated_floes.png")
    )

    ## Run function with Matlab inputs

    @time isolated_floes = IceFloeTracker.segmentation_F(
        segmentation_C_ice_mask[ice_floe_test_region...],
        segmentation_B_not_ice_mask[ice_floe_test_region...],
        watershed_intersect[ice_floe_test_region...],
        cloudmask[ice_floe_test_region...],
        landmask[ice_floe_test_region...],
        ice_labels,
    )

    isolated_floes_filename =
        "$(test_output_dir)/isolated_floes_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist isolated_floes isolated_floes_filename

    @test typeof(isolated_floes) == typeof(matlab_isolated_floes)
    @test test_similarity(
        isolated_floes, matlab_isolated_floes[ice_floe_test_region...], 0.07
    )
end
