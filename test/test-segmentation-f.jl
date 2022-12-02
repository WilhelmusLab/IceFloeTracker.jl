@testset "Segmentation-F" begin
    println("------------------------------------------------")
    println("------------ Create Segmentation-F Test --------------")

    ## Load inputs for comparison
    segmentation_B_not_ice_mask = load("$(test_data_dir)/matlab_not_ice_mask.png") .> 0.5
    segmentation_C_ice_mask = convert(BitMatrix, load(segmented_c_test_file))
    cloudmask = convert(BitMatrix, load(cloudmask_test_file))
    watershed_intersect = load("$(test_data_dir)/matlab_watershed_intersect.png") .> 0.5
    ice_labels =
        Int64.(vec(DelimitedFiles.readdlm("$(test_data_dir)/ice_labels_matlab.csv", ',')))
    strel_disk_2 = readdlm(strel_file_2, ',', Bool)
    strel_disk_4 = readdlm(strel_file_4, ',', Bool)
    matlab_isolated_floes = convert(
        BitMatrix, load("$(test_data_dir)/matlab_isolated_floes.png")
    )

    ## Run function with Matlab inputs

    @time isolated_floes = IceFloeTracker.segmentation_F(
        segmentation_C_ice_mask,
        segmentation_B_not_ice_mask,
        watershed_intersect,
        cloudmask,
        ice_labels,
        strel_disk_2,
        strel_disk_4,
    )

    isolated_floes_filename =
        "$(test_output_dir)/isolated_floes_" *
        Dates.format(Dates.now(), "yyyy-mm-dd-HHMMSS") *
        ".png"
    IceFloeTracker.@persist isolated_floes isolated_floes_filename

    @test typeof(isolated_floes) == typeof(matlab_isolated_floes)
    @test test_similarity(isolated_floes, matlab_isolated_floes, 0.5)
end
