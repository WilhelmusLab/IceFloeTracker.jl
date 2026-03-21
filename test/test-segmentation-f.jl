@testitem "Segmentation-F" begin
    using DelimitedFiles: readdlm
    using Images: complement, float64, load

    include("config.jl")
    include("test_error_rate.jl")

    ## Load inputs for comparison
    brightened_gray = float64.(load("$(test_data_dir)/matlab_I.png"))
    ice_intersect = convert(BitMatrix, load(segmented_c_test_file))
    matlab_BW7 = load("$(test_data_dir)/matlab_BW7.png") .> 0.499
    fc_image = load(
        "$(test_data_dir)/beaufort-chukchi-seas_falsecolor.2020162.aqua.250m.tiff"
    )[test_region...]
    ## Load function arg files

    cloudmask = .!convert(BitMatrix, load(cloudmask_test_file))
    # convert ocean mask into land mask, so land=1
    landmask = convert(BitMatrix, load(current_landmask_file)[test_region...])
    watershed_intersect = load(watershed_test_file) .> 0.499
    # ice_labels = Int64.(vec(readdlm("$(test_data_dir)/ice_labels_floe_region.csv", ',')))

    # New method
    morphed_grayscale = LopezAcosta2019.reconstruct_and_mask(
        brightened_gray[ice_floe_test_region...],
        watershed_intersect[ice_floe_test_region...],
        ice_intersect[ice_floe_test_region...],
    )
    # kmeans binarization, again
    segF_binarized =
        kmeans_binarization(
            morphed_grayscale,
            apply_cloudmask(
                fc_image[ice_floe_test_region...], cloudmask[ice_floe_test_region...]
            );
            k=3,
            cluster_selection_algorithm=LopezAcosta2019.IceDetectionLopezAcosta2019(),
        ) .* .!watershed_intersect[ice_floe_test_region...]

    @info "Splitting floes"
    isolated_floes = LopezAcosta2019.morph_split_floes(
        segF_binarized, cloudmask[ice_floe_test_region...]
    )

    @test test_similarity(matlab_BW7[ice_floe_test_region...], isolated_floes, 0.013)

    @persist isolated_floes "./test_outputs/isolated_floes.png" true
    @persist matlab_BW7[ice_floe_test_region...] "./test_outputs/matlab_isolated_floes.png" true
end
