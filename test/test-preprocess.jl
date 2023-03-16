@testset "preprocess" begin
    landmask_imgs = deserialize(joinpath(test_data_dir, "pipeline/preprocess", "landmasks.jls"))
    truecolor_img = IceFloeTracker.loadimg(dir=".", fname=truecolor_test_image_file)[test_region...]
    reflectance_img = IceFloeTracker.loadimg(dir=".", fname=reflectance_test_image_file)[test_region...]

    segmented_floes = IceFloeTracker.preprocess(
        truecolor_img,
        reflectance_img,
        landmask_imgs)

    segmented_floes_expected = load("$(test_data_dir)/matlab_BW7.png") .> 0.499

    @test test_similarity(segmented_floes[ice_floe_test_region...], segmented_floes_expected[ice_floe_test_region...], 0.071)
end