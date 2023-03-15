# @testset "preprocess" begin
    landmask_imgs = deserialize(joinpath(test_data_dir,"pipeline/preprocess","landmasks.jls"))
    reflectance_img = IceFloeTracker.load(dir=".",fname=reflectance_test_image_file)[test_region...];
    truecolor_img = load(dir=".", fname=truecolor_test_image_file)[test_region...];

    segmented_expected = IceFloeTracker.load(raw"test\test_inputs\pipeline\preprocess\segmented_floes_expected.png") .> 0
    IceFloeTracker.Gray.(segmented_expected)

    segmented_floes = IceFloeTracker.preprocess(
        truecolor_img,
        reflectance_img,
        landmask_imgs)
# end