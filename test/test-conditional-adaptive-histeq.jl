@testset "Conditional adaptivehisteq" begin
    datadir = joinpath(@__DIR__, "test_inputs/")
    path_ref_image = joinpath(datadir, "NE_Greenland_reflectance.2020162.aqua.250m.tiff")
    ref_image = float64.(load(path_ref_image))
    dilated_landmask = BitMatrix(load(joinpath(datadir, "matlab_landmask.png")))

    redchannel_cahe = IceFloeTracker._get_red_channel_cloud_cae(
        false_color_image=ref_image,
        landmask=dilated_landmask,
        prelim_threshold=110.0,
        band_7_threshold=200.0,
        band_2_threshold=190.0,
    )

    @test sum(redchannel_cahe) == 1320925065
end
